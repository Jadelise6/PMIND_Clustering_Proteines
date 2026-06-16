# Script effectuant la boucle principale de la deuxième approche

import pandas as pd
from collections import Counter
import pickle
from sklearn.model_selection import train_test_split

# Chemins
global_path = "/tempory/21234701"

split_dir = f"{global_path}/filtered_clusters"
embed_clusters = f"{split_dir}/embed_filtered_scores.tsv"
blast_clusters = f"{split_dir}/blat_filtered_scores.tsv"

# Annotations compatibles - format : {id_annot: set(annotations_compatibles)}
with open("../pfam_compatible_labels.pkl", 'rb') as f:
    compatibilites = pickle.load(f)
    
FUNCTION_FILE = f"{global_path}/pfam_analysis/pfam_annotations.tsv"     

def prepare_train_test(map_annotations, test_size=0.2):
    """ Sépare les annotations en train (visibles) et test (cachées) """
    prot_ids = list(map_annotations.keys())
    train_ids, test_ids = train_test_split(prot_ids, test_size=test_size, random_state=42)
    
    train_annots = {pid: map_annotations[pid] for pid in train_ids}
    test_annots = {pid: map_annotations[pid] for pid in test_ids}
    
    with open(f"{global_path}/2_app/test_set_ground_truth.pkl", 'wb') as f:
        pickle.dump(test_annots, f)
        
    print(f"Train set : {len(train_annots)} ; Test set : {len(test_annots)}")
    return train_annots

def read_annotations_with_fallback(path):
    try:
        return pd.read_csv(path, sep="\t", encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep="\t", encoding="latin-1")
    
def get_annotations(annot_file):
    print("Récupération des annotations")
    df = read_annotations_with_fallback(annot_file)
    map_annotations = df[["prot_id", "annot_id"]].dropna(subset=["prot_id", "annot_id"]).groupby("prot_id")["annot_id"].apply(list).to_dict()
    return map_annotations

def check_compatibility(list_a, list_b, compatibilites):
    """ Vérifie si toutes les étiquettes de list_a sont compatibles avec celles de list_b """
    for a in list_a:
        autorises = compatibilites.get(a, set())
        # Si une seule annotation de b n'est pas dans les compatibles de a, il y a conflit, on ne propage pas
        for b in list_b:
            if b not in autorises:
                return False
    return True

def resolve_internal_conflicts(top_labels, counts, compatibilites):
    """ Résout les conflits internes au sein des étiquettes majoritaires du cluster """
    clean_labels = set(top_labels)
    to_remove = set()
    
    # On cherche les paires incompatibles dans les étiquettes retenues
    list_labels = list(clean_labels)
    for i in range(len(list_labels)):
        for j in range(i + 1, len(list_labels)):
            a, b = list_labels[i], list_labels[j]
            if b not in compatibilites.get(a, set()):
                # Conflit, l'annotation avec la proportion la plus grande l'emporte
                if counts[a] > counts[b]:
                    to_remove.add(b)
                elif counts[b] > counts[a]:
                    to_remove.add(a)
                else:
                    # Égalité parfaite : par prudence, on supprime les deux
                    to_remove.add(a)
                    to_remove.add(b)
                    
    return list(clean_labels - to_remove)

def propagate_step(clusters_file, annotations, compatibilites, threshold=0.7):
    """ Une itération de propagation sur un graphe (BLAST ou Embeddings) """
    new_annotations = annotations.copy()
    columns = ["id1", "id2", "score", "cluster_id1"]
    df_filtered = pd.read_csv(clusters_file, sep="\t", encoding="utf-8", header=None, names=columns)
    
    for cluster_id, group in df_filtered.groupby('cluster_id1'):
        updates = cluster_propagate_step(group, new_annotations, compatibilites, threshold)
        new_annotations.update(updates)
        
    return new_annotations
    
def cluster_propagate_step(cluster_group, annotations, compatibilites, threshold=0.7):
    prots_in_cluster = set(cluster_group['id1']).union(set(cluster_group['id2']))
    
    cluster_labels = []
    for p in prots_in_cluster:
        if p in annotations:
            cluster_labels.extend(annotations[p])
    
    if not cluster_labels:
        return {} 

    counts = Counter(cluster_labels)
    total_labels = len(cluster_labels)
    
    # Seuil à 30% d'occurrence dans le cluster
    raw_top_labels = [label for label, count in counts.items() if count / total_labels > 0.3]
    
    # Nettoyage des conflits internes au profil du cluster avant propagation
    top_labels = resolve_internal_conflicts(raw_top_labels, counts, compatibilites)
    if not top_labels:
        return {}

    updates = {}
    for _, row in cluster_group.iterrows():
        if row['score'] < threshold:
            continue
        
        p1, p2 = row['id1'], row['id2']
        
        # Cas 1 : p1 est annoté, p2 ne l'est pas
        if p1 in annotations and p2 not in annotations:
            intersection = set(annotations[p1]).intersection(top_labels)
            if intersection:
                updates[p2] = list(intersection)
                
        # Cas 2 : p2 est annoté, p1 ne l'est pas
        elif p2 in annotations and p1 not in annotations:
            intersection = set(annotations[p2]).intersection(top_labels)
            if intersection:
                updates[p1] = list(intersection)
                
        # Cas 3 : Les deux sont déjà annotés
        elif p1 in annotations and p2 in annotations:
            # On vérifie si les étiquettes majoritaires du cluster sont compatibles avec p1 et p2
            if check_compatibility(top_labels, annotations[p1], compatibilites):
                # Union sans doublons
                union_p1 = list(set(annotations[p1]).union(top_labels))
                if len(union_p1) > len(annotations[p1]):
                    updates[p1] = union_p1
                    
            if check_compatibility(top_labels, annotations[p2], compatibilites):
                union_p2 = list(set(annotations[p2]).union(top_labels))
                if len(union_p2) > len(annotations[p2]):
                    updates[p2] = union_p2
                
    return updates

def update_annotations(annots_source, annots_cible, compatibilites):
    """ Transfert sécurisé des annotations d'une partition à l'autre (b et d) """
    for prot_id, labels_source in annots_source.items():
        # Si la protéine n'est pas du tout annotée dans la cible
        if prot_id not in annots_cible:
            annots_cible[prot_id] = labels_source
        else:
            # Si elle est déjà annotée, on applique la tolérance zéro conflit
            labels_cible = annots_cible[prot_id]
            if check_compatibility(labels_source, labels_cible, compatibilites):
                # Tout est compatible -> Union
                annots_cible[prot_id] = list(set(labels_cible).union(labels_source))
            # Sinon -> paas de propagation, on conserve l'état d'origine de la cible
                
    return annots_cible

def main_co_propagation(blast_clusters, embed_clusters, train_annots, compatibilites, max_iter=10):
    annots_blast = train_annots.copy()
    annots_embed = train_annots.copy()
        
    for i in range(max_iter):
        print(f"Itération {i+1}")
        
        count_blast_start = len(annots_blast)
        count_embed_start = len(annots_embed)
        
        print("Propagation BLAST")
        annots_blast = propagate_step(blast_clusters, annots_blast, compatibilites, threshold=0.7)
        
        print("Reporter BLAST sur Embeddings")
        annots_embed = update_annotations(annots_blast, annots_embed, compatibilites)
        
        print("Propagation Embeddings")
        annots_embed = propagate_step(embed_clusters, annots_embed, compatibilites, threshold=0.7)
        
        print("Reporter Embeddings sur BLAST")
        annots_blast = update_annotations(annots_embed, annots_blast, compatibilites)
        
        diff_blast = len(annots_blast) - count_blast_start
        diff_embed = len(annots_embed) - count_embed_start
        
        print(f"Gain BLAST: {diff_blast} ; Gain Embeddings: {diff_embed}")
        print(f"Total annoté dans BLAST : {len(annots_blast)}")
        
        if diff_blast == 0 and diff_embed == 0:
            print("Convergence atteinte : aucune nouvelle modification détectée.")
            break
                    
    return annots_blast, annots_embed

# Exécution
map_annotations = get_annotations(FUNCTION_FILE)
train_annots = prepare_train_test(map_annotations)
    
final_blast, final_embed = main_co_propagation(blast_clusters, embed_clusters, train_annots, compatibilites)

results = {'blast_results': final_blast, 'embed_results': final_embed}
with open(f"{global_path}/2_app/propagation_results2.pkl", 'wb') as f:
    pickle.dump(results, f)

print("Terminé avec succès !")