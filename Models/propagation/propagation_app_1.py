import pandas as pd
from collections import Counter
import pickle
from sklearn.model_selection import train_test_split
import gc  # Pour forcer le nettoyage de la RAM entre les modèles

# Chemins
global_path = "/tempory/21234701"
split_dir = f"{global_path}/filtered_clusters"

model = "leiden"
param_vals = [1.0] # gamma

# model = "mcl"
# param_vals = [2.0] # r

ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]

FUNCTION_FILE = f"{global_path}/pfam_analysis/pfam_annotations.tsv"

with open("./pfam_compatible_labels.pkl", 'rb') as f:
    compatibilites = pickle.load(f)

def prepare_train_test(map_annotations, test_size=0.2):
    prot_ids = list(map_annotations.keys())
    train_ids, test_ids = train_test_split(prot_ids, test_size=test_size, random_state=42)
    
    train_annots = {pid: map_annotations[pid] for pid in train_ids}
    test_annots = {pid: map_annotations[pid] for pid in test_ids}
    
    with open(f"{global_path}/1_app/test_set_ground_truth.pkl", 'wb') as f:
        pickle.dump(test_annots, f)
        
    print(f"Train set : {len(train_annots)} ; Test set : {len(test_annots)}")
    return train_annots

def check_compatibility(list_a, list_b, compatibilites):
    for a in list_a:
        autorises = compatibilites.get(a, set())
        for b in list_b:
            if b not in autorises:
                return False
    return True

def resolve_internal_conflicts(top_labels, counts, compatibilites):
    clean_labels = set(top_labels)
    to_remove = set()
    list_labels = list(clean_labels)
    for i in range(len(list_labels)):
        for j in range(i + 1, len(list_labels)):
            a, b = list_labels[i], list_labels[j]
            if b not in compatibilites.get(a, set()):
                if counts[a] > counts[b]:
                    to_remove.add(b)
                elif counts[b] > counts[a]:
                    to_remove.add(a)
                else:
                    to_remove.add(a)
                    to_remove.add(b)
    return list(clean_labels - to_remove)

def cluster_propagate_step(cluster_prots, cluster_rows, annotations, compatibilites, threshold=0.7):
    cluster_labels = []
    for p in cluster_prots:
        if p in annotations:
            cluster_labels.extend(annotations[p])
    
    if not cluster_labels:
        return {}

    counts = Counter(cluster_labels)
    total_labels = len(cluster_labels)
    
    raw_top_labels = [label for label, count in counts.items() if count / total_labels > 0.3]
    top_labels = resolve_internal_conflicts(raw_top_labels, counts, compatibilites)
    
    if not top_labels:
        return {}

    updates = {}
    # cluster_rows est une liste de tuples (id1, id2, score)
    for p1, p2, score in cluster_rows:
        if score < threshold:
            continue
        
        if p1 in annotations and p2 not in annotations:
            intersection = set(annotations[p1]).intersection(top_labels)
            if intersection:
                updates[p2] = list(intersection)
                
        elif p2 in annotations and p1 not in annotations:
            intersection = set(annotations[p2]).intersection(top_labels)
            if intersection:
                updates[p1] = list(intersection)
                
        elif p1 in annotations and p2 in annotations:
            if check_compatibility(top_labels, annotations[p1], compatibilites):
                union_p1 = list(set(annotations[p1]).union(top_labels))
                if len(union_p1) > len(annotations[p1]):
                    updates[p1] = union_p1
                    
            if check_compatibility(top_labels, annotations[p2], compatibilites):
                union_p2 = list(set(annotations[p2]).union(top_labels))
                if len(union_p2) > len(annotations[p2]):
                    updates[p2] = union_p2
                
    return updates

def run_single_graph_propagation(filtered_scores_path, train_annots, compatibilites, max_iter=10):
    annots_current = train_annots.copy()
    
    print("   Chargement du fichier de scores en mémoire...")
    columns = ["id1", "id2", "score", "cluster_id1"]
    df_filtered = pd.read_csv(filtered_scores_path, sep="\t", header=None, names=columns)
    
    # Pour chaque cluster, on garde le set de ses protéines et la liste de ses arêtes
    print("Pré-indexation des clusters")
    cluster_data = {}
    for cluster_id, group in df_filtered.groupby('cluster_id1'):
        prots = set(group['id1']).union(set(group['id2']))
        
        # Convertion du sous-dataframe en tableau numpy
        rows = list(zip(group['id1'], group['id2'], group['score']))
        cluster_data[cluster_id] = (prots, rows)
    
    # Libération de la mémoire du DataFrame
    del df_filtered
    gc.collect()
    
    # Boucle de propagation
    for i in range(max_iter):
        count_start = len(annots_current)
        
        for cluster_id, (prots, rows) in cluster_data.items():
            updates = cluster_propagate_step(prots, rows, annots_current, compatibilites, threshold=0.7)
            annots_current.update(updates)
            
        diff = len(annots_current) - count_start
        print(f"   Itération {i+1} -> Nouvelles protéines : {diff} | Total : {len(annots_current)}")
        
        if diff == 0:
            print("   Convergence atteinte.")
            break
            
    return annots_current

# Chargement initial des annotations
print("Chargement des annotations globales")
df_annots = pd.read_csv(FUNCTION_FILE, sep="\t", encoding="utf-8").dropna(subset=["prot_id", "annot_id"])
map_annotations = df_annots.groupby("prot_id")["annot_id"].apply(list).to_dict()
train_annots = prepare_train_test(map_annotations)

for a in ALPHA_LIST:
    for param in param_vals:
        filtered_scores_file = f"{split_dir}/{model}_filtered_scores{a}_{param}.tsv"
            
        print(f"Lancement : alpha = {a} | param = {param}")
        final_annotations = run_single_graph_propagation(filtered_scores_file, train_annots, compatibilites)

        # Sauvegarde
        output_pkl = f"{global_path}/1_app/{model}_propagation_results{a}_{param}.pkl"
        with open(output_pkl, 'wb') as f:
            pickle.dump(final_annotations, f)
            
        del final_annotations
        gc.collect()
            
print(f"Terminé !")