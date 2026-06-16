# Script pour séparer des blocs de clusters dans différents fichiers afin d'alléger le poids en RAM pour la deuxième approche

# Points importants du script :
# - découpage des fichiers de clusters en petits fichiers
# - enlever les clusters sans aucune annotation et avec que des annotations (0% et 100%) pour libérer de la place
# - retirer les protéine de ces clusters sans annotation de nos fichiers de score

import pandas as pd

# Chemins
global_path = "/tempory/21234701"

# BLAST et embeddings avec Leiden et gamma = 1.0 (meilleur modèle et meilleur paramètre)
# Format : prot_id cluster_id
blast_clusters_file = f"{global_path}/leiden/combined_graph_alpha_princ_comp_leiden1.0_1.0.tsv"
embed_clusters_file = f"{global_path}/leiden/combined_graph_alpha_princ_comp_leiden0.0_1.0.tsv"

# Fichiers des graphes
# Format : prot_id1 prot_id2 score (pas d'entête)
blast_scores_file = f"{global_path}/combined_graph_alpha_1.0.tsv"
embed_scores_file = f"{global_path}/combined_graph_alpha_0.0.tsv"

# Répertoire pour stocker les clusters et scores filtrés
split_dir = f"{global_path}/filtered_clusters"

FUNCTION_FILE = f"{global_path}/pfam_analysis/pfam_annotations.tsv"   # Fichier des annotations, format : prot_id source  annot_id  desc

def get_useful_clusters(clusters_path, annotations_dict):
    """ Identifie les clusters qui ont au moins une annotation et au moins une protéine sans annotation
    """
    print(f"Analyse des clusters dans {clusters_path}")
    df = pd.read_csv(clusters_path, sep="\t")
    
    # Marquer si chaque protéine est annotée
    df['is_annotated'] = df['prot_id'].apply(lambda x: x in annotations_dict)
    
    # Grouper par cluster pour compter les prot annotées et le nombre de prot total
    stats = df.groupby('cluster_id')['is_annotated'].agg(['sum', 'count'])
    stats['nb_not_annotated'] = stats['count'] - stats['sum']
    
    # On veut sum > 0 et nb_not_annotated > 0
    useful_ids = stats[(stats['sum'] > 0) & (stats['nb_not_annotated'] > 0)].index.tolist()
    
    print(f"Clusters totaux : {len(stats)}")
    print(f"Clusters utiles : {len(useful_ids)}")
    
    # On retourne le DataFrale filtré et le set des prot_id à conserver pour les scores
    df_useful = df[df['cluster_id'].isin(useful_ids)]
    
    return df_useful, set(df_useful['prot_id'].unique())

def filter_scores_and_split(scores_path, useful_prots, prot_to_cluster, file_prefix):
    """ Parcourt le fichier de scores, ne garde que les arêtes entre protéines utiles du même cluster, et écrit par blocs
    """
    print(f"Filtrage du fichier de scores : {scores_path}")
    output_scores = f"{file_prefix}filtered_scores.tsv"
    
    count_kept = 0
    with open(scores_path, 'r', encoding="utf-8", errors="replace") as f_in, open(output_scores, 'w') as f_out:
        for line in f_in:
            parts = line.split()
            if len(parts) < 3:
                continue
            
            p1, p2, score = parts[0], parts[1], parts[2]
            
            # Les deux protéines sont dans la liste "utile"
            if p1 in useful_prots and p2 in useful_prots:
                # Elles sont dans le même cluster
                c1 = prot_to_cluster.get(p1)
                c2 = prot_to_cluster.get(p2)
                
                if c1 == c2 and c1 is not None:
                    f_out.write(f"{p1}\t{p2}\t{score}\t{c1}\n")
                    count_kept += 1
                    
    print(f"Arêtes conservées : {count_kept}")

# Charger les annotations en dictionnaire - on n'a besoin que des clés pour savoir qui est annoté
print("Récupération des annotations")
annot_df = pd.read_csv(FUNCTION_FILE, sep="\t", usecols=["prot_id"]).dropna()
annot_set = set(annot_df["prot_id"].unique())

# Traitement BLAST
blast_useful_df, blast_prot_set = get_useful_clusters(blast_clusters_file, annot_set)
blast_map = dict(zip(blast_useful_df.prot_id, blast_useful_df.cluster_id))
filter_scores_and_split(blast_scores_file, blast_prot_set, blast_map, f"{split_dir}/blat_")

# Traitement Embeddings
embed_useful_df, embed_prot_set = get_useful_clusters(embed_clusters_file, annot_set)
embed_map = dict(zip(embed_useful_df.prot_id, embed_useful_df.cluster_id))
filter_scores_and_split(embed_scores_file, embed_prot_set, embed_map, f"{split_dir}/embed_")

print("Terminé !")