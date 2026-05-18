# Ce script permet de nettoyer les clusters avant la propagation d'annotations pour ne garder que ce qui est nécessaire et optimiser la propagation

import pandas as pd

# Choix du modèle considéré
model = "leiden"
param_vals = [0.1, 0.5, 1.0, 2.0, 5.0] # gamma

# model = "mcl"
# param_vals = [1.3, 1.5, 2.0, 2.5, 3.0] # r

# Chemins
global_path = "/tempory/21234701"
split_dir = f"{global_path}/filtered_clusters"

# BLAST & Embeddings
# Format : prot_id \t cluster_id
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
CLUSTERS_FILES = f"{global_path}/{model}/combined_graph_alpha_princ_comp_{model}"

# Fichier du graphe de score combiné correspondant
# Format attendu : prot_id1 \t prot_id2 \t score
SCORES_FILES =  f"{global_path}/combined_graph_alpha_"

FUNCTION_FILE = f"{global_path}/pfam_analysis/pfam_annotations.tsv"

def get_useful_clusters(clusters_path, annotations_dict):
    print(f"Analyse des clusters dans {clusters_path}")
    df = pd.read_csv(clusters_path, sep="\t")
    
    df['is_annotated'] = df['prot_id'].apply(lambda x: x in annotations_dict)
    stats = df.groupby('cluster_id')['is_annotated'].agg(['sum', 'count'])
    stats['nb_not_annotated'] = stats['count'] - stats['sum']
    
    # On garde les clusters qui ont du train et du test
    useful_ids = stats[(stats['sum'] > 0) & (stats['nb_not_annotated'] > 0)].index.tolist()
    
    print(f"Clusters totaux : {len(stats)} | Clusters utiles : {len(useful_ids)}")
    df_useful = df[df['cluster_id'].isin(useful_ids)]
    return df_useful, set(df_useful['prot_id'].unique())

def filter_scores_and_split(scores_path, useful_prots, prot_to_cluster, output_path):
    print(f"Filtrage du fichier de scores : {scores_path}")
    count_kept = 0
    with open(scores_path, 'r', encoding="utf-8", errors="replace") as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            parts = line.split()
            if len(parts) < 3:
                continue
            p1, p2, score = parts[0], parts[1], parts[2]
            
            if p1 in useful_prots and p2 in useful_prots:
                c1 = prot_to_cluster.get(p1)
                c2 = prot_to_cluster.get(p2)
                if c1 == c2 and c1 is not None:
                    f_out.write(f"{p1}\t{p2}\t{score}\t{c1}\n")
                    count_kept += 1
    print(f"Arêtes conservées : {count_kept}")

# Exécution du filtrage
print("Récupération des annotations")
annot_df = pd.read_csv(FUNCTION_FILE, sep="\t", usecols=["prot_id"]).dropna()
annot_set = set(annot_df["prot_id"].unique())

for a in ALPHA_LIST:
    scores_file = f"{SCORES_FILES}{a}.tsv"
    for param in param_vals:
        print(f"alpha = {a} et param = {param}")
        clusters_file = f"{CLUSTERS_FILES}{a}_{param}.tsv"
        
        useful_df, useful_prot_set = get_useful_clusters(clusters_file, annot_set)
        prot_to_cluster_map = dict(zip(useful_df.prot_id, useful_df.cluster_id))

        output_filtered_scores = f"{split_dir}/{model}_filtered_scores{a}_{param}.tsv"
        filter_scores_and_split(scores_file, useful_prot_set, prot_to_cluster_map, output_filtered_scores)
        
print(f"Terminé !")