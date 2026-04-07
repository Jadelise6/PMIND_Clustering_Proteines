import pandas as pd
from sklearn.metrics import adjusted_rand_score
import csv

# BLAST
CLUSTERS_BLAST = "/tempory/21234701/reduced_blast_princ_comp_leiden_clusters" # Clusters de sortie

# BLAST & Embeddings
ALPHA_LIST = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.93, 0.96, 0.99, 1.0]

EMB_BLAST_CLUSTERS = "/tempory/21234701/combined_graph_alpha_princ_comp_leiden"
emb_blast_files = [f"{EMB_BLAST_CLUSTERS}{a}" for a in ALPHA_LIST] 

# Analyse des différences entre les clusters de BLAST et BLAST & Embeddings pour les mêmes valeurs de gamma
gamma_vals = [0.1, 0.5, 1.0, 2.0, 5.0]
rand_scores = {} # gamma: [<randscores>]

for gamma in gamma_vals:
    rand_scores[gamma] = []
    
    # Compare BLAST avec Emb & BLAST pour tous les alphas
    for emb_blast_file in emb_blast_files:
        # Récupérer les dataframes
        emb_blast_df = pd.read_csv(emb_blast_file + f"_{gamma}.tsv", sep='\t')
        emb_blast_df = emb_blast_df.rename(columns={'cluster_id': 'cluster_id_emb_blast'})
        blast_df = pd.read_csv(CLUSTERS_BLAST + f"_{gamma}.tsv", sep='\t')
        blast_df = blast_df.rename(columns={'cluster_id': 'cluster_id_blast'})
        
        # Intersection sur les IDs
        merged_df = pd.merge(emb_blast_df, blast_df, how='inner', on=['prot_id'])
        
        rand_scores[gamma].append(adjusted_rand_score(merged_df["cluster_id_emb_blast"].to_numpy(), merged_df["cluster_id_blast"].to_numpy()))

# Sauvegarder les rand scores dans un fichier CSV
with open("./rand_score.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    # Header avec les gammas
    writer.writerow(['alpha'] + [f"gamma_{gamma}" for gamma in gamma_vals])
    # Lignes avec les alphas et leurs scores
    for i, alpha in enumerate(ALPHA_LIST):
        writer.writerow([alpha] + [rand_scores[gamma][i] for gamma in gamma_vals])
