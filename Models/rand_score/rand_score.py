import pandas as pd
from sklearn.metrics import adjusted_rand_score
import csv

# Choix du modèle considéré
#model = "leiden"
#param_vals = [0.1, 0.5, 1.0, 2.0, 5.0] # gamma
#param_name = "gamma"

model = "mcl"
param_vals = [1.3, 1.5, 2.0, 2.5, 3.0] # mcl
param_name = "inflation"

# Chemins
global_path = "/tempory/21234701"

# BLAST
CLUSTERS_BLAST = f"{global_path}/{model}/reduced_blast_princ_comp_{model}_clusters" # Clusters de sortie

# BLAST & Embeddings
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]

EMB_BLAST_CLUSTERS = f"{global_path}/{model}/combined_graph_alpha_princ_comp_{model}"
emb_blast_files = [f"{EMB_BLAST_CLUSTERS}{a}" for a in ALPHA_LIST] 

# Analyse des différences entre les clusters de BLAST et BLAST & Embeddings pour les mêmes valeurs de param
rand_scores = {} # param: [<randscores>]

for param in param_vals:
    rand_scores[param] = []
    
    # Compare BLAST avec Emb & BLAST pour tous les alphas
    for emb_blast_file in emb_blast_files:
        # Récupérer les dataframes
        emb_blast_df = pd.read_csv(emb_blast_file + f"_{param}.tsv", sep='\t')
        emb_blast_df = emb_blast_df.rename(columns={'cluster_id': 'cluster_id_emb_blast'})
        blast_df = pd.read_csv(CLUSTERS_BLAST + f"_{param}.tsv", sep='\t')
        blast_df = blast_df.rename(columns={'cluster_id': 'cluster_id_blast'})
        
        # Intersection sur les IDs
        merged_df = pd.merge(emb_blast_df, blast_df, how='inner', on=['prot_id'])
        
        rand_scores[param].append(adjusted_rand_score(merged_df["cluster_id_emb_blast"].to_numpy(), merged_df["cluster_id_blast"].to_numpy()))

# Sauvegarder les rand scores dans un fichier CSV
with open(f"./rand_score_{model}.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    # Header avec les params
    writer.writerow(['alpha'] + [f"{param_name}_{param}" for param in param_vals])
    # Lignes avec les alphas et leurs scores
    for i, alpha in enumerate(ALPHA_LIST):
        writer.writerow([alpha] + [rand_scores[param][i] for param in param_vals])

print("terminé !")