# Affiche les scores des métriques d'évaluation pour les clusters passés en paramètre
# Les métriques sont implémentées pour des points dans sklearn, mais pas pour des graphes.
# Nous allons donc devoir les implémenter à la main...

BLAST_COMP_FILE = "/tempory/21234701/blast_components.tsv"                    # Fichier prévu pour contenir les composantes de BLAST
PRIN_COMP_CLUSTERS = "/tempory/21234701/blast_princ_comp_leiden_clusters.tsv" # Clusters de la composante principale
OUTPUT_CSV = "/tempory/21234701/blast_cluster_validation_scores.csv"          # Fichier de sortie, csv avec les scores

import pandas as pd
import numpy as np

# Chargement des clusters
print("Chargement des clusters")
node_to_cluster = {}

# Chargement des résultats de Leiden (Composante 0)
leiden_df = pd.read_csv(PRIN_COMP_CLUSTERS, sep='\t')
for _, row in leiden_df.iterrows():
    node_to_cluster[row['prot_id']] = f"C0_{row['cluster_id']}"

# Initialisation des variables permettant de calculer les métriques
""" Métriques considérées
- Co-distance : somme des distances intra-cluster
- Index de Dunn : distance minimale entre deux clusters / diamètre interne du cluster
- Coupe : somme des poids des arêtes reliant un cluster au reste de la partition
- Coupe normalisée : divise la valeur de la coupe par le volume total du cluster (somme de tous les poids de ses arêtes)

Le coefficient de silhouette a été exclu en raison de sa trop grande complexité.

Chaque métrique est calculée pour un cluster.
Afin d'obtenir les résultats pour l'ensemble de la partition, il suffit de faire la somme ou la moyenne selon la métrique considérée.
"""

# On stocke les sommes pour calculer Co-distance, Coupe (cut), Vol
intra_weights = {}    # Somme w_ij internes
inter_weights = {}    # Somme w_ij externes
cluster_volume = {}   # Somme totale w par cluster
max_intra_dist = {}   # Pour le diamètre (Dunn)
min_inter_dist = {}   # Pour la distance inter (Dunn)

# Statistiques des clusters
nb_edges = {}         # Nombre d'arêtes
nb_nodes = {}         # Nombre de nœuds (protéines)

print("Parcours du graphe pour le calcul des métriques")
with open(BLAST_COMP_FILE, "r") as f:
    next(f) # Header
    for line in f:
        parts = line.split()
        u, v, w = parts[0], parts[1], float(parts[2])
        comp_id = parts[3]
        
        # Déterminer le cluster
        # Soit la protéine est dans un cluster de la composante principale, soit dans une autre composante
        c_u = node_to_cluster.get(u, f"Comp_{comp_id}")
        c_v = node_to_cluster.get(v, f"Comp_{comp_id}")
        
        # Conversion similarité -> distance
        d = 1.0 - min(w, 1.0) # w entre 0 et 1.0
        
        # Initialisation à 0 des métriques des clusters
        for c in [c_u, c_v]:
            if c not in cluster_volume:
                cluster_volume[c] = 0.0
                intra_weights[c] = 0.0
                inter_weights[c] = 0.0
                max_intra_dist[c] = 0.0
                nb_edges[c] = 0
                nb_nodes[c] = 0
        
        cluster_volume[c_u] += w
        cluster_volume[c_v] += w
        
        # Si on a une arête dans un même cluster
        if c_u == c_v:
            nb_edges[c_u] += 1
            nb_nodes[c_u] += 2
            
            # Intra-cluster
            intra_weights[c_u] += w
            # Maj diamètre
            if d > max_intra_dist[c_u]:
                max_intra_dist[c_u] = d
        else:
            nb_nodes[c_u] += 1
            nb_nodes[c_u] += 1
            
            # Inter-cluster - Coupe
            inter_weights[c_u] += w
            inter_weights[c_v] += w
            
            pair = tuple(sorted((c_u, c_v)))
            if pair not in min_inter_dist or d < min_inter_dist[pair]:
                min_inter_dist[pair] = d
                
# Calcul des scores par cluster
results = []

# Distance minimale globale entre clusters pour Dunn
global_min_inter = min(min_inter_dist.values()) if min_inter_dist else 1.0

print("Finalisation des calculs")
for c in cluster_volume:
    # Coupe et NCut
    cut = inter_weights[c]
    vol = cluster_volume[c]
    ncut = cut / vol if vol > 0 else 0
    
    # Dunn local (Estimation : min_inter / diamètre_propre)
    diam = max_intra_dist[c]
    dunn = global_min_inter / diam if diam > 0 else 0
    
    results.append({"cluster_id": c, "co_distance": intra_weights[c], "dunn_index": dunn,
                    "cut": cut, "ncut": ncut, "nb_nodes": nb_nodes[c], "nb_edges": nb_edges[c]})

df_final = pd.DataFrame(results)
df_final.to_csv(OUTPUT_CSV, index=False)
print(f"Terminé !")