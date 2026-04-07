# Affiche les scores des métriques d'évaluation pour les clusters passés en paramètre
# Les métriques sont implémentées pour des points dans sklearn, mais pas pour des graphes.
# Nous allons donc devoir les implémenter à la main...
import pandas as pd
import numpy as np

# BLAST
BLAST_COMP_FILE = "/tempory/21234701/reduced_blast_components.tsv"                # Fichier prévu pour contenir les composantes de BLAST
PRIN_COMP_CLUSTERS = "/tempory/21234701/reduced_blast_princ_comp_leiden_clusters" # Clusters de la composante principale
OUTPUT_CSV = "/tempory/21234701/reduced_blast_cluster_validation_scores"          # Fichier de sortie, csv avec les scores

# BLAST & Embeddings
ALPHA_LIST = [0.2, 0.4, 0.5, 0.6, 0.8]
COMP_FILES = "/tempory/21234701/components_combined_graph_alpha_"    # Racine fichier prévu pour contenir les composantes des fichiers
com_files = [f"{COMP_FILES}{a}.tsv" for a in ALPHA_LIST]

OUTPUT_CLUSTERS = "/tempory/21234701/combined_graph_alpha_princ_comp_leiden"
prin_comp_files = [f"{OUTPUT_CLUSTERS}{a}" for a in ALPHA_LIST]

ALPHA_OUTPUT_CSV = "/tempory/21234701/combined_graph_alpha_cluster_validation_scores"     # Fichier de sortie, csv avec les scores
output_csvs = [f"{ALPHA_OUTPUT_CSV}{a}" for a in ALPHA_LIST]

def traitement_graphe(com_file, prin_comp_file, output_csv):
    gamma_vals = [0.1, 0.5, 1.0, 2.0, 5.0]
    for gamma in gamma_vals:
        print(f"gamma = {gamma}")
        
        # Chargement des clusters
        print("Chargement des clusters")
        node_to_cluster = {}

        # Chargement des résultats de Leiden (Composante principale de BLAST)
        leiden_df = pd.read_csv(prin_comp_file + f"_{gamma}.tsv", sep='\t')
        for _, row in leiden_df.iterrows():
            node_to_cluster[row['prot_id']] = f"{row['cluster_id']}"

        # Initialisation des variables permettant de calculer les métriques
        """ Métriques considérées
        - Co-distance : somme des distances intra-cluster
        - Coupe normalisée : divise la valeur de la coupe par le volume total du cluster (somme de tous les poids de ses arêtes)

        Le coefficient de silhouette a été exclu en raison de sa trop grande complexité.

        Chaque métrique est calculée pour un cluster.
        Afin d'obtenir les résultats pour l'ensemble de la partition, il suffit de faire la somme ou la moyenne selon la métrique considérée.
        """

        # On stocke les sommes pour calculer Co-distance, Vol
        intra_weights = {}    # Somme w_ij internes
        inter_weights = {}    # Somme w_ij externes
        cluster_volume = {}   # Somme totale w par cluster

        # Statistiques des clusters
        nb_edges = {}         # Nombre d'arêtes
        nb_nodes = {}         # Nombre de nœuds (protéines)

        print("Parcours du graphe pour le calcul des métriques")
        with open(com_file, "r") as f:
            next(f) # Header
            for line in f:
                parts = line.split()
                u, v, w = parts[0], parts[1], float(parts[2])
                comp_id = parts[3]
                
                # Déterminer le cluster
                # Soit la protéine est dans un cluster de la composante principale, soit dans une autre composante
                c_u = node_to_cluster.get(u, "")
                c_v = node_to_cluster.get(v, "")
                
                # On ne garde que les protéines qui sont dans la composante principale
                if c_u == "" or c_v == "":
                    continue
                
                # Conversion similarité -> distance
                d = 1.0 - min(w, 1.0) # w entre 0 et 1.0
                
                # Initialisation à 0 des métriques des clusters
                for c in [c_u, c_v]:
                    if c not in cluster_volume:
                        cluster_volume[c] = 0.0
                        intra_weights[c] = 0.0
                        inter_weights[c] = 0.0
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

                else:
                    nb_nodes[c_u] += 1
                    nb_nodes[c_u] += 1
                    
                    # Inter-cluster - Coupe
                    inter_weights[c_u] += w
                    inter_weights[c_v] += w    
                        
        # Calcul des scores par cluster
        results = []

        print("Finalisation des calculs")
        for c in cluster_volume:
            # Coupe et NCut
            cut = inter_weights[c]
            vol = cluster_volume[c]
            ncut = cut / vol if vol > 0 else 0

            results.append({"cluster_id": c, "co_distance": intra_weights[c], "ncut": ncut, "nb_nodes": nb_nodes[c], "nb_edges": nb_edges[c]})

        df_final = pd.DataFrame(results)
        df_final.to_csv(output_csv + f"_{gamma}.csv", index=False)
        print(f"Terminé !")
        
     
# BLAST
traitement_graphe(BLAST_COMP_FILE, PRIN_COMP_CLUSTERS, OUTPUT_CSV)

# # BLAST & embeddings
# for com_file, prin_comp_file, output_csv in zip(com_files, prin_comp_files, output_csvs):
#     traitement_graphe(com_file, prin_comp_file, output_csv)