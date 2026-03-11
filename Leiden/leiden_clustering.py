# Script appliquant l'algorithme de clustering Leiden sur le fichier passé en paramètre
# -> Ici que sur la composante principale car les autres composantes sont déjà de bons clusters
BLAST_COMP_FILE = "/tempory/21234701/blast_components.tsv"                 # Fichier des composantes
OUTPUT_CLUSTERS = "/tempory/21234701/blast_princ_comp_leiden_clusters.tsv" # Clusters de sortie

import igraph as ig
import leidenalg as la
import numpy as np

# Récupère la composante principale (id 0)
print("Chargement de la composante")
edges_with_weights = []
with open(BLAST_COMP_FILE, "r") as f:
    next(f) # Sauter le header
    for line in f:
        parts = line.split()
        if len(parts) >= 4 and parts[3] == "0": # composante principale
            try:
                id1, id2 = parts[0], parts[1]
                score = float(parts[2])

                if score > 0:
                    edges_with_weights.append((id1, id2, score))
            except (ValueError, IndexError):
                continue

# Création du graphe igraph - indispensable pour Leiden
print(f"Création du graphe ({len(edges_with_weights)} arêtes)")
g = ig.Graph.TupleList(edges_with_weights, weights=True)

print(f"Exécution de l'algorithme de Leiden")
partition = la.find_partition(g, la.RBConfigurationVertexPartition, weights=g.es['weight'], resolution_parameter=1.0)

with open(OUTPUT_CLUSTERS, "w") as f_out:
    f_out.write("prot_id\tcluster_id\n")
    for cluster_id, nodes in enumerate(partition):
        for node_idx in nodes:
            prot_id = g.vs[node_idx]["name"]
            f_out.write(f"{prot_id}\t{cluster_id}\n")
                
print(f"Terminé ! {len(partition)} clusters créés")