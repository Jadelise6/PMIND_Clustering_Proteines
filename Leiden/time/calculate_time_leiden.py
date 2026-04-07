# Réexécute l'algorithme de Leiden, cette fois-ci pour toutes les composantes, et calcule le temps passé
# pour afficher une courbe
# Reprend donc le code de leiden_clustering.py

import igraph as ig
import leidenalg as la
import time
import csv

# BLAST
BLAST_COMP_FILE = "/tempory/21234701/reduced_blast_components.tsv"       # Fichier des composantes
output_file = "/tempory/21234701/time_temporary_leiden"                  # Fichier temporaire
time_file = "./time_leiden.csv"                                          # Fichier CSV des temps
gamma = 1.0

# Listes des temps d'exécution et tailles
component_sizes = []
loading_times = []
graph_times = []
leiden_times = []
writing_times = []

def traitement_graphe(input_file, output_file, comp_id):
    # Récupère la composante
    loading_start_time = time.time()
    print(f"Chargement de la composante {comp_id}")
    edges_with_weights = []
    with open(input_file, "r") as f:
        next(f) # Sauter le header
        for line in f:
            parts = line.split()
            if len(parts) >= 4 and parts[3] == str(comp_id):
                try:
                    id1, id2 = parts[0], parts[1]
                    score = float(parts[2])

                    if score > 0:
                        edges_with_weights.append((id1, id2, score))
                except (ValueError, IndexError):
                    continue
    loading_times.append(time.time() - loading_start_time)
               
    # Création du graphe igraph - indispensable pour Leiden
    graph_start_time = time.time()
    g = ig.Graph.TupleList(edges_with_weights, weights=True)
    graph_times.append(time.time() - graph_start_time)
    
    # Taille de la composante (nombre de sommets)
    comp_size = g.vcount()
    component_sizes.append(comp_size)

    # Application de Leiden avec gamma=1.0
    leiden_start_time = time.time()
    partition = la.find_partition(g, la.RBConfigurationVertexPartition, weights=g.es['weight'], resolution_parameter=gamma, seed=42)
    leiden_times.append(time.time() - leiden_start_time)

    writing_start_time = time.time()
    with open(output_file + f"_{comp_id}.tsv", "w") as f_out:
        f_out.write("prot_id\tcluster_id\n")
        for cluster_id, nodes in enumerate(partition):
            for node_idx in nodes:
                prot_id = g.vs[node_idx]["name"]
                f_out.write(f"{prot_id}\t{cluster_id}\n")
    writing_times.append(time.time() - writing_start_time)
                
    print(f"Terminé pour composante {comp_id} ! Taille: {comp_size} sommets, {len(partition)} clusters créés")
        

for comp_id in range(20): # 20 première composantes pour avoir une bonne estimation
    traitement_graphe(BLAST_COMP_FILE, output_file, comp_id)

# Sauvegarder tous les temps dans un fichier CSV
with open(time_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['component_size', 'loading_time', 'graph_time', 'leiden_time', 'writing_time'])
    for i in range(len(component_sizes)):
        writer.writerow([
            component_sizes[i],
            loading_times[i],
            graph_times[i],
            leiden_times[i],
            writing_times[i]
        ])

print(f"Données sauvegardées dans {time_file}")