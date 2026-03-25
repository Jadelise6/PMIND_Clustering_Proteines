import networkx as nx

# Chemins
# BLAST
BLAST_FILE = "/tempory/21234701/graph_edges.tsv"                   # Graphe considéré : fichier d'entrée
BLAST_COMP_FILE = "/tempory/21234701/blast_components.tsv"         # Fichier prévu pour contenir les composantes de BLAST

# BLAST & Embeddings
ALPHA_LIST = [0.2, 0.4, 0.5, 0.6, 0.8]

INPUT_FILES = "/tempory/21234701/combined_graph_alpha_"              # Racine graphes considérés : fichiers d'entrée
input_files = [f"{INPUT_FILES}{a}.tsv" for a in ALPHA_LIST]

COMP_FILES = "/tempory/21234701/components_combined_graph_alpha_"    # Racine fichier prévu pour contenir les composantes des fichiers
output_files = [f"{COMP_FILES}{a}.tsv" for a in ALPHA_LIST]

def traitement_graphe(input_file, output_file):
    # Construction du graphe dans Networkx
    print("Construction du graphe pour identifier les composantes")
    G = nx.Graph()
    with open(input_file, "r") as f:
        for line in f:
            parts = line.split()
            
            # Ajout des arêtes
            if len(parts) >= 2:
                G.add_edge(parts[0], parts[1])

    # Composantes du graphe (le graphe n'étant pas connexe, il a plusieurs composantes...)
    print("Identification des composantes")
    components = list(nx.connected_components(G))

    # Trie des composantes par taille (la plus grande en ID 0)
    components.sort(key=len, reverse=True)

    # Dictionnaire associant à chaque nœud sa composante
    node_to_comp = {}
    for comp_id, nodes in enumerate(components):
        for node in nodes:
            node_to_comp[node] = comp_id

    # Écriture des composantes dans le fichier de sortie
    print("Écriture du fichier des composantes")
    # Relecture du fichier pour ajouter l'ID de la composante
    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        f_out.write("id1\tid2\tscore\tcomp_id\n")
        for line in f_in:
            parts = line.split()
            if len(parts) < 3:
                continue
            
            id1, id2, score = parts[0], parts[1], parts[2]
            comp_id = node_to_comp.get(id1) # id1 et id2 sont forcément dans la même comp
            f_out.write(f"{id1}\t{id2}\t{score}\t{comp_id}\n")

    print(f"Graphe complet: {G.number_of_nodes()} nœuds, {G.number_of_edges()} arêtes")
    print(f"Nombre de composantes : {nx.number_connected_components(G)}")
    print(f"Densité: {nx.density(G):.6f}")


# BLAST
# traitement_graphe(BLAST_FILE, BLAST_COMP_FILE)

# BLAST & embeddings
for input_file, output_file in zip(input_files, output_files):
    traitement_graphe(input_file, output_file)