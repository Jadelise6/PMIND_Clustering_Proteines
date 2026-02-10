# Fichier : run_mcl.py

import subprocess

# Chemin du graphe généré précédemment
graph_file = "graph_edges.tsv"

# Paramètres MCL
mcl_output = "graph_mcl.out"
inflation = 2.0  # paramètre classique à ajuster si besoin

# Commande MCL
# --abc : edge list au format a b w
# -I : inflation
# -o : output
cmd = [
    "mcl", graph_file,
    "--abc",
    "-I", str(inflation),
    "-o", mcl_output
]

print("Lancement de MCL...")
subprocess.run(cmd)
print(f"Clusters générés dans {mcl_output}")
