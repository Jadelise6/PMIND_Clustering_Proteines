# Réexécute l'algorithme de MCL, cette fois-ci pour toutes les composantes, et calcule le temps passé
# pour afficher une courbe
# Reprend donc le code de mcl_clustering.py
import multiprocessing
import subprocess
import os
import time
import csv
import pandas as pd

# Chemin vers le mcl installé avec conda
mcl_path = "/tempory/21234701/pmind_conda/bin/mcl"

# Nombre de cœurs de la machine
threads = multiprocessing.cpu_count()

# BLAST
BLAST_COMP_FILE = "/tempory/21234701/reduced_blast_components.tsv"       # Fichier des composantes
output_file = "/tempory/21234701/time_mcl/time_temporary_mcl"            # Fichier temporaire
time_file = "./time_mcl.csv"                                             # Fichier CSV des temps
top_components_file = "./top_components_blast.csv"                       # Fichier temporaire des 20 plus grandes composantes
inflation = 2.0

# Listes des temps d'exécution et tailles
writing_input_times = []
mcl_times = []
writing_output_times = []

df_components = pd.read_csv(top_components_file) # comp_id, nb_nodes

def traitement_graphe(input_file, output_file, comp_id):
    # On crée un fichier temporaire sans la 4ème colonne
    writing_input_start_time = time.time()
    
    print(f"Réécriture du fichier de la composante {comp_id}")
    clean_abc = f"{input_file}_comp{comp_id}.clean.abc"
    
    print(f"Filtrage de la composante {comp_id} vers {clean_abc}")
    comp_size = 0
    
    with open(input_file, 'r') as f_in, open(clean_abc, 'w') as f_out:
        next(f_in) # Skip header
        for line in f_in:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[3] == str(comp_id):
                id1, id2, score = parts[0], parts[1], parts[2]
                # On écrit l'arête dans les deux sens pour MCL
                f_out.write(f"{id1}\t{id2}\t{score}\n")
                f_out.write(f"{id2}\t{id1}\t{score}\n")
                comp_size += 1
                
    writing_input_times.append(time.time() - writing_input_start_time)
               
    final_tsv = f"{output_file}_comp{comp_id}.tsv"
    temp_mcl_out = f"{output_file}_comp{comp_id}_raw.tmp"
    
    # Application de MCL avec inflation=2.0
    mcl_start_time = time.time()
    print(f"Lancement MCL - Inflation {inflation}")
    cmd = [
        mcl_path, clean_abc,
        "--abc",
        "-I", str(inflation),
        "-te", str(threads),
        "-abc-tf", "mul(10)",
        "-o", temp_mcl_out
    ]
    subprocess.run(cmd, check=True)
    mcl_times.append(time.time() - mcl_start_time)

    writing_output_start_time = time.time()
    # Conversion une ligne par cluster -> prot_id \t cluster_id
    if os.path.exists(temp_mcl_out):
        with open(temp_mcl_out, 'r') as f_in, open(final_tsv, 'w') as f_out:
            f_out.write("prot_id\tcluster_id\n")
            for cluster_id, line in enumerate(f_in):
                proteins = line.strip().split('\t')
                for prot in proteins:
                    f_out.write(f"{prot}\t{cluster_id}\n")
        
        os.remove(temp_mcl_out) # Nettoyage
    writing_output_times.append(time.time() - writing_output_start_time)
                
    print(f"Terminé pour composante {comp_id} !")

for comp_id in df_components["comp_id"]: # 20 premières composantes présentes dans le fichier
    traitement_graphe(BLAST_COMP_FILE, output_file, comp_id)

# Sauvegarder tous les temps dans un fichier CSV
with open(time_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['nb_nodes_comp', 'writing_input_time', 'mcl_time', 'writing_output_time'])
    for i in range(len(mcl_times)):
        writer.writerow([
            df_components["nb_nodes"].iloc[i],
            writing_input_times[i],
            mcl_times[i],
            writing_output_times[i]
        ])

print(f"Données sauvegardées dans {time_file}")