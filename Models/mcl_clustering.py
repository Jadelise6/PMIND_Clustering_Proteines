# Script appliquant l'algorithme de clustering MCL sur le fichier passé en paramètre
# -> Ici que sur la composante principale car les autres composantes sont déjà de bons clusters
import multiprocessing
import subprocess
import os

# Chemins
global_path = "/tempory/21234701"

# Chemin vers le mcl installé avec conda
mcl_path = f"{global_path}/pmind_conda/bin/mcl"

# Nombre de cœurs de la machine
threads = multiprocessing.cpu_count()

# BLAST
BLAST_COMP_FILE = f"{global_path}/reduced_blast_components.tsv"                      # Fichier des composantes
OUTPUT_CLUSTERS_BLAST = f"{global_path}/mcl/reduced_blast_princ_comp_mcl_clusters"   # Clusters de sortie

# BLAST & Embeddings
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 0.93, 0.96, 0.99, 1.0]
COMP_FILES = f"{global_path}/components_combined_graph_alpha_"    # Racine fichier prévu pour contenir les composantes des fichiers
input_files = [f"{COMP_FILES}{a}.tsv" for a in ALPHA_LIST]

OUTPUT_CLUSTERS = f"{global_path}/mcl/combined_graph_alpha_princ_comp_mcl"
output_files = [f"{OUTPUT_CLUSTERS}{a}" for a in ALPHA_LIST]

def traitement_graphe(input_file, output_base):
    # On crée un fichier temporaire sans la 4ème colonne
    clean_abc = input_file + ".clean.abc"
    print(f"Filtrage de la composante 0 vers {clean_abc}")
    
    with open(input_file, 'r') as f_in, open(clean_abc, 'w') as f_out:
        next(f_in) # Skip header
        for line in f_in:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[3] == "0":
                id1, id2, score = parts[0], parts[1], parts[2]
                # On écrit l'arête dans les deux sens pour MCL
                f_out.write(f"{id1}\t{id2}\t{score}\n")
                f_out.write(f"{id2}\t{id1}\t{score}\n")

    # Appel de MCL
    inflation_vals = [1.3, 1.5, 2.0, 2.5, 3.0]
    for inflation in inflation_vals:
        final_tsv = f"{output_base}_{inflation}.tsv"
        temp_mcl_out = f"{output_base}_{inflation}_raw.tmp"
        
        print(f"Lancement MCL - Inflation {inflation}")
        cmd = [
            mcl_path, clean_abc,
            "--abc",
            "-I", str(inflation),
            "-te", str(threads),
            "-abc-tf", "mul(10)",
            "-o", temp_mcl_out
        ]
        
        subprocess.run(cmd)

        # Conversion une ligne par cluster -> prot_id \t cluster_id
        if os.path.exists(temp_mcl_out):
            with open(temp_mcl_out, 'r') as f_in, open(final_tsv, 'w') as f_out:
                f_out.write("prot_id\tcluster_id\n")
                for cluster_id, line in enumerate(f_in):
                    proteins = line.strip().split('\t')
                    for prot in proteins:
                        f_out.write(f"{prot}\t{cluster_id}\n")
            
            os.remove(temp_mcl_out) # Nettoyage
            print(f"Terminé pour inflation={inflation}")
        
# BLAST
#traitement_graphe(BLAST_COMP_FILE, OUTPUT_CLUSTERS_BLAST)

# BLAST & embeddings
for input_file, output_file in zip(input_files, output_files):
    traitement_graphe(input_file, output_file)