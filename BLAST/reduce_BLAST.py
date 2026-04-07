# Calcule l'intersection entre BLAST et BLAST & Embeddings, afin de ne garder dans BLAST que les protéines de BLAST & Embeddings
import os

ID_SIZE = 100
RECORD_SIZE = ID_SIZE + (512 * 4)

# Chemins
# BLAST
BLAST_FILE = "/tempory/21234701/graph_edges.tsv"                   # Graphe considéré : fichier d'entrée
REDUCED_BLAST_FILE = "/tempory/21234701/reduced_graph_edges.tsv"   # Blast réduit

# Embeddings
EMB_DAT = "/tempory/21234701/output_proteinbert/embeddings_final_deduplicated.dat"

def get_valid_ids():
    print("Chargement de l'index des IDs des embeddings")
    valid_ids = set()
    with open(EMB_DAT, "rb") as f:
        file_size = os.path.getsize(EMB_DAT)
        num_records = file_size // RECORD_SIZE
        for i in range(num_records):
            chunk = f.read(RECORD_SIZE)
            prot_id = chunk[:ID_SIZE].decode('utf-8', errors='ignore').strip()
            valid_ids.add(prot_id)
            
    return valid_ids

def reduce_blast(blast_file, new_blast_file):
    valid_ids = get_valid_ids()
    
    # Création fichier de sortie
    f_blast_new = open(new_blast_file, "w")
    
    with open(blast_file, "r") as f_blast:
        for prot in f_blast:
            prot = prot.split()
            id1 = prot[0]
            id2 = prot[1]
            value = float(prot[2])
            
            if id1 in valid_ids and id2 in valid_ids:
                f_blast_new.write(f"{id1}\t{id2}\t{value:.6f}\n")

reduce_blast(BLAST_FILE, REDUCED_BLAST_FILE)