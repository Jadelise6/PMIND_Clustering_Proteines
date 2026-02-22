# Script utilisé pour générer les graphes résultants de la combinaison linéaire pondérée entre BLAST et la similarité cosinus
# des embeddings, selon plusieurs alpha

import os

# Paramètres
EMB_DAT = "/tempory/21234701/output_proteinbert/embeddings_final_deduplicated.dat"
BLAST_FILE = "/tempory/21234701/graph_edges.tsv"
ALPHA_LIST = [0.2, 0.4, 0.5, 0.6, 0.8]
COSINE_FILE = "/tempory/21234701/darkdino_cosine_graph.tsv"

ID_SIZE = 100
RECORD_SIZE = ID_SIZE + (512 * 4)

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

def combine_and_filter(valid_ids):
    print(f"Lancement pour alpha = {ALPHA_LIST}")
    
    # Création des fichiers de sortie
    outputs = {a: open(f"/tempory/21234701/combined_graph_alpha_{a}.tsv", "w") for a in ALPHA_LIST}
    
    with open(BLAST_FILE, "r") as f_blast, open(COSINE_FILE, "r") as f_cos:
        # Sauter les headers si nécessaire
        next(f_cos) 
        
        count = 0
        written = 0
        for line_b, line_c in zip(f_blast, f_cos):
            parts_b = line_b.split()
            parts_c = line_c.split()
            
            try:
                id1, id2 = parts_b[0], parts_b[1]
                
                if id1 in valid_ids and id2 in valid_ids:
                    s_blast = float(parts_b[2])
                    s_cosine = float(parts_c[2])
                    
                    for a in ALPHA_LIST:
                        score_final = (a * s_blast) + ((1 - a) * s_cosine)
                        outputs[a].write(f"{id1}\t{id2}\t{score_final:.6f}\n")
                    written += 1
                
                count += 1
                if count % 1000000 == 0:
                    print(f"Traitées : {count/1e6:.1f}M | Conservées : {written/1e6:.1f}M")
                    
            except (IndexError, ValueError):
                continue

    for f in outputs.values():
        f.close()
    print(f"Terminé ! Graphe final : {written} arêtes conservées.")

if __name__ == "__main__":
    valid_ids = get_valid_ids()
    combine_and_filter(valid_ids)