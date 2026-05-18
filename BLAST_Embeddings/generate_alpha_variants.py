# Script utilisé pour générer les graphes résultants de la combinaison linéaire pondérée entre BLAST et la similarité cosinus
# des embeddings, selon plusieurs alpha

import os

# Chemins
global_path = "/tempory/21234701"

EMB_DAT = f"{global_path}/output_proteinbert/embeddings_final_deduplicated.dat"
BLAST_FILE = f"{global_path}/graph_edges.tsv"
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 0.93, 0.96, 0.99, 1.0]
COSINE_FILE = f"{global_path}/darkdino_cos_graph_threshold.tsv"

ID_SIZE = 100
RECORD_SIZE = ID_SIZE + (512 * 4)

def get_valid_ids():
    # Crée un ensemble des IDs des embeddings, qui sont donc les Ids valides pour BLAST
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
    
    # Chargement du fichier de la similaité cosinus des embeddings
    print("Chargement du fichier cosinus")
    cosine_scores = {}
    
    with open(COSINE_FILE, "r") as f_cos:
        next(f_cos)  # Skip header
        for line in f_cos:
            parts = line.split()
            if len(parts) >= 3:
                
                try:
                    id1, id2 = parts[0], parts[1]
                    score = float(parts[2])
                    cosine_scores[(id1, id2)] = score
                    cosine_scores[(id2, id1)] = score
                    
                except (ValueError, IndexError):
                    continue
        
    # Création des fichiers de sortie
    outputs = {a: open(f"/tempory/21234701/combined_graph_alpha_{a}.tsv", "w") for a in ALPHA_LIST}
    
    # Parcourir BLAST et chercher les scores cosinus correspondants
    print("Chargement du fichier BLAST et combinaison linéaire pondérée")
    count = 0
    written = 0
    
    with open(BLAST_FILE, "r") as f_blast:
        for line_b in f_blast:
            parts_b = line_b.split()
            
            try:
                id1, id2 = parts_b[0], parts_b[1]
                s_blast = float(parts_b[2])
                
                if id1 in valid_ids and id2 in valid_ids:
                    # Chercher le score cosinus - 0 si absent
                    s_cosine = cosine_scores.get((id1, id2), 0.0)
                    
                    # Combinaison linéaire pondérée avec alpha
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
    print(f"Terminé !")

valid_ids = get_valid_ids()
combine_and_filter(valid_ids)