# Script utilisé pour générer le graphe avec la similarité cosinus des embeddings
# Pour chaque arête de BLAST, il calcule l'arête reliant les même sommets dans les embeddings

import numpy as np
import os

# Paramètres
DAT_FILE = "/tempory/21234701/output_proteinbert/embeddings_final_deduplicated.dat"
BLAST_FILE = "/tempory/21234701/diamond_alignments.filter.tsv"
OUTPUT_FILE = "/tempory/21234701/darkdino_cosine_graph.tsv"

ID_SIZE = 100           # Taille de l'ID (même ID que les .fasta)
VEC_SIZE = 512
RECORD_SIZE = ID_SIZE + (VEC_SIZE * 4)

def load_normalized_embeddings():
    """ Charge les vecteurs et les normalise pour transformer le produit scalaire en cosinus """
    print("Chargement et normalisation des vecteurs en RAM")
    emb_dict = {}
    file_size = os.path.getsize(DAT_FILE)
    num_records = file_size // RECORD_SIZE
    
    with open(DAT_FILE, "rb") as f:
        for i in range(num_records):
            chunk = f.read(RECORD_SIZE)
            if not chunk: break
            
            prot_id = chunk[:ID_SIZE].decode('utf-8', errors='ignore').strip()
            vec = np.frombuffer(chunk[ID_SIZE:], dtype=np.float32).copy()
            
            # Normalisation L2
            norm = np.linalg.norm(vec)
            if norm > 0:
                vec /= norm
            
            emb_dict[prot_id] = vec
            
            if i % 1000000 == 0:
                print(f"Chargé : {i}/{num_records} vecteurs")
    
    print("Chargement terminé")
    return emb_dict

def process_graph(emb_dict):
    """ Parcourt le fichier BLAST et calcule le cosinus pour chaque arête (-> ainsi on aura les mêmes arêtes dans BLAST et cosinus embeddings) """
    print("Traitement du graphe BLAST")
    count = 0
    found = 0
    
    with open(BLAST_FILE, "r") as f_in, open(OUTPUT_FILE, "w") as f_out:
        # Entête
        f_out.write("id1\tid2\tcosine_similarity\n")
        
        for line in f_in:
            if not line.strip(): continue
            parts = line.split()
            
            try:
                # Dans BLAST, ID1 est à l'index 0, ID2 est à l'index 4
                id1 = parts[0]
                id2 = parts[4]
                
                if id1 in emb_dict and id2 in emb_dict:
                    # Comme les vecteurs sont normalisés, Cosinus = Produit Scalaire
                    cos_sim = np.dot(emb_dict[id1], emb_dict[id2])
                    f_out.write(f"{id1}\t{id2}\t{cos_sim:.6f}\n")
                    found += 1
                
                count += 1
                if count % 1000000 == 0:
                    print(f"Lignes traitées : {count}... (Arêtes calculées : {found})", flush=True)
                    
            except IndexError:
                # Cas où la ligne ne respecte pas le format attendu
                continue

    print(f"\nTerminé !")
    print(f"Fichier généré : {OUTPUT_FILE}")
    print(f"Total arêtes traitées : {count}")
    print(f"Arêtes avec score cosinus : {found}")

if __name__ == "__main__":
    embeddings = load_normalized_embeddings()
    process_graph(embeddings)
