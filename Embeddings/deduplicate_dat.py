# Script utilisé pour supprimer les doublons dans le graphe final des embeddings
# (des doublons ont pu être créés s'il y a eu des erlances selon les checkpoints)

import os

# Paramètres
INPUT_DAT = "/tempory/21234701/output_proteinbert/embeddings_final_clean.dat"
OUTPUT_DAT = "/tempory/21234701/output_proteinbert/embeddings_final_deduplicated.dat"

ID_SIZE = 100
VEC_SIZE = 512
RECORD_SIZE = ID_SIZE + (VEC_SIZE * 4)

def deduplicate():
    if not os.path.exists(INPUT_DAT):
        print("Erreur : Fichier source introuvable.")
        return

    file_size = os.path.getsize(INPUT_DAT)
    num_records = file_size // RECORD_SIZE
    
    # Stockage {id: position_dans_le_fichier}
    # Si un ID apparaît deux fois, la nouvelle position écrasera l'ancienne
    print(f"Analyse des doublons sur {num_records} enregistrements")
    last_positions = {}
    
    with open(INPUT_DAT, "rb") as f:
        for i in range(num_records):
            f.seek(i * RECORD_SIZE)
            # Lecture de l'ID
            prot_id = f.read(ID_SIZE).decode('utf-8', errors='ignore').strip()
            last_positions[prot_id] = i
            
            if i % 500000 == 0:
                print(f"Lu : {i}/{num_records}...")

    unique_count = len(last_positions)
    
    # Trie des positions
    sorted_positions = sorted(last_positions.values())

    with open(INPUT_DAT, "rb") as f_in, open(OUTPUT_DAT, "wb") as f_out:
        for idx, pos in enumerate(sorted_positions):
            f_in.seek(pos * RECORD_SIZE)
            record = f_in.read(RECORD_SIZE)
            f_out.write(record)
            
            if idx % 500000 == 0:
                print(f"Copié : {idx}/{unique_count}")

    print("\nTerminé !")

if __name__ == "__main__":
    deduplicate()