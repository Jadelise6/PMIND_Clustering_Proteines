# Ceci est le script qui a été lancé sur les ordis de la fac pour calculer les vecteurs embeddings avec ESM-2

import torch, os, gc
import numpy as np
import pandas as pd
from transformers import AutoTokenizer, EsmModel
from tqdm import tqdm

# Configuration
bin_file = "/tempory/21234701/embeddings_7M.bin"
csv_file = "all_darkdino.csv"
CHUNK_SIZE = 1024               # longueur d'une séquence
EMB_DIM = 320                   # dimension des embeddings
BYTES_PER_EMB = EMB_DIM * 2     # bytes par embedding

# Utiliser cuda est préférable s'il est présent -> plus rapide
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Chargement du modèle
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = EsmModel.from_pretrained("facebook/esm2_t6_8M_UR50D").to(device)
model.eval()

# Chargement des séquences
print("Chargement des séquences")
df = pd.read_csv(csv_file, header=None, usecols=[0])
seqs = df[0].astype(str).tolist()
del df
indices_a_patcher = [i for i, s in enumerate(seqs) if len(s) > 512]

# Déterminer où reprendre la lecture
# On regarde combien de protéines sont déjà dans le fichier .bin
if os.path.exists(bin_file):
    file_size = os.path.getsize(bin_file)
    checkpoint_path = "progression_patch.txt"
    
    if os.path.exists(checkpoint_path):
        with open(checkpoint_path, 'r') as f_prog:
            content = f_prog.read().strip()
            if content:
                start_list_idx = int(content)
                print(f"Reprise à l'index : {start_list_idx}")
            else:
                start_list_idx = 0
    else:
        start_list_idx = 0
        
print(f"Reprise à l'index de liste : {start_list_idx} / {len(indices_a_patcher)}")

def get_safe_embedding(sequence, model, tokenizer):
    """ Calcul des embeddings """
    with torch.no_grad():
        # Si la séquence est courte, calcul direct
        if len(sequence) <= CHUNK_SIZE:
            inputs = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=CHUNK_SIZE).to(device)
            outputs = model(**inputs)
            return outputs.last_hidden_state.mean(dim=1).cpu().numpy().astype(np.float16)
        else:
            # Si la séquence est longue (>1024), on utilise la fenêtre glissante
            chunks = [sequence[i : i + CHUNK_SIZE] for i in range(0, len(sequence), CHUNK_SIZE)]
            chunk_embeddings = []
            for chunk in chunks:
                inputs = tokenizer(chunk, return_tensors="pt", truncation=True, max_length=CHUNK_SIZE).to(device)
                outputs = model(**inputs)
                chunk_embeddings.append(outputs.last_hidden_state.mean(dim=1).cpu().numpy())
                
            # Moyenne des morceaux pour obtenir l'embedding global
            return np.mean(chunk_embeddings, axis=0).astype(np.float16)

with open(bin_file, "r+b") as f:
    for i in tqdm(range(start_list_idx, len(indices_a_patcher))):
        real_idx = indices_a_patcher[i]
        sequence = seqs[real_idx]
        
        try:
            # Calcul des embeddings
            new_emb = get_safe_embedding(sequence, model, tokenizer)
            
            # Injection dans le fichier
            position = real_idx * 320 * 2
            f.seek(position)
            f.write(new_emb.tobytes())
            
            # Sauvegarde (incrémentation de l'index pour la reprise...)
            with open(checkpoint_path, 'w') as f_prog:
                f_prog.write(str(i + 1))
            
            # Forçage de l'écriture disque
            f.flush()
            os.fsync(f.fileno())

        except Exception as e:
            # Si une protéine fait planter le modèle, on log l'erreur et on avance
            print(f"\nErreur sur la protéine à la liste-index {i} (index CSV {real_idx}): {e}")
            with open(checkpoint_path, 'w') as f_prog:
                f_prog.write(str(i + 1))
            continue

        if i % 100 == 0:
            gc.collect()
            torch.cuda.empty_cache()
