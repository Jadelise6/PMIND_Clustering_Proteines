# Ceci est le script qui a été lancé sur les ordis de la fac pour calculer les vecteurs embeddings


import torch, gc, os
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from transformers import AutoTokenizer, EsmModel

# Configuration
base_path = "/tempory/21234701" # stockage des résultats dans tempory/21234701
os.makedirs(base_path, exist_ok=True)

os.environ['HF_HOME'] = os.path.join(base_path, "hf_cache")

input_file = "all_darkdino.csv"
output_file = os.path.join(base_path, "embeddings.bin")
checkpoint_file = os.path.join(base_path, "progression.txt")

BATCH_SIZE = 256    # batch de protéines, on en traite plusieurs à la fois pour que ce soit plus rapide mais pas trop pour ne pas surcharger la RAM
MAX_LEN = 512       # longueur maximale prise en compte pour éviter la surcharge pour la RAM (sinon ça plante...)

# Chargement du modèle
# cuda = GPU NVIDIA, utilisée si présente car beaucoup plus rapide que CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Calcul sur : {device}")

# Modèle ESM_2
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = EsmModel.from_pretrained("facebook/esm2_t6_8M_UR50D").to(device)
model.eval()

# Chargement des données
print("Lecture du CSV")
df = pd.read_csv(input_file, header=None, usecols=[0])
seqs = df[0].astype(str).tolist()
total_seqs = len(seqs)
del df 
gc.collect()

print(f"Nombre de séquences chargées : {total_seqs}")

# Gestion de la reprise (si le script s'arrête il reprend où il s'était arrêté)
start_idx = 0
if os.path.exists(checkpoint_file):
    with open(checkpoint_file, 'r') as f:
        content = f.read().strip()
        if content:
            start_idx = int(content)
    print(f"Reprise du calcul à l'index : {start_idx}")

# Calcul
print(f"Traitement des protéines...")

with open(output_file, "ab") as f_out:
    with torch.no_grad():
        for i in tqdm(range(start_idx, total_seqs, BATCH_SIZE)):
            batch = [s[:MAX_LEN] for s in seqs[i : i + BATCH_SIZE]]
            
            # Tokenization nécessaire pour le modèle
            inputs = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=MAX_LEN).to(device)
            outputs = model(**inputs) # unpacking (déballage) du dictionnaire et utilisation du modèle
            
            emb = outputs.last_hidden_state.mean(dim=1).cpu().numpy().astype(np.float16)
            
            # Sauvegarde du vecteur
            f_out.write(emb.tobytes())
            
            # Sauvegarde la position (checkpoint) tous les 10 batches
            if (i // BATCH_SIZE) % 10 == 0:
                with open(checkpoint_file, 'w') as f_check:
                    f_check.write(str(i + BATCH_SIZE))
            
            del inputs, outputs
            
            # Gestion de la mémoire
            # Tous les 5 batches (donc ici tous les 1280 protéines), nettoie la RAM et libère la mémoire GPU
            if i % (BATCH_SIZE * 5) == 0:
                gc.collect()
                if device.type == 'cuda':
                    torch.cuda.empty_cache()

print(f"Terminé !")