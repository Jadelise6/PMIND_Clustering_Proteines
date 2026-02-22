# Script utilisé pour créer les Embeddings avec ProteinBERT

import os
import sys
import pandas as pd
import numpy as np
import pickle
import re
from sklearn.decomposition import IncrementalPCA

# Configuration du CPU
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

# Importation du modèle ProteinBERT
BASE_PATH = "/tempory/21234701/protein_bert-master"
if BASE_PATH not in sys.path:
    sys.path.insert(0, BASE_PATH)

from proteinbert.conv_and_global_attention_model import GlobalAttention
import proteinbert.conv_and_global_attention_model
proteinbert.conv_and_global_attention_model.ConvBertLayer = GlobalAttention
from proteinbert import load_pretrained_model

# Paramètres
INPUT_CSV = "/tempory/21234701/darkdino_cleaned_for_bert.csv" 
OUTPUT_DIR = "/tempory/21234701/output_proteinbert"
CHECKPOINT_FILE = os.path.join(OUTPUT_DIR, "checkpoint_pca_final.pkl")
FINAL_DAT = os.path.join(OUTPUT_DIR, "embeddings_final_clean.dat")

WINDOW_SIZE = 1024      # Taille max d'une protéine
ANNOTATION_DIM = 8943   # Dimension exigée par le modèle
REDUCED_DIM = 512       # Dimension réduite après le PCA
BATCH_SIZE_GPU = 32     # Taille des batchs envoyés à ProteinBERT (limité en mémoire)
ID_SIZE = 100           # Taille de l'ID (même ID que les .fasta)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Initialisation
print("Chargement de ProteinBERT", flush=True)
pretrained_model_generator, input_encoder = load_pretrained_model()
model = pretrained_model_generator.create_model(WINDOW_SIZE)

# Reprise des checkpoints et du PCA
if os.path.exists(CHECKPOINT_FILE):
    with open(CHECKPOINT_FILE, 'rb') as f:
        ckpt = pickle.load(f)
        pca = ckpt['pca']
        processed_count = ckpt['processed_count']
    print(f"Reprise au checkpoint : {processed_count}", flush=True)
    write_mode = "ab"
else:
    pca = IncrementalPCA(n_components=REDUCED_DIM)
    processed_count = 0
    write_mode = "wb"

def clean_seq(seq):
    """ Nettoie la séquence au préalable pour éviter les erreurs de ProteinBERT """
    s = re.sub(r'[^A-Z]', '', str(seq).upper())
    return s[:WINDOW_SIZE] if len(s) >= 5 else None

# Boucle principale
df_iter = pd.read_csv(INPUT_CSV, chunksize=1000, skiprows=range(1, processed_count + 1))

print("Démarrage du traitement", flush=True)

with open(FINAL_DAT, write_mode) as f_out:
    for chunk in df_iter:
        raw_ids = chunk['id'].tolist()
        raw_seqs = chunk['sequence'].tolist()
        
        chunk_all_embs = []
        chunk_valid_ids = []

        for i in range(0, len(raw_seqs), BATCH_SIZE_GPU):
            sub_ids = raw_ids[i : i + BATCH_SIZE_GPU]
            sub_seqs_raw = raw_seqs[i : i + BATCH_SIZE_GPU]
            
            batch_seqs = []
            batch_ids = []
            for s, tid in zip(sub_seqs_raw, sub_ids):
                cleaned = clean_seq(s)
                if cleaned:
                    batch_seqs.append(cleaned)
                    batch_ids.append(str(tid))

            if batch_seqs:
                try:
                    # Encodage brut
                    encoded_tuple = input_encoder.encode_X(batch_seqs, seq_len=WINDOW_SIZE)
                    
                    # Extraction et nettoyage manuel
                    indices_list = []
                    for item in encoded_tuple[0]:
                        # On force l'élément à être une liste plate d'entiers
                        if len(item) > 0 and isinstance(item[0], (list, np.ndarray)):
                             flat_item = list(item[0])
                        else:
                             flat_item = list(item)
                        
                        # On s'assure qu'elle fait exactement 1024
                        flat_item = flat_item[:WINDOW_SIZE]
                        if len(flat_item) < WINDOW_SIZE:
                            flat_item += [0] * (WINDOW_SIZE - len(flat_item))
                        
                        indices_list.append(flat_item)

                    X_seq = np.array(indices_list, dtype=np.int32)
                    
                    X_ann = np.zeros((len(batch_seqs), ANNOTATION_DIM), dtype=np.float32)
                    if len(encoded_tuple) > 1 and encoded_tuple[1] is not None:
                        try:
                            # On tente la conversion simple, sinon zéros
                            ann_data = np.array(encoded_tuple[1], dtype=np.float32)
                            if ann_data.shape[1] == ANNOTATION_DIM:
                                X_ann = ann_data
                        except:
                            pass
                    
                    # Prédiction
                    _, global_embs = model.predict([X_seq, X_ann], batch_size=len(batch_seqs), verbose=0)
                    
                    chunk_all_embs.extend(global_embs)
                    chunk_valid_ids.extend(batch_ids)
                    
                except Exception as e:
                    print(f"Erreur sous-batch index {processed_count + i}: {e}", flush=True)
                    
        if chunk_all_embs:
            if processed_count < 100000:
                pca.partial_fit(chunk_all_embs)
            
            reduced_embs = pca.transform(chunk_all_embs)
            
            for j in range(len(chunk_valid_ids)):
                header = chunk_valid_ids[j].ljust(ID_SIZE)[:ID_SIZE].encode('utf-8')
                f_out.write(header + reduced_embs[j].astype('float32').tobytes())

        processed_count += len(chunk)
        if processed_count % 5000 == 0:
            with open(CHECKPOINT_FILE, 'wb') as f_ck:
                pickle.dump({'pca': pca, 'processed_count': processed_count}, f_ck)
            print(f"Checkpoint : {processed_count} protéines traitées.", flush=True)