import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import train_test_split
from multiprocessing import Pool, cpu_count
import time

global_path = "/tempory/21234701"
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
knn_files = [f"{global_path}/kNN/kNN_combined_graph_alpha_{a}.tsv" for a in ALPHA_LIST]
output_files = [f"{global_path}/kNN/metrics_kNN_combined_graph_alpha_{a}.csv" for a in ALPHA_LIST]
FUNCTION_FILE = f"{global_path}/pfam_analysis/pfam_annotations.tsv"

GROUPED_DATA = {}
TRAIN_PROTS = set()
MAP_ANNOTATIONS = {}
K_VALUES = [3, 5, 7, 11]

def read_annotations_with_fallback(path):
    """ Résoud des problèmes d'encodage """
    try:
        return pd.read_csv(path, sep="\t", encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep="\t", encoding="latin-1")

def get_annotations(annot_file):
    print("Récupération des annotations")
    df = read_annotations_with_fallback(annot_file)
    return df[["prot_id", "annot_id"]].dropna(subset=["prot_id", "annot_id"]).groupby("prot_id")["annot_id"].apply(list).to_dict()

def process_protein_chunk(prot_ids):
    """ Calcule les prédictions pour un lot de protéines de test et pour tous les k d'un coup """
    local_results = {k: {"y_true": [], "y_pred": []} for k in K_VALUES}
    
    # Parcours des ids de protéines
    for prot_id in prot_ids:
        if prot_id not in GROUPED_DATA or prot_id not in MAP_ANNOTATIONS:
            continue
            
        true_annot = MAP_ANNOTATIONS[prot_id]
        neighbors_df = GROUPED_DATA[prot_id]
        
        # Filtrer pour ne garder que les voisins présents dans le train
        ids2 = neighbors_df["id2"].values
        scores = neighbors_df["score"].values
        annots_list = neighbors_df["annotations_list"].values
        
        valid_idx = [i for i, pid in enumerate(ids2) if pid in TRAIN_PROTS]
        if not valid_idx: continue
        
        # Trier par score décroissant
        valid_idx = sorted(valid_idx, key=lambda i: scores[i], reverse=True)
        max_k = max(K_VALUES)
        top_indices = valid_idx[:max_k]
        
        # Pour cette protéine, on calcule les votes cumulés au fur et à mesure pour chaque k
        current_votes = {}
        for idx_count, idx in enumerate(top_indices):
            score = float(scores[idx])
            for annot in annots_list[idx]:
                current_votes[annot] = current_votes.get(annot, 0) + score
                
            current_k = idx_count + 1
            if current_k in K_VALUES:
                if current_votes:
                    max_vote = max(current_votes.values())
                    pred_annot = [a for a, v in current_votes.items() if v >= 0.5 * max_vote]
                else:
                    pred_annot = []
                    
                local_results[current_k]["y_true"].append(true_annot)
                local_results[current_k]["y_pred"].append(pred_annot)
                
    return local_results

def run_knn_optimized(input_file, output_file):
    global GROUPED_DATA, TRAIN_PROTS, MAP_ANNOTATIONS
    print(f"Traitement de {input_file}")
    start_time = time.time()

    # Chargement et nettoyage
    df_knn = pd.read_csv(input_file, sep="\t", encoding="utf-8")
    df_knn["score"] = pd.to_numeric(df_knn["score"], errors="coerce")
    df_knn = df_knn.dropna(subset=["score"])
    df_knn["annotations_list"] = df_knn["annotations"].apply(lambda x: [a.strip() for a in str(x).split(";") if a.strip()])

    print("Création du dictionnaire")
    GROUPED_DATA = {prot_id: group for prot_id, group in df_knn.groupby("id1")}
    
    all_prots = list(GROUPED_DATA.keys())
    train_prots, test_prots = train_test_split(all_prots, test_size=0.2, random_state=42, shuffle=True)
    TRAIN_PROTS = set(train_prots)

    print(f"Protéines train: {len(TRAIN_PROTS)}, Protéines test: {len(test_prots)}")

    # Découper les protéines de test en lots pour les cœurs de calcul
    n_cores = cpu_count()
    print(f"Utilisation de {n_cores} cores pour paralléliser les protéines")
    chunks = np.array_split(test_prots, n_cores)

    with Pool(n_cores) as pool:
        chunk_results = pool.map(process_protein_chunk, chunks)

    # Fusionner les listes y_true et y_pred récoltées par les différents cœurs
    final_results = []
    for k in K_VALUES:
        combined_y_true = []
        combined_y_pred = []
        for res in chunk_results:
            combined_y_true.extend(res[k]["y_true"])
            combined_y_pred.extend(res[k]["y_pred"])
            
        if combined_y_true:
            mlb = MultiLabelBinarizer()
            mlb.fit(combined_y_true + combined_y_pred)
            y_true_bin = mlb.transform(combined_y_true)
            y_pred_bin = mlb.transform(combined_y_pred)

            metrics = {
                "k": k,
                "Recall": recall_score(y_true_bin, y_pred_bin, average="micro", zero_division=0),
                "Precision": precision_score(y_true_bin, y_pred_bin, average="micro", zero_division=0),
                "F1": f1_score(y_true_bin, y_pred_bin, average="micro", zero_division=0)
            }
        else:
            metrics = {"k": k, "Recall": 0, "Precision": 0, "F1": 0}
        final_results.append(metrics)

    pd.DataFrame(final_results).to_csv(output_file, index=False)
    elapsed = time.time() - start_time
    print(f"Temps de traitement : {elapsed:.1f}s ({elapsed/60:.1f}m)")

if __name__ == "__main__":
    MAP_ANNOTATIONS = get_annotations(FUNCTION_FILE)

    for input_file, output_file in zip(knn_files, output_files):
        run_knn_optimized(input_file, output_file)
        # Nettoyage de la mémoire
        GROUPED_DATA = {} 

    print("Terminé !")