# Ce script va créer un fichier csv utilisable pour le kNN
import pandas as pd

# Chemins
global_path = "/tempory/21234701"

# BLAST & Embeddings
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]

INPUT_FILES = f"{global_path}/combined_graph_alpha_"                    # Racine graphes considérés : fichiers d'entrée
input_files = [f"{INPUT_FILES}{a}.tsv" for a in ALPHA_LIST]

KNN_FILES = f"{global_path}/kNN/kNN_combined_graph_alpha_"              # Racine fichier prévu pour contenir les fichiers utilisés par le KNN
output_files = [f"{KNN_FILES}{a}.tsv" for a in ALPHA_LIST]

FUNCTION_FILE = f"{global_path}/pfam_analysis/pfam_annotations.tsv"     # Fichier des annotations, format : prot_id source  annot_id  desc


def read_annotations_with_fallback(path):
    """ Résolution des problèmes d'encodage """
    try:
        return pd.read_csv(path, sep="\t", encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep="\t", encoding="latin-1")

def get_annotations(annot_file):
    print("Récupération des annotations")
    # Lecture du fichier
    df_annotations = read_annotations_with_fallback(annot_file)
    
    # Transformation en dictionnaire {idProt: [annotation1, annotation2, ...]}
    map_annotations = (
        df_annotations[["prot_id", "annot_id"]]
        .dropna(subset=["prot_id", "annot_id"])
        .groupby("prot_id")["annot_id"]
        .apply(lambda series: list(dict.fromkeys(series.astype(str))))
        .to_dict()
    )
    
    return map_annotations
    
def traitement_graphe(input_file, output_file, map_annotations):
    print("Lecture et écriture des données")
    with open(input_file, "r", encoding="utf-8", errors="replace") as f_in, open(output_file, "w", encoding="utf-8") as f_out:
        f_out.write("id1\tid2\tscore\tannotations\n")
        for line in f_in:
            parts = line.split()
            
            if len(parts) >= 2:
                prot_1 = parts[0]
                prot_2 = parts[1]
                score = parts[2] if len(parts) >= 3 else ""
                
                # Récupérer les annotations des deux protéines
                annotation_1 = map_annotations.get(prot_1)
                annotation_2 = map_annotations.get(prot_2)
                
                # Écriture dans le nouveau fichier
                if annotation_2:
                    f_out.write(f"{prot_1}\t{prot_2}\t{score}\t{';'.join(annotation_2)}\n")
                if annotation_1:
                    f_out.write(f"{prot_2}\t{prot_1}\t{score}\t{';'.join(annotation_1)}\n")
          
# Récupère les annotations et les stocke dans un dictionnaire en RAM      
map_annotations = get_annotations(FUNCTION_FILE)

# BLAST & embeddings
for input_file, output_file in zip(input_files, output_files):
    traitement_graphe(input_file, output_file, map_annotations)
    
print("Terminé !")