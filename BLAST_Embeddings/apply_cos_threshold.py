# Script utilisé pour appliquer un seuil au graphe avec la similarité cosinus des embeddings

# Chemins
global_path = "/tempory/21234701"

COS_FILE = f"{global_path}/darkdino_cosine_graph.tsv"
OUTPUT_FILE = f"{global_path}/darkdino_cos_graph_threshold.tsv"

threshold = 0.1         # Seuil identique à celui de BLAST

with open(COS_FILE, "r") as f_in, open(OUTPUT_FILE, "w") as f_out:
    # Entête
    f_out.write("id1\tid2\tcos\n")
    
    lignes = f_in.readlines()
    for line in lignes[1:]:
        if not line.strip():
            continue
        parts = line.split()
        
        try:
            id1 = parts[0]
            id2 = parts[1]
            cos = float(parts[2])
            
            # Application du seuil
            if cos > 0.1:
                f_out.write(f"{id1}\t{id2}\t{cos:.6f}\n")
         
        except IndexError:
            # Cas où la ligne ne respecte pas le format attendu
            continue

print(f"\nTerminé !")