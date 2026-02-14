# Fichier : histogram_positions.py

import matplotlib.pyplot as plt
import numpy as np

align_file = "protein_top_tranches.tsv"  # ton fichier
output_image = "positions_histogram.png"

# Définir les bins pour l'histogramme (ajuster selon longueur max des protéines)
bins = np.arange(0, 2001, 10)
counts = np.zeros(len(bins) - 1, dtype=int)

# Lecture ligne par ligne et extraction de toutes les positions
with open(align_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # Supposons que les positions soient séparées par espace ou tab
        positions = line.replace(",", " ").split()
        for p in positions:
            try:
                pos = int(p)
                idx = np.searchsorted(bins, pos, side='right') - 1
                if 0 <= idx < len(counts):
                    counts[idx] += 1
            except ValueError:
                continue  # ignore les valeurs non numériques

# Tracer et sauvegarder
plt.figure(figsize=(10, 6))
plt.bar(bins[:-1], counts, width=10, edgecolor='black')
plt.xlabel("Position dans la protéine")
plt.ylabel("Nombre de fois où la position est importante")
plt.title("Distribution des positions importantes dans les alignements")
plt.tight_layout()
plt.savefig(output_image)
print(f"Histogramme enregistré dans : {output_image}")


"""# Fichier : positions_par_intervalle.py

align_file = "protein_top_tranches.tsv"  # ton fichier
output_file = "positions_par_intervalle.tsv"

interval_size = 1024  # taille de la fenêtre
counts = {}

with open(align_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # positions séparées par espace/tab/virgule
        positions = line.replace(",", " ").split()
        positions = [int(p) for p in positions if p.isdigit()]
        if not positions:
            continue
        start = min(positions)
        end = start + interval_size - 1
        interval = f"{start}-{end}"
        counts[interval] = counts.get(interval, 0) + 1

# Écrire le fichier
with open(output_file, "w") as f_out:
    for interval, n in sorted(counts.items(), key=lambda x: int(x[0].split("-")[0])):
        if n > 0:
            f_out.write(f"{interval}\t{n}\n")

print(f"Fichier des intervalles enregistré : {output_file}")
"""