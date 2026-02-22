# Fichier : positions_in_align.py
import gzip
from collections import defaultdict

ALIGNMENTS_FILE = "../../diamond_alignments.filter/diamond_alignments.filter.tsv"
OUTPUT_FILE = "protein_top_tranches.tsv"
TRANCHE_SIZE = 1024

# Si ton fichier est compressé : utiliser gzip.open
open_func = open
if ALIGNMENTS_FILE.endswith(".gz"):
    open_func = gzip.open

# Fonction pour calculer la tranche la plus couverte
def best_tranche(coverage_list, tranche_size=TRANCHE_SIZE):
    if not coverage_list:
        return None, 0
    min_pos = min(start for start, end in coverage_list)
    max_pos = max(end for start, end in coverage_list)
    length = max_pos + 1
    counts = [0] * length
    for start, end in coverage_list:
        for i in range(start-1, end):
            counts[i] += 1
    # Sliding window
    max_cov = sum(counts[:tranche_size])
    best_start = 1
    window_sum = max_cov
    for i in range(1, length - tranche_size):
        window_sum += counts[i + tranche_size - 1] - counts[i - 1]
        if window_sum > max_cov:
            max_cov = window_sum
            best_start = i + 1
    best_end = best_start + tranche_size - 1
    return best_start, best_end

# Stockage temporaire des positions pour une protéine à la fois
current_prot = None
current_cov = []
results = []

print("Lecture ligne par ligne et calcul des tranches...")

with open_func(ALIGNMENTS_FILE, "rt") as f:
    header = f.readline()  # ignorer l'en-tête
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        qseqid, qstart, qend = parts[0].strip(), int(parts[2]), int(parts[3])
        
        if current_prot is None:
            current_prot = qseqid

        if qseqid != current_prot:
            # calcul tranche précédente
            start, end = best_tranche(current_cov)
            results.append((current_prot, start, end))
            # reset
            current_prot = qseqid
            current_cov = []

        current_cov.append((qstart, qend))

# Dernière protéine
if current_prot is not None:
    start, end = best_tranche(current_cov)
    results.append((current_prot, start, end))

# Sauvegarde progressive
print(f"Sauvegarde dans {OUTPUT_FILE} ...")
with open(OUTPUT_FILE, "w") as fout:
    fout.write("protein\tstart\tend\n")
    for prot, start, end in results:
        fout.write(f"{prot}\t{start}\t{end}\n")

print("Terminé ! Résultats prêts pour les embeddings.")
