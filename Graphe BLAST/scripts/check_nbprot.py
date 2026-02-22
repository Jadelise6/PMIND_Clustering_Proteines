# Fichier : check_nbprot.py

from pathlib import Path

fasta_dir = Path("../../PMIND2026_cluster/darkdino_fasta")  # dossier contenant tes .fasta
align_file = Path("../../diamond_alignments.filter/diamond_alignments.filter.tsv")

# --- Étape 1 : extraire tous les identifiants des FASTA ---
fasta_proteins = set()

for fasta_path in fasta_dir.glob("*.fasta"):
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                prot_id = line[1:].split()[0].replace(" ", "")  # retirer espace
                fasta_proteins.add(prot_id)

print(f"Nombre de protéines uniques dans FASTA : {len(fasta_proteins)}")

# --- Étape 2 : extraire identifiants de Diamond ligne par ligne ---
align_proteins = set()

with open(align_file) as f:
    header = next(f)  # sauter la ligne d'entête
    for line in f:
        cols = line.strip().split("\t")
        qseqid = cols[0].replace(" ", "")
        sseqid = cols[4].replace(" ", "")
        align_proteins.update([qseqid, sseqid])

print(f"Nombre de protéines uniques dans Diamond : {len(align_proteins)}")

# --- Étape 3 : comparer ---
missing_in_align = fasta_proteins - align_proteins
missing_in_fasta = align_proteins - fasta_proteins

print(f"Protéines FASTA absentes dans Diamond : {len(missing_in_align)}")
print(f"Protéines Diamond absentes dans FASTA : {len(missing_in_fasta)}")

# --- Étape 4 (optionnel) : montrer 10 exemples ---
print("Exemples FASTA absentes dans Diamond :", list(missing_in_align)[:10])
print("Exemples Diamond absentes dans FASTA :", list(missing_in_fasta)[:10])

# Stocker
with open("fasta_to_skip.txt", "w") as f:
    for p in missing_in_align:
        f.write(p + "\n")