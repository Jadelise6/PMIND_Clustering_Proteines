import glob
from collections import defaultdict
import os

# -----------------------
# CONFIG
# -----------------------

# Bases d'annotations fonctionnelles conservées
VALID_SOURCES = {
    "Pfam",
    "PANTHER",
    "SMART",
    "Gene3D",
    "SUPERFAMILY",
    "ProSiteProfiles",
    "ProSitePatterns",
    "PRINTS",
    "CDD",
    "NCBIfam",
    "PIRSF",
    "FunFam",
    "Hamap",
    "SFLD"
}

ANNOTATION_FILES = "data/annotations/darkdino_annotation/*.tsv"
OUTPUT_FILE = "outputs/annotations_filtered_full.tsv"
os.makedirs("outputs", exist_ok=True)

# -----------------------
# LECTURE ET FILTRAGE
# -----------------------

protein_annotations = defaultdict(list)  # liste pour garder plusieurs infos

files = glob.glob(ANNOTATION_FILES)

for file in files:
    print("Reading:", file)
    with open(file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue  # ignore les lignes mal formées

            protein = parts[0]
            source = parts[1]
            annot_id = parts[2]
            desc = parts[3] if len(parts) > 3 else ""

            if source in VALID_SOURCES:
                protein_annotations[protein].append({
                    "source": source,
                    "annot_id": annot_id,
                    "desc": desc
                })

# -----------------------
# ÉCRITURE DU FICHIER FILTRÉ
# -----------------------

with open(OUTPUT_FILE, "w") as out:
    out.write("prot_id\tsource\tannot_id\tdesc\n")  # header
    for protein, annots in protein_annotations.items():
        for annot in annots:
            out.write(f"{protein}\t{annot['source']}\t{annot['annot_id']}\t{annot['desc']}\n")

print("Filtering finished.")
print("Proteins annotated:", len(protein_annotations))