import pandas as pd

diamond_file = "diamond_alignments.filter.tsv"
output_file = "graph_edges.tsv"

cols = [
    "qseqid", "qlen", "qstart", "qend",
    "sseqid", "slen", "sstart", "send",
    "length", "pident", "ppos",
    "score", "evalue", "bitscore"
]

chunksize = 200_000

with open(output_file, "w") as f_out:
    for chunk in pd.read_csv(diamond_file, sep="\t", names=cols, chunksize=chunksize):

        # filtres essentiels
        chunk = chunk[chunk["qseqid"] != chunk["sseqid"]]
        chunk = chunk[chunk["evalue"] < 1e-5]

        # couverture
        chunk["coverage"] = chunk["length"] / chunk[["qlen", "slen"]].min(axis=1)

        # poids final
        chunk["weight"] = (chunk["pident"] / 100) * chunk["coverage"]

        # filtrage optionnel pour alléger
        chunk = chunk[chunk["weight"] > 0.1]

        # écrire edge list
        for _, r in chunk.iterrows():
            f_out.write(f"{r['qseqid']}\t{r['sseqid']}\t{r['weight']}\n")

print("Terminé.")
