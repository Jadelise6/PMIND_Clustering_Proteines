# PMIND_Clustering_Proteines

## Structure du projet
```
├───Embeddings
│   │   create_all_darkdino.ipynb   # crée le fichier csv contenant les protéines à partir des génomes .fasta et affiche quelques statistiques
│   │   deduplicate_dat.py          # supprime les doublons dans le graphe final des embeddings
│   │   esm2_script.py              # crée les embeddings des protéines en utilisant le modèle ESM_2
│   │   generate_alpha_variants.py  # crée les graphes résultants de la combinaison linéaire pondérée entre BLAST et la similarité cosinus
│   │   generate_cos_graph.py       # crée le graphe avec la similarité cosinus des embeddings 
|   |   proteinBERT_script.py       # crée les embeddings des protéines en utilisant le modèle ProteinBERT
│   │
│   └───Data
│           darkdino_cleaned_for_bert.csv   # fichier csv contenant toutes les protéines sans X utilisé pour ProteinBERT
│           darkdino_cos_graph.tsv          # fichier tsv contenant le graph de la similarité cosinus des embeddings de ProteinBERT
│           embeddings_esm2.bin             # fichier binaire contenant les embeddings créé par ESM_2
│           embeddings_proteinbert.dat      # fichier dat contenant les embeddings créé par ProteinBERT
```
