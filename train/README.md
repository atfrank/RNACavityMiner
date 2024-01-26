# RNACavityMiner: Mining For Druggable Cavities in RNA

## Install Dependencies

### Python Modules
```
conda create -n cavityminer_train --file environment.yml -c conda-forge -c bioconda
conda activate cavityminer_train
pip install --user scikit-learn==0.22.2.post1
```

## (A) t-SNE Analysis of the Feature Space
```
# E.g., using features generated with cutoff=20.0 and eta=4
conda activate cavityminer
mkdir tsne/
python tsne.py --infile features/rna-ligand-pocket-20-4.csv --outfile tsne/info-20-4.csv
```

## (B) Train and Test Cavity Classifiers
```
# E.g., using features generated with cutoff=20.0 and eta=4
conda activate cavityminer_train
mkdir models/
python train.py --infile features/rna-ligand-pocket-20-4.csv --tag final_rc_20-eta-4
```
