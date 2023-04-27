# RNACavityMiner
RNACavityMiner: Classifiers to Mine For Ligandable Binding Cavity in RNA

## Prerequisite
* [GCC 8.2.0]()

* [openbabel](http://openbabel.org/wiki/Category:Installation)

* [RxDock](https://www.rxdock.org/)

* [PyMOL](https://pymol.org/). For first time PyMOL users: you will need a [PyMOL license file](https://pymol.org/2/buy.html?q=buy), as PyMOL is a commercial software.


## Quick Start
```
git clone https://github.com/atfrank/RNACavityMiner.git
cd RNACavityMiner/
```
### Install Dependencies

#### Python Modules
```
conda create -n cavityminer python=3.7
conda activate cavityminer
conda install -c schrodinger pymol=2.4 -y
conda install -c conda-forge openbabel
conda install -c bioconda rxdock -y
conda install pandas numba tqdm -y
pip install --user --force scikit-learn==0.22.2.post1
pip install --user --force xgboost==0.90
```

#### RNAPosers
```
git clone https://github.com/atfrank/RNAPosers.git
cd RNAPosers/
make clean
make
cd ..
cp RNAPosers/bin/featurize bin/
rm -rf RNAPosers
```

### Using RNACavityMiner
Main script is `src/miner_grid.sh`. It requires a pdb file containing RNA 3D structure as input.
Example:
```
conda activate cavityminer
export CAVITYMINER="/path/to/RNACavityMiner/"
cd test/
./../src/miner_grid.sh 1ANR_1.pdb
```
The predicted cavities and their corresponding scores are written in file `receptor/predicted_cavities.txt`.

#### Output
```
pdb     tag   cavityID      x       y          z      pred_MLP  pred_XGB   pred_RF       pred_LR  pred_Extra                                                                                                      
1ANR_1  none         1   3.74   1.254   0.801000  4.322632e-02  0.106187  0.460681  3.785127e-01    0.555957
1ANR_1  none         2   3.74   1.254  12.801000  1.347721e-04  0.010728  0.250420  6.755266e-02    0.321917
1ANR_1  none         3   3.74   1.254  14.801000  1.475043e-04  0.022874  0.121178  1.792128e-02    0.225286
1ANR_1  none         4   3.74   1.254  16.800999  1.469739e-02  0.011660  0.088704  1.346384e-02    0.158415
1ANR_1  none         5   5.74   1.254  -9.199000  1.133171e-03  0.015931  0.071120  2.053839e-02    0.089329
...      ...       ...    ...     ...        ...           ...       ...       ...           ...         ...
1ANR_1  none      1766 -18.26  -2.746  -7.199000  2.709290e-03  0.003831  0.070224  1.418785e-02    0.088595
1ANR_1  none      1767 -18.26  -0.746  14.801000  7.025161e-09  0.003085  0.105943  2.861976e-09    0.171420
1ANR_1  none      1768 -22.26   1.254   6.801000  1.536864e-03  0.004555  0.076359  4.109725e-02    0.125883
1ANR_1  none      1769   3.74   1.254  -1.199000  1.534997e-01  0.118597  0.289568  3.087298e-01    0.392483
1ANR_1  none      1770 -12.26 -16.746  -5.199000  1.741724e-01  0.009482  0.088870  2.845870e-02    0.107240
```
## Publications
* Jingru Xie, and Aaron T. Frank. "Mining For Ligandable Cavities in RNA" (https://doi.org/10.1021/acsmedchemlett.1c00068)

## Data and ML-Code
* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4049068.svg)](https://doi.org/10.5281/zenodo.4049068)



