# RNACavityMiner
RNACavityMiner: Classifiers to Mine For Druggable Binding Cavity in RNA

## Prerequisite
* [GCC 8.2.0]()

* [openbabel](http://openbabel.org/wiki/Category:Installation)

* [RxDock](https://www.rxdock.org/)

* [PyMOL](https://pymol.org/). For first time PyMOL users: you will need a [PyMOL license file](https://pymol.org/2/buy.html?q=buy), as PyMOL is a commercial software.

* see `environment.yml`

## Quick Start
```
git clone git@github.com:karoka/RNACavityMiner.git
cd RNACavityMiner/
```
### Install Dependencies

#### Python Modules
```
conda create -n cavityminer --file environment.yml
conda activate cavityminer
```

#### RNAPosers
```
git clone --depth=1 git@github.com:atfrank/RNAPosers.git
cd RNAPosers/
make clean
make
cd ..
cp RNAPosers/bin/featurize bin/
rm -rf RNAPosers
```

### Using RNACavityMiner
Main script is `src/miner.sh`. It requires a pdb file containing RNA 3D structure as input.
Example:
```
cd test/
./../src/miner.sh receptor.pdb
```
The predicted cavities and their corresponding scores are written in file `receptor/predicted_cavities.txt`.

#### Output
| pdb          | tag   | cavityID | x     | y     | z      | pred_MLP | pred_XGB | pred_RF | pred_LR | pred_Extra|
|--------------|-------|----------|-------|-------|--------|----------|----------|---------|----------|----------|
| receptor.pdb | decoy | 1        | 2.835 | 0.475 | -0.482 | 1.000    | 0.886    | 0.705   | 1.000    | 0.886    |

#### Optional Arguments
```

## Publications

* (In preparation): Jingru Xie, and Aaron T. Frank. "Mining For Druggable Cavities in RNA"
