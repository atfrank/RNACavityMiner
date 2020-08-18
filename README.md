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
git clone https://github.com/atfrank/RNACavityMiner.git
cd RNACavityMiner/
```
### Install Dependencies

#### Python Modules
```
conda create -n cavityminer --file environment.yml -c conda-forge -c bioconda -c schrodinger
conda activate cavityminer
pip install --user scikit-learn==0.22.2.post1
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
conda activate cavityminer
export CAVITYMINER="/path/to/RNACavityMiner/"
cd test/
./../src/miner.sh 1ANR_1.pdb
```
The predicted cavities and their corresponding scores are written in file `receptor/predicted_cavities.txt`.

#### Output
```
pdb,tag,cavityID,x,y,z,pred_MLP,pred_XGB,pred_RF,pred_LR,pred_Extra
1ANR_1,none,1,-6.161,7.090,2.111,0.111,0.204,0.548,0.600,0.490
1ANR_1,none,2,-5.423,-12.140,12.179,0.115,0.769,0.672,0.499,0.442
1ANR_1,none,3,-2.887,-5.340,7.370,0.294,0.665,0.800,0.915,0.673
1ANR_1,none,4,-11.601,2.035,13.648,0.000,0.003,0.161,0.016,0.265
1ANR_1,none,5,-20.377,-13.995,-0.740,0.000,0.010,0.098,0.005,0.151
1ANR_1,none,6,-14.886,2.647,7.492,0.000,0.001,0.137,0.002,0.249
1ANR_1,none,7,-7.620,7.864,17.366,0.000,0.001,0.085,0.004,0.107
1ANR_1,none,8,-21.181,0.314,7.673,0.000,0.002,0.083,0.027,0.145
1ANR_1,none,9,-13.190,-8.324,15.641,0.000,0.046,0.762,0.042,0.527
1ANR_1,none,10,5.000,3.750,12.750,0.000,0.002,0.128,0.019,0.243
1ANR_1,none,11,9.750,-3.000,5.750,0.000,0.001,0.131,0.001,0.215
```
#### Optional Arguments
```

## Publications

* (In preparation): Jingru Xie, and Aaron T. Frank. "Mining For Druggable Cavities in RNA"
