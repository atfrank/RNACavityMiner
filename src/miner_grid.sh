#!/bin/bash

if [[ $# -ne 1 ]]
then
    echo "Path to pdb file required."
    echo "Usage: $0 <path-to-pdb-file>"
    exit 1
fi

# initialize variable
cutoff=20.0
selatms=":ADE.C1' :ADE.C2 :ADE.C2' :ADE.C3' :ADE.C4 :ADE.C4' :ADE.C5 :ADE.C5' :ADE.C6 :ADE.C8 :ADE.N1 :ADE.N3 :ADE.N6 :ADE.N7 :ADE.N9 :ADE.O2' :ADE.O3' :ADE.O4' :ADE.O5' :ADE.OP1 :ADE.OP2 :ADE.P :CYT.C1' :CYT.C2 :CYT.C2' :CYT.C3' :CYT.C4 :CYT.C4' :CYT.C5 :CYT.C5' :CYT.C6 :CYT.N1 :CYT.N3 :CYT.N4 :CYT.O2 :CYT.O2' :CYT.O3' :CYT.O4' :CYT.O5' :CYT.OP1 :CYT.OP2 :CYT.P :GUA.C1' :GUA.C2 :GUA.C2' :GUA.C3' :GUA.C4 :GUA.C4' :GUA.C5 :GUA.C5' :GUA.C6 :GUA.C8 :GUA.N1 :GUA.N2 :GUA.N3 :GUA.N7 :GUA.N9 :GUA.O2' :GUA.O3' :GUA.O4' :GUA.O5' :GUA.O6 :GUA.OP1 :GUA.OP2 :GUA.P :URA.C1' :URA.C2 :URA.C2' :URA.C3' :URA.C4 :URA.C4' :URA.C5 :URA.C5' :URA.C6 :URA.N1 :URA.N3 :URA.O2 :URA.O2' :URA.O3' :URA.O4 :URA.O4' :URA.O5' :URA.OP1 :URA.OP2 :URA.P"
rowatms=":UNK."

# get source directory
DIRMODEL=${CAVITYMINER}/models
DIRPY=${CAVITYMINER}/py
DIRSCR=${CAVITYMINER}/src
DIRBIN=${CAVITYMINER}/bin
rdockParam=${CAVITYMINER}/params/blind_docking_de_novo.prm
pml=${CAVITYMINER}/py/fix_receptor.pml
predictPy=${CAVITYMINER}/py/predict.py
model=${DIRMODEL}/models.pkl
cavityID="none"

# get command line arguments
initialPDB=$1
rna=`basename $initialPDB | sed 's/.pdb//g'`

# setup working directory
workingDIR="`pwd`/${rna}"
rm -rfv ${workingDIR}
mkdir -p ${workingDIR}

# goto working directory
cp -v ${initialPDB} ${workingDIR}/.
cd ${workingDIR}
outFile=${workingDIR}/predicted_cavities.txt
featureFile=${workingDIR}/features.csv
cavityFile=${workingDIR}/receptor_grid.xyz

# generate grid file
echo "started making grids.."
python ${CAVITYMINER}/py/grid_maker.py -c ${initialPDB} -o receptor
echo "finished making grids.."

# generate MOL2
echo "started converting coordinate file to MOL2.."
obabel -ipdb ${initialPDB} -omol2 -O receptor.mol2 #&> /dev/null
echo "finished converting coordinate file to MOL2.."

# get center
echo "started determining center.."
pymol -cqr $pml &> receptor.log
coorx=`grep "CENTER X" receptor.log | grep -v PyMOL | awk '{print $4}'`
coory=`grep "CENTER Y" receptor.log | grep -v PyMOL | awk '{print $4}'`
coorz=`grep "CENTER Z" receptor.log | grep -v PyMOL | awk '{print $4}'`
# store cavity center
cavx=${coorx}
cavy=${coory}
cavz=${coorz}
echo "finished determining center.."

if [[ -s ${cavityFile} ]]
then
    echo "preparting to featurize.."
    ncav=`wc -l ${cavityFile} | awk '{print $1}'`
    echo "${ncav}         " > cavity_decoy.xyz
    echo "decoys" >> cavity_decoy.xyz
    cat ${cavityFile} >> cavity_decoy.xyz
    obabel -ixyz cavity_decoy.xyz -opdb -O cavity_decoy.pdb #&> /dev/null
    sed 's/UNL/UNK/g' cavity_decoy.pdb | sed 's/LIG/UNK/g' > cavity_decoy_tmp.pdb
    mv cavity_decoy_tmp.pdb cavity_decoy.pdb
    
    obabel -ipdb cavity_decoy.pdb -omol2 -O cavity_decoy.mol2 #&> /dev/null
    
    # create complex native
    grep -v 'TER' receptor.pdb | grep -v 'END' > complex_cavity_decoy.pdb
    echo "TER" >> complex_cavity_decoy.pdb
    grep 'HETATM' cavity_decoy.pdb >> complex_cavity_decoy.pdb
    echo "END" >> complex_cavity_decoy.pdb
    
    # featurize decoy cavity
    echo "started featurizing.."
    ${CAVITYMINER}/bin/featurize \
        -etaStartPow 1 \
        -numEta 4 \
        -cutoff ${cutoff} \
        -scalar 1 \
        -molecular 0 \
        -outfile raw_feature \
        -rowatm "`echo ${rowatms}`" \
        -selatm "`echo ${selatms}`"\
        -mol2 cavity_decoy.mol2 complex_cavity_decoy.pdb #> /dev/null
    
    # process features
    echo "processing features.."
    paste <(awk -v model=${cavityID} -v x=${cavx} -v y=${cavy} -v z=${cavz} -v rna=${rna} '{print rna, model, sqrt(($2-x)*($2-x)+($3-y)*($3-y)+($4-z)*($4-z))}' ${cavityFile}) <(cat raw_feature.txt | grep -v 'atmid') | sed 's/,/./g' >> ${featureFile}
    
    # compute ligandability
    echo "now estimating ligandability.."
    python ${predictPy} ${model} ${featureFile} ${cavityFile} ${outFile} 1 #&> /dev/null
    cat ${outFile}
    
    echo "starting to prune grids"
    python ${CAVITYMINER}/py/grid_pruner.py  -c ${initialPDB} -s ${outFile} -o receptor
    
fi
