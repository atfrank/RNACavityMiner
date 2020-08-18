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

# get command line arguments
initialPDB=$1
initialPDB=`realpath $initialPDB`
rna=`basename $initialPDB | sed 's/.pdb//g'`
model=${DIRMODEL}/models.pkl

# setup working directory
workingDIR="`pwd`/${rna}"
mkdir -p ${workingDIR}
cd ${workingDIR}
outFile=${workingDIR}/predicted_cavities.txt
featureFile=${workingDIR}/features.csv
cavityFile=${workingDIR}/cavity_centers.txt
cavityID="none"

rm -f *

# get pocket analysis files
babel -ipdb ${initialPDB} -omol2 receptor.mol2 &> /dev/null
pymol -cqr $pml | tee receptor.log
coorx=`grep "CENTER X" receptor.log | grep -v PyMOL | awk '{print $4}'`
coory=`grep "CENTER Y" receptor.log | grep -v PyMOL | awk '{print $4}'`
coorz=`grep "CENTER Z" receptor.log | grep -v PyMOL | awk '{print $4}'`

# get reference pocket
cp ${rdockParam} receptor.prm
echo "sed -i 's/XX/${coorx}/g' receptor.prm" | bash
echo "sed -i 's/YY/${coory}/g' receptor.prm" | bash
echo "sed -i 's/ZZ/${coorz}/g' receptor.prm" | bash

# geometerical detect cavities
rbcavity -W -d -r ${workingDIR}/receptor.prm &> cavities.out

# store cavity center
cavx=${coorx}
cavy=${coory}
cavz=${coorz}

# create blind cavity centers
cat cavities.out  | grep 'Center=' | sed 's/Center=/\nCenter=/g' | sed 's/Extent=/\nExtent=/g' | grep Center= | sed 's/(/( /g' | sed 's/)/ )/g' | sed 's/,/ /g' | awk '{print $2, $3, $4}' | tr ' ' '\n' | sed '3~3 s/$/\nMy Text/g' | tr '\n' ' ' | sed 's/My Text/\n/g' | sed '/^[[:space:]]*$/d' | awk '{printf "%s %4.3f %4.3f %4.3f\n", "C", $1,$2,$3}' > cavity_centers.txt
if [[ -s cavity_centers.txt ]]
then
	ncav=`wc -l cavity_centers.txt | awk '{print $1}'`
    echo "${ncav}         " > cavity_decoy.xyz
    echo "decoys" >> cavity_decoy.xyz
    cat cavity_centers.txt >> cavity_decoy.xyz
    babel -ixyz cavity_decoy.xyz -opdb cavity_decoy.pdb &> /dev/null
    sed -i 's/UNL/UNK/g' cavity_decoy.pdb
    sed -i 's/LIG/UNK/g' cavity_decoy.pdb
    babel -ipdb cavity_decoy.pdb -omol2 cavity_decoy.mol2 &> /dev/null

    # create complex native
    grep -v 'TER' receptor.pdb | grep -v 'END' > complex_cavity_decoy.pdb
    echo "TER" >> complex_cavity_decoy.pdb
    grep 'HETATM' cavity_decoy.pdb >> complex_cavity_decoy.pdb
    echo "END" >> complex_cavity_decoy.pdb

    # featurize decoy cavity
	${CAVITYMINER}/bin/featurize \
        -etaStartPow 1 \
        -numEta 4 \
        -cutoff ${cutoff} \
        -scalar 1 \
        -molecular 0 \
		-outfile raw_feature \
        -rowatm "`echo ${rowatms}`" \
        -selatm "`echo ${selatms}`"\
        -mol2 cavity_decoy.mol2 complex_cavity_decoy.pdb > /dev/null
    
    # process features
	paste <(awk -v model=${cavityID} -v x=${cavx} -v y=${cavy} -v z=${cavz} -v rna=${rna} '{print rna, model, sqrt(($2-x)*($2-x)+($3-y)*($3-y)+($4-z)*($4-z))}' cavity_centers.txt) <(cat raw_feature.txt | grep -v 'atmid') >> ${featureFile}
	
	# compute druggability
	python ${predictPy} ${model} ${featureFile} ${cavityFile} ${outFile} 1 &> /dev/null
	cat ${outFile}
	
fi
