#!/bin/bash
module load anaconda
module load gcc

source activate Jingru
module load openbabel
module load rdock

for i in {1..300}
do
	/home/afrankz/local_software/repo/CavityPoser/src/cavityPoser.sh /home/afrankz/demos/pocketPosers/cavities/IRES/${i}/receptor.pdb
	mv receptor receptor_${i}
done
