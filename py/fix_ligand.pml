import sys
import inspect
from glob import glob
from pymol import cmd
import numpy
import os
import subprocess

cmd.load("lig.pdb")
cmd.alter("lig", "resn = 'UNK'")
cmd.alter("lig", "resi = '1'")
cmd.alter("lig", "type = 'HETATM'")
cmd.save("lig.pdb")
