import sys
import inspect
from glob import glob
from pymol import cmd
import numpy
import os
import subprocess

cmd.load("receptor.mol2")
cmd.alter("receptor and resn C+RC+rC+C3+C5", "resn = 'CYT'")
cmd.alter("receptor and resn A+RA+rA+A3+A5", "resn = 'ADE'")
cmd.alter("receptor and resn U+RU+rU+U3+U5", "resn = 'URA'")
cmd.alter("receptor and resn G+RG+rG+G3+G5", "resn = 'GUA'")
cmd.alter("receptor", "type = 'ATOM'")
cmd.save("receptor.pdb")

center = cmd.get_position()

print("MY CENTER X: %s"%(int(center[0]))) 
print("MY CENTER Y: %s"%(int(center[1])))
print("MY CENTER Z: %s"%(int(center[2])))