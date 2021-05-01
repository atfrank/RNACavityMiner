import os
import sys
import argparse
import pandas as pd
import numpy as np
from glob import glob
from pymol import cmd, CmdException
import itertools as it

def logic(index, skip):
        """ check modulus """
        if index % skip == 0:
             return True
        return False

def check_file(file):
        """ check whether file is present """
        try:
            if os.path.getsize(file) > 0:
                return True
            else:
                return False
        except:
            return False

def prune_grid(rna, score_file, outname, cutoff = 0.99):
    k = 1
    df = pd.read_csv(score_file, header = 0, sep = ",")
    means = df[['pred_MLP','pred_XGB','pred_RF','pred_LR','pred_Extra']].apply(np.mean, 'columns')
    cutoff = np.quantile(a=means, q=cutoff)
    
    for pred, pos in zip(df[['pred_MLP','pred_XGB','pred_RF','pred_LR','pred_Extra']].values, df[['x','y','z']].values):
        if np.mean(pred) > cutoff:
            print(pred, pos, np.mean(pred))
            cmd.pseudoatom("tmpPoint2", hetatm = 1, b = np.mean(pred), name="C", resn = "UNK", resi=k, chain="ZZ", pos= [pos[0], pos[1], pos[2]])
            k += 1

    # write out grid file
    coor = "%s_pruned_grid.xyz"%(outname)
    xyz = cmd.get_coords('tmpPoint', 1)
    df = pd.DataFrame.from_records(xyz)
    df.insert(0, "element", "C")
    df.to_csv(coor, index = False, header = False, sep = " ")
    
    # write out complex
    cmd.create("complex", "%s tmpPoint"%rna)
    coor = "%s_pruned_grid.pdb"%(outname)
    cmd.save(coor, "complex")
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--coord_in", help="path to input coordinate file")
    parser.add_argument("-s","--score_in", help="path to input ligandability score file")
    parser.add_argument("-o","--coord_out", help="path to output (prefix)")
    
    # parse command line
    a = parser.parse_args()
    cmd.load(a.coord_in, "receptor")
    prune_grid(rna = "receptor", score_file = a.score_in, outname = a.coord_out)

