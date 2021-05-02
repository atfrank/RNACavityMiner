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

def prune_grid(rna, score_file, outname, quantile = 0.99, sasa_cutoff = 20.0):
    # make sure all atoms within an object occlude one another
    cmd.flag("ignore", "none")

    # use solvent-accessible surface with high sampling density
    cmd.set('dot_solvent', 1)
    cmd.set('dot_density', 3)
    cmd.set('solvent_radius', 2.0)

    k = 1
    df = pd.read_csv(score_file, header = 0, sep = ",")
    means = df[['pred_MLP','pred_XGB','pred_RF','pred_LR','pred_Extra']].apply(np.mean, 'columns')
    cutoff = np.quantile(a=means, q=quantile)
    
    for pred, pos in zip(df[['pred_MLP','pred_XGB','pred_RF','pred_LR','pred_Extra']].values, df[['x','y','z']].values):
        if np.mean(pred) > cutoff:
            # create tmp complex
            cmd.pseudoatom("tmpPoint3", hetatm = 1, name="C", resn = "UNK", pos= [pos[0], pos[1], pos[2]])
            cmd.create("complextmp", "%s tmpPoint3"%rna)
            sasa = cmd.get_area('resn UNK and not polymer and complextmp')
            cmd.delete("complextmp tmpPoint3")
            # remove really highly exposed points
            if sasa < sasa_cutoff:
                cmd.pseudoatom("tmpPoint", hetatm = 1, b = np.mean(pred), q = sasa, name="C", resn = "UNK", resi=k, chain="ZZ", pos= [pos[0], pos[1], pos[2]])
                print(pred, pos, np.mean(pred), cutoff, sasa)
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
    prune_grid(rna = "receptor", score_file = a.score_in, outname = a.coord_out, quantile = 0.95, sasa_cutoff = 10.0)

