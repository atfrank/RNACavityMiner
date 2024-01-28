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

def make_profile(profile_name, obj="resn UNK and name C", rna_obj = 'complex', cutoff = 5):
    """ 
    Generates a ligandability profile, i.e., ligandability as a function of residue position    
    """
    atoms1 = cmd.get_model(obj)
    profile = {"%s"%(at.resi): 0 for at in cmd.get_model("name C1' and polymer.nucleic and %s"%(rna_obj)).atom}
    print(profile)
    for at in atoms1.atom:
        sele = "br. polymer.nucleic and %s within %s of (resn %s and resi %s and name %s)"%(rna_obj, cutoff+2, at.resn, at.resi, at.name)       
        atoms2 = cmd.get_model(sele)
        ligandability = at.b
        sasa = at.q
        for bt in atoms2.atom:          
            resi = bt.resi
            distance = np.sqrt((at.coord[0]-bt.coord[0])**2 +  (at.coord[1]-bt.coord[1])**2 + (at.coord[2]-bt.coord[2])**2)
            sij = ligandability*np.exp(-(distance/cutoff)**2)
            profile[resi] += sij
            #print(at.q, sele)
            #if distance < cutoff: print(distance, sele)
    print(profile)  
    df = pd.DataFrame.from_records([{"residue": key, "ligandability": profile[key]} for key in profile.keys()])
    df.to_csv(profile_name, index=False)
    return(profile)

def remove_isolated(obj="resn UNK and name C", cutoff = 5):
    atoms = cmd.get_model(obj)
    for at in atoms.atom:
        distances = []
        for bt in atoms.atom:
            distance = np.sqrt((at.coord[0]-bt.coord[0])**2 +  (at.coord[1]-bt.coord[1])**2 + (at.coord[2]-bt.coord[2])**2)
            if distance != 0.0 and distance < cutoff: distances.append(distance)
        if len(distances) < 4:
            sele = "resn %s and resi %s and name %s"%(at.resn, at.resi, at.name)
            print(sele)
            cmd.remove(sele)
            atoms = cmd.get_model(obj)

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
    coor = "cavity_pruned_grid.sd"
    cmd.save(coor, "tmpPoint")

    # remove isolated
    remove_isolated()
        
    # write out grid file
    coor = "%s_pruned_grid_clusters.xyz"%(outname)
    xyz = cmd.get_coords('tmpPoint', 1)
    df = pd.DataFrame.from_records(xyz)
    df.insert(0, "element", "C")
    df.to_csv(coor, index = False, header = False, sep = " ")    
    
    # generate ligandability profile
    cmd.create("complex", "%s tmpPoint"%rna)
    profile = make_profile("%s_profile.csv"%(outname), rna_obj = 'complex')
    cmd.alter("(polymer.nucleic)", 'b=0.0')
    [cmd.alter("(polymer.nucleic and resi %s)"%(key), 'b=%6.3f'%(profile[key])) for key in profile.keys()]
    
    # write out complex
    coor = "%s_pruned_grid_clusters.pdb"%(outname)
    cmd.save(coor, "complex")
    coor = "cavity_pruned_grid.sd"
    cmd.save(coor, "tmpPoint")

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--coord_in", help="path to input coordinate file")
    parser.add_argument("-s","--score_in", help="path to input ligandability score file")
    parser.add_argument("-o","--coord_out", help="path to output (prefix)")
    
    # parse command line
    a = parser.parse_args()
    cmd.load(a.coord_in, "receptor")
    prune_grid(rna = "receptor", score_file = a.score_in, outname = a.coord_out, quantile = 0.80, sasa_cutoff = 10.0)

