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

def generate_grid(rna, outname = "test", grid_size = 2, padding = 2.0, mindist = 2.5, maxdist = 5.0, debug = False):
    # goto working directory

    # setup 3D grid
    extent = cmd.get_extent(selection=rna)
    x = np.arange(extent[0][0]-padding, extent[1][0]+padding, grid_size)
    y = np.arange(extent[0][1]-padding, extent[1][1]+padding, grid_size)
    z = np.arange(extent[0][2]-padding, extent[1][2]+padding, grid_size)
    xx, yy, zz = np.meshgrid(x, y, z)
    nx, ny, nz = xx.shape[0], xx.shape[1], xx.shape[2]
    
    # place pseudoatoms along grid
    k=1
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                cmd.pseudoatom("tmpPoint", hetatm = 1, name="C", resn = "UNK", resi=k, chain="ZZ", pos=[float(xx[ix, iy, iz]), float(yy[ix, iy, iz]), float(zz[ix, iy, iz])])
                # prune
                cmd.remove("resn UNK and resi %s within %s of polymer"%(k, mindist))
                cmd.remove("resn UNK and resi %s beyond %s of polymer"%(k, maxdist))

    # write out grid file
    coor = "%s_grid.xyz"%(outname)
    xyz = cmd.get_coords('tmpPoint', 1)
    df = pd.DataFrame.from_records(xyz)
    df.insert(0, "element", "C")
    df.to_csv(coor, index = False, header = False, sep = " ")
    
    # write out complex
    cmd.create("complex", "%s tmpPoint"%rna)
    coor = "%s_grid.pdb"%(outname)
    cmd.save(coor, "complex")
    
    if debug:
        cmd.show("surface", "polymer")
        cmd.show("spheres", "resn UNK")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--coord_in", help="path to input coordinate file")
    parser.add_argument("-o","--coord_out", help="path to output (prefix)")
    
    # parse command line
    a = parser.parse_args()
    cmd.load(a.coord_in, "receptor")
    generate_grid(rna = "receptor", outname = a.coord_out, grid_size = 2, padding = 2.0, mindist = 2.5, maxdist = 5.0, debug = False)
