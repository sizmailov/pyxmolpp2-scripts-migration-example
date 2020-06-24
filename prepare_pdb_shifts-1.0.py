from pyxmolpp2 import *
from glob import glob

# This scripts copy files for unwrapped nucleosome to pdb_shifts,
# then it creates three additional files:
#   1. Only protein file (for proton

path_to_pdb = 'pdb_shifts'

files = sorted(glob(path_to_pdb + '/run?????.pdb'))

for file in files:
    
    # load PDB-file
    frame = PdbFile(file).frames()[0]

    # save protein only chains
    frame.molecules[:8].to_pdb(file[:-4]+'_p.pdb')


    # save all chains with one cName for protein only
    for mol in frame.molecules[:8]:
        mol.name = "X"
    frame.to_pdb(file[:-4]+'_p_cn.pdb')


    # save all chains with one cName for protein + nucleic
    for mol in frame.molecules[8:]:
        mol.name = "X"
    frame.to_pdb(file[:-4]+'_cn.pdb')
