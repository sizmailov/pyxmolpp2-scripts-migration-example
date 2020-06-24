from bionmr_utils.md import *
from glob import glob
from subprocess import call


# This script process pdb files produced by pyrun scripts.
# Not protein chains are removed (only 8 first chains are retained).
# PDB files are formatted for processing with shiftx2:
#   - residues are renamed
#   - added two columns with useless data


path_to_pdb = 'pdb_shiftx2'

rename_map = {('HID', 'HIS'),
              ('HIE', 'HIS'),
              ('HIP', 'HIS'),
              ('LYN', 'LYS'),
              ('GLH', 'GLU'),
              ('CYM', 'CYS')}

files = sorted(glob(path_to_pdb + '/run?????.pdb'))

for file in files:
    print('Processing ', file)
    # load PDB-file, filter protein
    frame = PdbFile(file).frames()[0]
    
    for residue in frame.residues:
        residue.name = rename_map.get(residue.name, residue.name)

    # Note: occupancy/b_factors are not (yet) supported
    # 
    # for atom in frame.atoms:
    #     atom.occpancy = 1.0
    #     atom.b_factor = 0.0
    # 
    # frame.molecules[:8].to_pdb(file, write_occupancy=True, write_b_factors=True)

    frame.molecules[:8].to_pdb(file)
    # add 2 missing columns with useless data
    with open(file, 'r') as f:
        data = f.readlines()
    with open(file, 'w') as f:
        for line in data:
            if len(line) > 10:
                f.write(line[:54] + '  1.00  0.00' + line[65:])
