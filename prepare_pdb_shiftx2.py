from bionmr_utils.md import *
from glob import glob
from subprocess import call


# This script process pdb files produced by pyrun scripts.
# Not protein chains are removed (only 8 first chains are retained).
# PDB files are formatted for processing with shiftx2:
#   - residues are renamed
#   - added two columns with useless data


path_to_pdb = 'pdb_shiftx2'

rename_map = [('HID', 'HIS'),
              ('HIE', 'HIS'),
              ('HIP', 'HIS'),
              ('LYN', 'LYS'),
              ('GLH', 'GLU'),
              ('CYM', 'CYS')]

files = sorted(glob(path_to_pdb + '/run?????.pdb'))

for file in files:
    print('Processing ', file)
    # load PDB-file, filter protein
    frame = PdbFile(file).get_frame()
    new_frame = Frame(0)
    for chain in frame.asChains[:8]:
        new_frame.emplace(chain)
    # save protein to PDB-file
    new_frame.to_pdb(file)

    # add 2 missing columns with useless data
    with open(file, 'r') as f:
        data = f.readlines()
    with open(file, 'w') as f:
        for line in data:
            if len(line) > 10:
                f.write(line[:54] + '  1.00  0.00' + line[65:])

    # rename residues from amber to canonical
    for n1, n2 in rename_map:
        call('sed -i "s/%s/%s/g" ' % (n1, n2) + file, shell=True)
