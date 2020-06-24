from bionmr_utils.md import *
from glob import glob

# This scripts copy files for unwrapped nucleosome to pdb_shifts,
# then it creates three additional files:
#   1. Only protein file (for proton

path_to_pdb = 'pdb_shifts'

files = sorted(glob(path_to_pdb + '/run?????.pdb'))

for file in files:
    
    # load PDB-file
    frame = PdbFile(file).get_frame()

    # save protein only chains
    frame_p = Frame(0)
    protein = frame.asChains[:8]
    for chain in protein:
        frame_p.emplace(chain)
    frame_p.to_pdb(file[:-4]+'_p.pdb')

    # save all chains with one cName for protein + nucleic
    f_old = open(file, 'r')
    with open(file[:-4]+'_cn.pdb', 'w') as f:
        for line in f_old:
            if (line.split()[0] == 'TER') | ('CRYST' in line.split()[0]):
                continue
            f.write(line[0:21] + 'X' + line[22:])
    f.close()
    f_old.close()

    # save all chains with one cName for protein only
    f_old = open(file[:-4] + '_p.pdb', 'r')
    with open(file[:-4]+'_p_cn.pdb', 'w') as f:
        for line in f_old:
            if (line.split()[0] == 'TER') | ('CRYST' in line.split()[0]):
                continue
            f.write(line[0:21] + 'X' + line[22:])
    f.close()
    f_old.close()
