import numpy as np
from bionmr_utils.md import *
import sys


# setup parameters.
path_to_traj = "../.."
first_dat_file = 51
last_dat_file = 800
stride = 100
fout_1 = 'H4-1_contact_map_Ca.txt'
fout_2 = 'H4-2_contact_map_Ca.txt'

residues_of_interest_1 = set(sorted(list(range(136, 160))))
residues_of_interest_2 = set(sorted(list(range(623, 647))))

# read trajectory
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/5_run/run00001.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))

# create predicates
rid_ca_1 = rId.is_in(residues_of_interest_1) & (aName == "CA")
rid_ca_2 = rId.is_in(residues_of_interest_2) & (aName == "CA")

# create vatiables for accumulation of data
n_residues = len(residues_of_interest_1)
h4_1_map = np.zeros((n_residues, n_residues))
h4_2_map = np.zeros((n_residues, n_residues))

# run through trajectory and calculate vectors
print("Processing frames...")
for frame in traj[::stride]:
    sys.stdout.write("Frame %d of %d\r" % (frame.index/stride + 1, traj.size/stride))
    if frame.index == 0:
        h4_1_ca = frame.asAtoms.filter(rid_ca_1)
        h4_2_ca = frame.asAtoms.filter(rid_ca_2)

    for i in range(n_residues):
        for j in range(i, n_residues):
            h4_1_map[i][j] += distance(h4_1_ca[i].r, h4_1_ca[j].r)
            h4_2_map[i][j] += distance(h4_2_ca[i].r, h4_2_ca[j].r)
            h4_1_map[j][i] = h4_1_map[i][j]
            h4_2_map[j][i] = h4_2_map[i][j]

np.savetxt(fout_1, h4_1_map/(traj.size/stride), fmt='%.2f')
np.savetxt(fout_2, h4_2_map/(traj.size/stride), fmt='%.2f')
