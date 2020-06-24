import numpy as np
from math import ceil
from pyxmolpp2 import *
import os
from tqdm import tqdm
import sys
# from correlation_functions import cor

# setup parameters.
path_to_traj = "/home/seva/chromatin/5_solution_Widom_601/5_ff14SB_CUFIX_TIP4PEW_small_box/"
first_dat_file = 1
last_dat_file = 430
n_steps = last_dat_file - first_dat_file + 1
stride = 1
cut_autocorr_function = n_steps  # ns
HN_mask = "H"
align_dna = True
fft_acf = True
scaling = 1.005
ratio = 0.8
remove_first_point = True

residues_of_interest = set(list(range(136, 238)) + list(range(623, 725)))

# read trajectory
ref = PdbFile(path_to_traj + "/5_run/run00001.pdb").frames()[0]
traj = Trajectory(ref)
for i in tqdm(range(first_dat_file, last_dat_file+1), desc="reading .dat headers"):
    traj.extend(TrjtoolDatFile(path_to_traj + f"/5_run/run{i:05d}.dat"))

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (
    traj.n_frames(), ref.molecules.size, ref.residues.size, ref.atoms.size))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))

# create folders and open files
resids = []
resnames = []
for at in ref.atoms.filter((aName == HN_mask) & (rId.is_in(residues_of_interest))):
    resids.append(at.residue.id.serial)
    resnames.append(at.residue.name)
resi = tuple(zip(resids, resnames))
print("Autocorrelation functions will be calculated for following residues:")
print(resi)

if not os.path.exists("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)):
    os.makedirs("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file))
# vectors - list of VectXYZ
vectors = np.zeros((len(resids), traj.n_frames() // stride, 3), dtype=np.float64)
# vectors = []
skip_frames = []

first_time = True

# run through trajectory and calculate vectors
N, H, s1, s2, moved_atoms = None, None, None, None, None

with tqdm(traj[::stride], total=traj.n_frames() // stride, desc="Processing frames") as pbar:
    for frame_id, frame in enumerate(pbar):  # type: Trajectory.Frame
        # sys.stdout.write("Frame %d of %d\r" % (frame_id + 1, traj.n_frames() / stride))

        if first_time:
            s1 = frame.atoms.filter(aName == "P").coords
            s2 = ref.atoms.filter(aName == "P").coords
            moved_atoms = frame.atoms.filter(aName.is_in({"N",HN_mask,"P" })).coords
            assert len(s1) == len(s2)
            Nr = AtomSelection([ r["N"] for r in frame.residues.filter(rId.is_in(set(resids)))]).coords
            Hr = AtomSelection([ r[HN_mask] for r in frame.residues.filter(rId.is_in(set(resids)))]).coords
            assert Nr.size == len(resids)
            assert Nr.size == Hr.size
            delta = Hr.values - Nr.values
            vectors[:,0,:] = delta
            first_time = False
        else:
            moved_atoms.apply(s1.alignment_to(s2))

            delta = Hr.values - Nr.values
            vectors[:,frame_id,:] = delta

            if np.linalg.norm(delta, axis=1).max() < 0.5:
                skip_frames += [frame.index]
                pbar.set_postfix(skipped_frames=len(skip_frames))

print()

if len(skip_frames) > 0:
    fft_acf = False
    np.savetxt('skip_frames.txt', np.array(list(skip_frames)), delimiter=',', header='')

# calculate autocorrelation functions
steps = np.arange(int(cut_autocorr_function * 1000 / stride))
grid = []
nlim = ratio * len(vectors[0])
if not remove_first_point:
    grid.append(0)
tau = 1.0
while tau <= nlim:
    grid.append(int(tau))
    tau = ceil(tau * scaling)

print("Calculating autocorrelation functions...")
for i, (rid, rname) in enumerate(tqdm(zip(resids, resnames), "Calc autocorrelations", total=len(resids))):
    output_filename = f"cor_NH_{first_dat_file}-{last_dat_file}_diluted/{rid:04d}_{rname}.cor"
    if fft_acf:
        ac = calc_autocorr_order_2(vectors[i]) 
        np.savetxt(output_filename, np.vstack((steps[grid], ac[grid])).T, fmt="%14.6e")
    else:
        ac = cor(vectors[i], grid, skip_frames)
        np.savetxt(output_filename, np.vstack((steps[grid], ac)).T, fmt="%14.6e")

sys.stdout.write("\n")
print("Done!")
