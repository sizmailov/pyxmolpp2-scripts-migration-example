from bionmr_utils.md import *
from get_secondary_structure_residues import *
from tqdm import tqdm
import numpy as np
import json


#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
reference_pdb = pars["reference_pdb"]
rmsd_fnout = pars["rmsd_fnout"]
len_inner_turn = pars["len_inner_turn"]


#
#  set atom selections
#
protein_and_dna_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
probe = ((cName == "I") & (aName == "P"))
dna_align_pred = aName.is_in({"N1", "N9"})


#
#  load trajectory
#
path_to_traj = "../.."
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/1_build/ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" %
      (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))


#
#  get secondary structure information:
#  get residue ids for secondary structure from DSSP using Biopython and PDB
#
protein_chains = "ABCDEFGH"
ss_residues_by_chain = []
incr = 0
for c in protein_chains:
    ss = get_secondary_structure_residues(c, pdb_code="1KX5")
    ss_residues_by_chain.append(set([(x+incr) for x in ss]))
    incr += len(ref.asResidues.filter(cName == c))
# print("Secondary structure residues by chain")
# print(ss_residues_by_chain)
ss_residues = [r for chain_r in ss_residues_by_chain for r in chain_r]
# print("Secondary structure residues combined")
# print(ss_residues)
ss_residues = set(ss_residues)


#
#  get residue ids for inner and outer DNA turns
#
dna = ref.asAtoms.filter(aName == "N1")
dna_resids = []
for at in dna:
    dna_resids.append(int(at.rId.serial))
n_dna_resids = len(dna_resids)
n_dna_bp = int(n_dna_resids / 2)
idx_01 = int(dna_resids[0]-1 + n_dna_bp // 2)
idx_02 = int(dna_resids[0]-1 + n_dna_bp + n_dna_bp // 2)
inner_turn = set(list(range(idx_01-len_inner_turn, idx_01+len_inner_turn+1, 1)) +
                 list(range(idx_02-len_inner_turn, idx_02+len_inner_turn+1, 1)))

outer_turn = set(list(range(dna_resids[0], idx_01-len_inner_turn, 1)) +
                 list(range(idx_01+len_inner_turn+1, dna_resids[0]+n_dna_bp+1, 1)) +
                 list(range(dna_resids[0]+n_dna_bp+1, idx_02-len_inner_turn, 1)) +
                 list(range(idx_02+len_inner_turn+1, dna_resids[-1]+1, 1)))
print("Nucleotides:", n_dna_resids)
print("Base pairs:", n_dna_bp)
print("Middle index I:", idx_01)
print("Middle index J:", idx_02)
print("Inner turn:", idx_01-len_inner_turn, "-", idx_01+len_inner_turn, ",",
                     idx_02-len_inner_turn, "-", idx_02+len_inner_turn)
print("Outer turn:", dna_resids[0], "-", idx_01-len_inner_turn-1, ",",
                     idx_01+len_inner_turn+1, "-", dna_resids[0]+n_dna_bp, ",",
                     dna_resids[0]+n_dna_bp+1, "-", idx_02-len_inner_turn-1, ",",
                     idx_02+len_inner_turn+1, "-", dna_resids[-1])


#
#  reference atoms for alignment
#  use 3LZ0 structure with reconstructed tails and propka protonation
#
reference = PdbFile(reference_pdb).get_frame()
ref_ats = reference.asAtoms
ref_align_ca_ss = ref_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))
ref_align_dna = ref_ats.filter(dna_align_pred)
ref_align_union = ref_ats.filter(((aName == "CA") & (rId.is_in(ss_residues))) | dna_align_pred)
ref_align_inner_turn = ref_ats.filter(dna_align_pred).filter(rId.is_in(inner_turn))
ref_align_outer_turn = ref_ats.filter(dna_align_pred).filter(rId.is_in(outer_turn))


#
#  create data variables
#
rmsd_ref_align_ca_ss = np.zeros(int(traj.size / stride))
rmsd_ref_align_dna = np.zeros(int(traj.size / stride))
rmsd_ref_align_union = np.zeros(int(traj.size / stride))
rmsd_ref_align_inner_turn = np.zeros(int(traj.size / stride))
rmsd_ref_align_outer_turn = np.zeros(int(traj.size / stride))


#
#  Handle PBC:
#
def get_XST(path_to_traj):
    time, volume = np.genfromtxt(path_to_traj + "/5_run/summary.VOLUME", unpack=True)

    inpcrd = path_to_traj + "/1_build/box.inpcrd"

    with open(inpcrd) as f:
        unit_cell_line = f.readlines()[-1]

    unit_cell = UnitCell.from_rst7_line(unit_cell_line)

    result = []
    for t, frame_volume in zip(time, volume):
        cell = UnitCell(unit_cell)
        cell.scale_volume_to(frame_volume)
        result.append((t, cell))

    return result


# get PBC from XST file and inpcrd file
time_cells = get_XST(path_to_traj)
assert len(time_cells) >= traj.size, "PBC data is not enough"
inpcrd = path_to_traj + "/1_build/box.inpcrd"


#
#  run through trajectory
#
for frame in tqdm(traj[::stride], disable=False) | BuildCell("summary.VOLUME"):

    time, cell = time_cells[frame.index]

    if frame.index == 0:
        # create variables for selections from atoms of frame
        ref_pbc_ats = ref.atoms
        frame_ats = frame.atoms

        ref_probe = ref_pbc_ats.filter(probe)
        frame_probe = frame_ats.filter(probe)

        frame_align_ca_ss = frame_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))
        frame_align_dna = frame_ats.filter(dna_align_pred)
        frame_align_union = frame_ats.filter(((aName == "CA") & (rId.is_in(ss_residues))) | dna_align_pred)
        frame_align_inner_turn = frame_ats.filter(dna_align_pred).filter(rId.is_in(inner_turn))
        frame_align_outer_turn = frame_ats.filter(dna_align_pred).filter(rId.is_in(outer_turn))

        ats_ref = []
        ats_frame = []
        for cid in sorted(protein_and_dna_chains):
            ats_ref.append(ref_pbc_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))

    # unwrapping
    # align all reference atoms by one of DNA chains
    ref_pbc_ats.apply(ref_probe.alignment_to(frame_probe))
 
    for i, cid in enumerate(sorted(protein_and_dna_chains)):
        img = cell.closest_image_to(ats_ref[i].coords.mean(), ats_frame[i].coords.mean())
        ats_frame[i].apply(Translation(img.shift))

    # align nucleosome in frame by reference 3LZ0 structure using sec.str. CA atoms
    frame_ats.apply(frame_align_ca_ss.alignment_to(ref_align_ca_ss))

    # calculate RMSD
    k = int(frame.index // stride)
    rmsd_ref_align_ca_ss[k]      = ref_align_ca_ss.rmsd(frame_align_ca_ss)
    rmsd_ref_align_dna[k]        = ref_align_dna.rmsd(frame_align_dna)
    rmsd_ref_align_union[k]      = ref_align_union.rmsd(frame_align_union)
    rmsd_ref_align_inner_turn[k] = ref_align_inner_turn.rmsd(frame_align_inner_turn)
    rmsd_ref_align_outer_turn[k] = ref_align_outer_turn.rmsd(frame_align_outer_turn)


#
#  write RMSD to file
#
header = ["time [ps]",
          r"sec.str. C$\rm\alpha$ [A]",
          "DNA [A]",
          r"DNA & sec.str. C$\rm\alpha$ [A]",
          "DNA inner turn [A]",
          "DNA outer turn [A]"]
time = np.linspace((first_dat_file - 1) * 1000 + stride, last_dat_file * 1000, int(traj.size / stride))
rmsd = np.vstack((time,
                  rmsd_ref_align_ca_ss,
                  rmsd_ref_align_dna,
                  rmsd_ref_align_union,
                  rmsd_ref_align_inner_turn,
                  rmsd_ref_align_outer_turn)).T
np.savetxt(
    fname=rmsd_fnout,
    X=rmsd,
    fmt="%.5f",
    header=",".join(header),
    delimiter=","
)
