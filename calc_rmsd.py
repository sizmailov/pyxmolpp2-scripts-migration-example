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
    time, V = np.genfromtxt(path_to_traj + "/5_run/summary.VOLUME", unpack=True)

    inpcrd = path_to_traj + "/1_build/box.inpcrd"

    with open(inpcrd) as f:
        content = f.readlines()

    content = content[-1].split()[:6]

    a, b, c, alpha, beta, gamma = [float(x) for x in content]

    # https://en.wikipedia.org/wiki/Parallelepiped
    V0 = a * b * c * np.sqrt(
        1.0 + 2.0 * np.cos(alpha / 180 * np.pi) * np.cos(beta / 180 * np.pi) * np.cos(gamma / 180 * np.pi) - np.cos(
            alpha / 180 * np.pi) ** 2 - np.cos(beta / 180 * np.pi) ** 2 - np.cos(gamma / 180 * np.pi) ** 2)

    # XST file format
    # http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2003-2004/0234.html

    m_v1 = np.array([a, 0.0, 0.0])
    m_v2 = np.array([b * np.cos(gamma / 180 * np.pi), b * np.sin(gamma / 180 * np.pi), 0])
    m_v3 = np.zeros(3)

    m_v3[0] = c * np.cos(beta / 180 * np.pi)
    m_v3[1] = b / m_v2[1] * c * np.cos(alpha / 180 * np.pi) - m_v2[0] / m_v2[1] * m_v3[0]
    m_v3[2] = np.sqrt(c * c - m_v3[0] * m_v3[0] - m_v3[1] * m_v3[1])

    xst = np.zeros((len(V), 13))
    origin = np.zeros(3)

    for i, f_V in enumerate(V):
        scaling_factor = np.cbrt(f_V / V0)
        # print(scaling_factor)
        xst[i, 0] = time[i]
        xst[i, 1:10] = np.array(
            [m_v1[0], m_v1[1], m_v1[2], m_v2[0], m_v2[1], m_v2[2], m_v3[0], m_v3[1], m_v3[2]]) * scaling_factor
        xst[i, 10:13] = [origin[0], origin[1], origin[2]]

    return xst


# get PBC from XST file and inpcrd file
pbc = get_XST(path_to_traj)
v1 = VectorXYZ.from_numpy(pbc[:, 1:4])
v2 = VectorXYZ.from_numpy(pbc[:, 4:7])
v3 = VectorXYZ.from_numpy(pbc[:, 7:10])
assert len(v1) >= traj.size, "PBC data is not enough"
scaling_factors = [v1[i].len()/v1[0].len() for i in range(len(v1))]
lat_vec = LatticeVectors(v1[0], v2[0], v3[0])
shift_finder = BestShiftFinder(lat_vec)
inpcrd = path_to_traj + "/1_build/box.inpcrd"
with open(inpcrd) as f:
    content = f.readlines()
content = content[-1].split()[:6]
a, b, c, alpha, beta, gamma = [float(x) for x in content]


#
#  run through trajectory
#
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        # create variables for selections from atoms of frame
        ref_pbc_ats = ref.asAtoms
        frame_ats = frame.asAtoms

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
    ref_pbc_ats.transform(ref_probe.alignment_to(frame_probe))

    # cycle through chains and check first atom displacement
    shift_finder.scale_lattice_by(scaling_factors[int(frame.index / stride)])
    for i, cid in enumerate(sorted(protein_and_dna_chains)):
        dif, shift = shift_finder.find_best_shift(ats_ref[i].geom_center(),
                                                  ats_frame[i].geom_center())
        ats_frame[i].transform(Translation3d(shift))
    shift_finder.scale_lattice_by(1.0 / scaling_factors[int(frame.index / stride)])

    # align nucleosome in frame by reference 3LZ0 structure using sec.str. CA atoms
    frame_ats.transform(frame_align_ca_ss.alignment_to(ref_align_ca_ss))

    # calculate RMSD
    rmsd_ref_align_ca_ss[int(frame.index / stride)] = calc_rmsd(ref_align_ca_ss.toCoords, frame_align_ca_ss.toCoords)
    rmsd_ref_align_dna[int(frame.index / stride)] = calc_rmsd(ref_align_dna.toCoords, frame_align_dna.toCoords)
    rmsd_ref_align_union[int(frame.index / stride)] = calc_rmsd(ref_align_union.toCoords, frame_align_union.toCoords)
    rmsd_ref_align_inner_turn[int(frame.index / stride)] = calc_rmsd(ref_align_inner_turn.toCoords, frame_align_inner_turn.toCoords)
    rmsd_ref_align_outer_turn[int(frame.index / stride)] = calc_rmsd(ref_align_outer_turn.toCoords, frame_align_outer_turn.toCoords)


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
