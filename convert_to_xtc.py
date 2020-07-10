import os
from typing import Sequence
from pyxmolpp2 import *
from pyxmolpp2.pipe import ScaleUnitCell

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable


def to_xtc(traj: Sequence[Frame], output_dir, precision=1000, basename_format="run{:05d}", first_filename_number=1,
           frames_per_file=1000,
           allow_incomplete_files=False, ):
    if not allow_incomplete_files:
        assert len(traj) % frames_per_file == 0

    n_files = (len(traj) + frames_per_file - 1) // frames_per_file
    for i in tqdm(range(0, n_files)):
        xtc_filename = os.path.join(output_dir, basename_format.format(first_filename_number + i) + ".xtc")
        xtc_writer = XtcWriter(xtc_filename, precision=precision)
        for frame in traj[i:i + frames_per_file]:
            xtc_writer.write(frame)
        del xtc_writer
        del frame


path = os.path.expanduser("~/GB1/trj/30R1")
ref = PdbFile(os.path.join(path, "ref/box.pdb")).frames()[0]
traj = Trajectory(ref)

for i in range(1, 4):
    traj.extend(TrjtoolDatFile(os.path.join(path, f"5_run/run{i:05d}.dat")))

print(f"Trajectory: {traj.n_frames} frames, {traj.n_atoms}")

to_xtc(traj | ScaleUnitCell(os.path.join(path, "5_run/summary.VOLUME"), colnum=0, max_rows=traj.size),
       output_dir=os.path.join(path, "5_run"))
