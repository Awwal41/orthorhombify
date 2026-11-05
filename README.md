# Orthorhombify

Turn any crystal (**monoclinic**, **hexagonal**, **rhombohedral** or **triclinic**) into a
near-**orthorhombic** supercell (all angles close to 90 degrees), see which **lattice
directions** the new axes follow, write **LAMMPS** files with starting charges for a
ReaxFF + QEq equilibration, and go **back to the original cell** afterwards.

- No strain: the new box is built only from whole lattice vectors.
- Works on any setting (`P 21/c`, `P 21/n`, `P 1 1 21/b`, `P 63/mmc`, `R -3`, `P -1`, ...).
- Hexagonal and rhombohedral cells become exactly orthorhombic.
- Every result says how its axes and planes are made of the old ones (`a' = a + c`, ...).
- Three terminal commands and a Python API. Needs only NumPy and ASE.

| Cell shape | What you get |
|------------|--------------|
| cubic, tetragonal, orthorhombic | already orthogonal, nothing to do |
| monoclinic | best box around the unique axis (exact 2D search) |
| hexagonal | exact orthohexagonal box, 2 cells |
| rhombohedral | exact box, 6 cells |
| triclinic | best box from a 3D search (closer to 90 with more cells) |

Documentation: **https://awwal41.github.io/orthorhombify/**

---

## Install

```bash
pip install git+https://github.com/Awwal41/orthorhombify.git
```

No `git`? `pip install https://github.com/Awwal41/orthorhombify/archive/refs/heads/main.zip`.
More options (venv, conda, HPC, from source) in the [install guide](https://awwal41.github.io/orthorhombify/install.html).

## Use from the terminal

**1. Make the orthorhombic supercell**

```bash
orthorhombify my_crystal.cif --max-cells 8 --top 2 --outdir output
```

You get `output/my_crystal_orth_1.cif`, `.POSCAR` and `.json` (best), `_orth_2` (runner-up).
The terminal (and the `.json`) shows the lattice directions, e.g. for beta-HMX:

```text
new a' = a + c               along old [1 0 1]         55.05 deg from old a
new b' = b                   along old [0 1 0]          0.00 deg from old b
new c' = -4a + 3c            along old [-4 0 3]        42.99 deg from old c
old a  = (3/7)a' - (1/7)c'   along new [3 0 -1]
old c  = (4/7)a' + (1/7)c'   along new [4 0 1]
new (0 0 1) = old (-1 0 1)       old (0 0 1) = new (1 0 3)
```

**2. Make the LAMMPS files**

```bash
orthorhombify-lammps output/my_crystal_orth_1.cif --outdir lammps --ffield ffield.reax
```

You get `lammps/my_crystal_orth_1/my_crystal_orth_1.data` and `in.my_crystal_orth_1`.

**3. Run LAMMPS** (put your ReaxFF force-field file in that folder first)

```bash
cd lammps/my_crystal_orth_1
lmp -in in.my_crystal_orth_1
```

**4. Back to the original cell** (from the starting folder)

```bash
orthorhombify-transform lammps/my_crystal_orth_1/my_crystal_orth_1_equil.data --inverse output/my_crystal_orth_1.json
```

The same command goes from any cell to any other with `--matrix`, e.g. orthohexagonal
to hexagonal: `--matrix "1 0 0; -1/2 1/2 0; 0 0 1"`.

## Use from Python

```python
from ase.io import read
import orthorhombify as orth

atoms = read("my_crystal.cif")               # or read(orth.example_path()) for beta-HMX
best = orth.find_supercells(atoms, max_cells=8, top=1)[0]
supercell = best.build(atoms)
print(best.mapping.summary(atoms.cell))      # lattice directions old <-> new
orth.write_lammps_inputs(supercell, "lammps/my_crystal", "my_crystal", ffield="ffield.reax")
```

See [Using it from Python](https://awwal41.github.io/orthorhombify/python-api.html).

## Starting charges

| Element | Charge (e) |
|---------|-----------|
| O | -0.3668 |
| N | +0.3668 |
| C | -0.4320 |
| H | +0.2160 |

If these do not add up to zero for your crystal, every charge is shifted by the same
small amount so the box is neutral. QEq then updates the charges during the run.
Change a value with `--charge O=-0.35`, or keep the raw values with `--no-neutralize`.

## Example results

| Crystal | Input angle | Result | Unit cells | Atoms |
|---------|-------------|--------|-----------|-------|
| beta-HMX | beta = 102.06 | beta = 89.999 | 7 | 392 |
| MDNT | gamma = 116.66 | gamma = 89.90 | 7 | 840 |
| TNT (monoclinic) | gamma = 111.01 | gamma = 90.13 | 4 | 672 |
| FOX-7 | beta = 90.065 | already the best cell | 1 | 56 |

More in the [docs](https://awwal41.github.io/orthorhombify/examples.html).

## Development

```bash
git clone https://github.com/Awwal41/orthorhombify.git
cd orthorhombify
pip install -e ".[test]"
pytest
```

## License

MIT. See [LICENSE](LICENSE).
