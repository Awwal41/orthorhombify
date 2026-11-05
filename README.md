# Orthorhombify

Turn a **monoclinic** or **triclinic** crystal into a near-**orthorhombic** supercell
(all angles close to 90 degrees), then write **LAMMPS** files with starting charges
for a ReaxFF + QEq equilibration.

- No strain: the new box is built only from whole lattice vectors.
- Works on any monoclinic setting (`P 21/c`, `P 21/n`, `P 1 1 21/b`, ...).
- Two small scripts. Needs only Python, NumPy and ASE.

Documentation: **https://awwal41.github.io/orthorhombify/**

---

## Install

```bash
git clone https://github.com/Awwal41/orthorhombify.git
cd orthorhombify
pip install -r requirements.txt
```

## Use

**1. Make the orthorhombic supercell**

```bash
python orthorhombify_v2.py my_crystal.cif --max-cells 8 --top 2 --outdir output
```

You get `output/my_crystal_orth_1.cif` and `.POSCAR` (best), `_orth_2` (runner-up).

**2. Make the LAMMPS files**

```bash
python make_lammps_inputs.py output/my_crystal_orth_1.cif --outdir lammps --ffield ffield.reax
```

You get `lammps/my_crystal_orth_1/my_crystal_orth_1.data` and `in.my_crystal_orth_1`.

**3. Run LAMMPS** (put your ReaxFF force-field file in that folder first)

```bash
cd lammps/my_crystal_orth_1
lmp -in in.my_crystal_orth_1
```

## Starting charges

| Element | Charge (e) |
|---------|-----------|
| O | -0.3668 |
| N | +0.3668 |
| C | -0.4320 |
| H | +0.2160 |

If these do not add up to zero for your crystal, the script shifts every charge by the
same small amount so the box is neutral. QEq then updates the charges during the run.
Change a value with `--charge O=-0.35`, or keep the raw values with `--no-neutralize`.

## Example results

| Crystal | Input angle | Result | Unit cells | Atoms |
|---------|-------------|--------|-----------|-------|
| beta-HMX | beta = 102.06 | beta = 89.999 | 7 | 392 |
| MDNT | gamma = 116.66 | gamma = 89.90 | 7 | 840 |
| TNT (monoclinic) | gamma = 111.01 | gamma = 90.13 | 4 | 672 |
| FOX-7 | beta = 90.065 | already the best cell | 1 | 56 |

More in the [docs](https://awwal41.github.io/orthorhombify/examples.html).

## License

MIT. See [LICENSE](LICENSE).
