# Getting started

From a CIF file to a running LAMMPS equilibration in a few steps.

**You need:** Python 3.9 or newer, and a CIF file of your crystal. For the last step
you also need LAMMPS (built with the `REAXFF` package) and a ReaxFF force-field file.

Open a terminal:

- Windows: **PowerShell** (Start menu → type `powershell`)
- macOS: **Terminal**
- Linux: your usual terminal

---

## 1. Copy Orthorhombify

```bash
git clone https://github.com/Awwal41/orthorhombify.git
cd orthorhombify
```

Stay in this `orthorhombify` folder for every later command.

---

## 2. Install the two libraries

```bash
pip install -r requirements.txt
```

This installs **NumPy** and **ASE** (the Atomic Simulation Environment).

Check:

```bash
python -c "import ase; print(ase.__version__)"
```

You want a version number, not an error.

---

## 3. Make the orthorhombic supercell

Put your CIF in the folder, then:

```bash
python orthorhombify_v2.py my_crystal.cif --max-cells 8 --top 2 --outdir output
```

**Good sign:** the script says what it found, for example

```text
Detected MONOCLINIC, unique axis b (beta = 102.058 deg is the only non-90 angle).
Best possible angle deviation by supercell size N (x unit cell):
  N=1: 12.058  N=2: 6.647  N=3: 5.177  N=4: 0.833  ...  N=7: 0.001
[1] N=7 (392 atoms)  ...  angles = 90.000, 89.999, 90.000 deg
```

The `N=...` line shows how close to 90 degrees you can get for each box size. Pick the
trade-off you want (see [Orthorhombify](orthorhombify.md)).

You now have:

```text
output/my_crystal_orth_1.cif      best result
output/my_crystal_orth_1.POSCAR
output/my_crystal_orth_2.cif      runner-up
output/my_crystal_orth_2.POSCAR
```

Open the CIF in VESTA, Ovito or Mercury to look at it.

---

## 4. Make the LAMMPS files

```bash
python make_lammps_inputs.py output/my_crystal_orth_1.cif --outdir lammps --ffield ffield.reax
```

The script prints the net charge and fixes it to zero if needed:

```text
net charge with initial charges = -15.7472 e
neutralized: shifted every charge by +0.018747 e -> C = -0.4133, H = +0.2347, N = +0.3855, O = -0.3481
NOTE: ... shortest is below 20.0 A ... Consider --replicate 3 1 2
```

If you see the `NOTE`, run the same command again with the suggested `--replicate`.

You now have:

```text
lammps/my_crystal_orth_1/my_crystal_orth_1.data   atoms, box, charges
lammps/my_crystal_orth_1/in.my_crystal_orth_1     LAMMPS input script
```

---

## 5. Run LAMMPS

Copy your ReaxFF force-field file into `lammps/my_crystal_orth_1/` and name it
`ffield.reax` (or pass `--ffield yourfile` in step 4). Then:

```bash
cd lammps/my_crystal_orth_1
lmp -in in.my_crystal_orth_1
```

Your LAMMPS program may be called `lmp_serial`, `lmp_mpi` or `lmp`. In parallel:
`mpirun -np 8 lmp -in in.my_crystal_orth_1`.

The log prints the average C, H, N and O charge every 100 steps so you can watch QEq
update them. At the end you get `my_crystal_orth_1_equil.data` with the equilibrated
structure and charges.

---

## Many crystals at once

Both scripts accept several files or a wildcard (quote it on Windows too):

```bash
python orthorhombify_v2.py "*.cif" --max-cells 8 --top 1 --outdir output
python make_lammps_inputs.py "output/*_orth_1.cif" --outdir lammps
```

---

## If something goes wrong

| What you see | What to do |
|--------------|------------|
| `python` is not recognized | Install Python and tick **Add to PATH**, then open a **new** terminal |
| `No module named 'ase'` | `pip install -r requirements.txt` |
| `No candidates found` | Raise `--max-cells`, or remove `--anis-max` |
| Very long, thin box | Use the `_orth_2` result, or add `--anis-max 4` |
| `WARNING: shortest interatomic distance ...` | The CIF symmetry did not expand correctly; check the CIF has a `_symmetry_equiv_pos_as_xyz` loop |
| `no charge given for ['S']` | Add the element: `--charge S=-0.1` |
| LAMMPS: `Unknown pair style reaxff` | Rebuild LAMMPS with the `REAXFF` package |
| LAMMPS: `Cannot open ReaxFF potential file` | Put the force-field file in the run folder, or fix `--ffield` |
