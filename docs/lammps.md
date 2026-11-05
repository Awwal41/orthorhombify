# LAMMPS files and charges

`orthorhombify-lammps` turns a structure (CIF or POSCAR) into a LAMMPS data file with
starting charges, plus an input script that equilibrates it with ReaxFF + QEq.

```bash
orthorhombify-lammps output/my_crystal_orth_1.cif --outdir lammps --ffield ffield.reax
```

For each structure you get a folder:

```text
lammps/my_crystal_orth_1/
    my_crystal_orth_1.data     atom_style charge, real units
    in.my_crystal_orth_1       LAMMPS input script
```

---

## Starting charges

| Element | Charge (e) | LAMMPS type |
|---------|-----------|-------------|
| C | -0.4320 | 1 |
| H | +0.2160 | 2 |
| N | +0.3668 | 3 |
| O | -0.3668 | 4 |

Only elements in your crystal get a type, always in the order C, H, N, O.

Change or add a value with `--charge` (repeat it for more elements):

```bash
orthorhombify-lammps crystal.cif --charge O=-0.35 --charge S=-0.10
```

These are only **starting** values. QEq recomputes the charges at every step of the run.

---

## Neutralizing (on by default)

The charges above add up to zero for some formulas and not for others:

| Crystal | Formula per molecule | Net charge of the box | After neutralizing |
|---------|---------------------|-----------------------|--------------------|
| FOX-7 | C2H4N4O4 | 0.000 e | unchanged |
| beta-HMX | C4H8N8O8 | 0.000 e | unchanged |
| MDNT (840 atoms) | C3H3N5O4 | -15.75 e | every charge +0.0187 |
| TNT (672 atoms) | C7H5N3O6 | -97.42 e | every charge +0.1450 |

A charged periodic box is a problem for the electrostatics. So by default the script
shifts **every** charge by the same amount until the total is exactly zero, and prints
the new values:

```text
neutralized: shifted every charge by +0.144971 e -> C = -0.2870, H = +0.3610, N = +0.5118, O = -0.2218
```

A large shift (like TNT) is a hint that these charges were made for a different
molecule. You can give better ones with `--charge`, or keep the raw values with
`--no-neutralize`.

---

## What the input script does

| Stage | What runs | Default |
|-------|-----------|---------|
| Setup | `pair_style reaxff`, `fix qeq/reaxff` (charges updated every step) | `ffield.reax` |
| 1. Minimize | Conjugate-gradient energy minimization, fixed box | on (`--no-minimize` skips it) |
| 2. NVT | Constant volume, Nosé-Hoover thermostat | 10,000 steps, 300 K |
| 3. NPT | Constant pressure, each axis on its own (`aniso`) | 50,000 steps, 1 atm |
| Save | `write_data` after NVT and after NPT | `*_nvt.data`, `*_equil.data` |

Timestep is 0.1 fs. With `aniso` the box angles stay fixed, so the box stays
orthorhombic. Use `--barostat tri` to also let the angles relax.

All the settings sit at the top of the input file, so you can edit them there:

```text
variable        ffield     string  ffield.reax
variable        T          equal   300.0        # K
variable        P          equal   1.0          # atm
variable        dt         equal   0.1          # fs
variable        nvt_steps  equal   10000
variable        npt_steps  equal   50000
```

### Watching the charges

The thermo output has one column per element with its average charge:

```text
Step Temp Press PotEng TotEng Volume Lx Ly Lz c_qC c_qH c_qN c_qO
```

The trajectory `*_traj.lammpstrj` stores the charge of every atom every 1000 steps
(open it in Ovito and colour by `q`). The run ends by printing the final averages.

---

## Box size

ReaxFF and QEq use a 10 Å cutoff. If any box edge is shorter than 20 Å the script
prints a suggestion:

```text
NOTE: simulated box edges [6.93, 6.62, 11.31] A; shortest is below 20.0 A ... Consider --replicate 3 4 2
```

Add it to the command. This puts `replicate 3 4 2` in the input script, so LAMMPS
builds the bigger box at the start of the run:

```bash
orthorhombify-lammps output/FOX7_orth_1.cif --replicate 3 4 2
```

---

## Using ReaxFF-lg

For ReaxFF-lg force fields (common for energetic materials) add `--lgvdw`. This writes
`pair_style reaxff NULL lgvdw yes`.

## Running

```bash
cd lammps/my_crystal_orth_1
lmp -in in.my_crystal_orth_1
```

LAMMPS must be built with the `REAXFF` package. Put the force-field file in the same
folder.
