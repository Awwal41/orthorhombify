# Going back and other cells

`orthorhombify-transform` changes a structure from one cell to another cell of the same
lattice, in **either direction**:

- **Back to the original cell** after an MD run, to compare with experiment.
- **From orthorhombic to another setting**: hexagonal, monoclinic, a smaller cell, or any
  matrix you give it.

```bash
orthorhombify-transform FILE --inverse <name>_orth_1.json     # back to the original cell
orthorhombify-transform FILE --matrix "1 0 0; -1/2 1/2 0; 0 0 1"   # any transformation
```

It writes `<name>_original` (or `<name>_transformed`) as `.cif`, `.POSCAR` and `.json`,
and prints the [lattice directions](lattice-directions.md) between the two cells.

---

## After LAMMPS: back to the original cell

The usual workflow:

```bash
orthorhombify beta_HMX.cif --max-cells 8 --top 1 --outdir output
orthorhombify-lammps output/beta_HMX_orth_1.cif --outdir lammps --replicate 3 2 1
# ... run LAMMPS in lammps/beta_HMX_orth_1/ ...
orthorhombify-transform lammps/beta_HMX_orth_1/beta_HMX_orth_1_equil.data \
    --inverse output/beta_HMX_orth_1.json
```

The last command:

1. Reads the LAMMPS data file (elements come from the masses).
2. Works out the `replicate 3 2 1` from the atom count and box size (or pass `--replicate`).
3. Undoes the replicate and the supercell.
4. **Averages** all copies of each atom (3 x 2 x 1 x 7 = 42 copies here) into one unit cell.
5. Prints the equilibrated cell next to the original.

For a test box stretched by +1.2 %, -0.5 % and +0.8 % along x, y, z, with 0.08 Å of
random noise on every atom (a stand-in for an NPT run):

```text
input : a,b,c = 26.4668, 21.4144, 37.7111 A   alpha,beta,gamma = 90.000, 90.000, 90.000 deg
output: a,b,c = 6.5817, 10.7072, 7.3782 A   alpha,beta,gamma = 90.000, 101.838, 90.000 deg
56 atoms in the new cell (det = 1/42); method: average, copies of each atom agree within 0.420 A
undid LAMMPS replicate 3 2 1
change vs the original cell in beta_HMX_orth_1.json:
  a 6.5209 -> 6.5817 (+0.93%)   b 10.7610 -> 10.7072 (-0.50%)   c 7.3062 -> 7.3782 (+0.99%)
  alpha 90.0000 -> 90.0000 (+0.00%)   beta 102.0580 -> 101.8380 (-0.22%)   gamma 90.0000 -> 90.0000 (+0.00%)
```

Stretching the box along x, y, z changes the monoclinic a, c **and beta** of the original
cell, because the old axes are mixtures of the new ones.

So you get the **equilibrated monoclinic lattice parameters** (a, b, c, beta) to compare
with the experimental CIF, and an averaged structure in the original setting.

This works because Orthorhombify and LAMMPS store the copies of each atom one after
another. The script checks this: copies of an atom must agree within `--tol` (0.5 Å).
If they do not, it falls back to merging (below).

---

## From orthorhombic to other cells

Give the transformation with `--matrix`, as `new cell = P @ old cell`. Rows are
separated by `;`, and fractions are allowed (they make the cell smaller).

| Goal | `--matrix` |
|------|------------|
| Orthohexagonal (a, a+2b, c) → hexagonal | `"1 0 0; -1/2 1/2 0; 0 0 1"` |
| Hexagonal → orthohexagonal | `"1 0 0; 1 2 0; 0 0 1"` |
| C-centred orthorhombic → primitive (monoclinic shape) | `"1/2 -1/2 0; 1/2 1/2 0; 0 0 1"` |
| Swap axes: a b c → b c a | `"0 1 0; 0 0 1; 1 0 0"` |
| Plain 2 x 2 x 1 supercell | `"2 2 1"` |
| Exactly undo any `orthorhombify` result | `--inverse <name>_orth_1.json` |

Example, an orthorhombic magnesium cell (4 atoms) to hexagonal:

```text
orthorhombify-transform Mg_ortho.cif --matrix "1 0 0; -1/2 1/2 0; 0 0 1"
  input : a,b,c = 3.2100, 5.5599, 5.2100 A   alpha,beta,gamma = 90.000, 90.000, 90.000 deg
  output: a,b,c = 3.2100, 3.2100, 5.2100 A   alpha,beta,gamma = 90.000, 90.000, 120.000 deg
  2 atoms in the new cell (det = 1/2); method: average, copies of each atom agree within 0.000 A
  Lattice directions (old cell a, b, c -> new cell a', b', c'):
    new a' = a                  along old [1 0 0]          0.00 deg from old a
    new b' = -(1/2)a + (1/2)b   along old [-1 1 0]        30.00 deg from old b
    new c' = c                  along old [0 0 1]          0.00 deg from old c
    ...
```

### How atoms are placed in a smaller cell

| `--reduce` | What it does | Use it for |
|------------|--------------|------------|
| `auto` (default) | `average` if the atom order allows it, otherwise `dedup` | Almost always |
| `average` | Averages the copies of each atom | Orthorhombify / LAMMPS output, MD snapshots |
| `dedup` | Folds every atom into the new cell and merges atoms of the same element closer than `--tol` | Structures from elsewhere, with any atom order |

The new cell must hold a whole number of atoms (`atoms x det(P)`). If the matrix does not
describe a cell of this lattice, you get an error, or a warning that the atom count is off.

---

## From Python

```python
from ase.io import read
import orthorhombify as orth

md = read("beta_HMX_orth_1_equil.data", format="lammps-data", atom_style="charge", units="real")
unit = orth.back_to_original(md, "output/beta_HMX_orth_1.json")     # replicate worked out
print(unit.cell.cellpar())

hexa = orth.transform(read("Mg_ortho.cif"), "1 0 0; -1/2 1/2 0; 0 0 1")
```
