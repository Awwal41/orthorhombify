# Orthorhombify (the supercell search)

`orthorhombify_v2.py` finds a new box made of whole lattice vectors whose angles are as
close to 90 degrees as possible. Because it only recombines lattice vectors, the crystal
is **not strained**: every atom keeps its exact environment. The new box is just bigger.

```bash
python orthorhombify_v2.py my_crystal.cif --max-cells 8 --top 2 --outdir output
```

---

## What it does with each kind of cell

The script reads the three cell angles and picks a method by itself.

| Your cell | What happens |
|-----------|--------------|
| **Monoclinic** (one angle is not 90) | Exact 2D search. See below. |
| **Triclinic** (two or three angles are not 90) | 3D search over integer supercells and axis shears (the original v1 method). |
| **Orthorhombic** (all 90) | Nothing to do. The cell is written out unchanged. |

### Monoclinic

A monoclinic cell has one **unique axis** that is already at 90 degrees to the other two.
Only the angle between the other two vectors (call them **p** and **q**) needs fixing.

The script tries every pair of new in-plane vectors

```text
v1 = i1*p + j1*q
v2 = i2*p + j2*q        (i, j are whole numbers)
```

keeps the unique axis as it is, and scores how close `v1` and `v2` are to 90 degrees.
The new box holds `N = |i1*j2 - j1*i2|` unit cells.

It works for every monoclinic setting:

| Space-group symbol | Unique axis | Angle that is fixed |
|--------------------|-------------|---------------------|
| `P 21/c`, `P 21/n`, `C 2/c` | b | beta |
| `P 1 1 21/b`, `P 1 1 21/a` | c | gamma |
| `P 21/b 1 1` | a | alpha |

### Triclinic

All three angles can be off, so the search is in 3D: every Hermite normal form supercell
up to `--max-cells`, combined with a small set of axis permutations and shears. It is
slower and less exhaustive than the monoclinic search. The triclinic-only flags
(`--unimod`, `--fast`, `--gate`, `--stop-on-tol`) control it. See [CLI](cli.md).

---

## Choosing a result

For monoclinic cells the script prints the best angle you can get for each box size:

```text
Best possible angle deviation by supercell size N (x unit cell):
  N=1: 20.182  N=2: 0.470  N=3: 0.470  N=4: 0.130  N=5: 0.130 ...
```

Bigger boxes give better angles but more atoms. Three ways to choose:

| You want | Use |
|----------|-----|
| The best angle up to some size | `--max-cells 8` (default ranking) |
| The smallest box that is "good enough" | `--tol 0.5` (smallest N within 0.5 degrees) |
| No long, thin boxes | `--anis-max 4` (longest edge at most 4x the shortest) |

`--top 2` writes the two best results so you can compare them.

---

## Exactly 90 degrees

The result is usually a few hundredths of a degree off 90. Most MD codes are fine with
that: LAMMPS handles the small tilt. If your code needs an exactly orthogonal box:

```bash
python orthorhombify_v2.py my_crystal.cif --force-orthogonal
```

This squares the box and keeps the fractional coordinates. It introduces a tiny strain,
which is printed (for example `max Green-Lagrange strain = 5.67e-04` for FOX-7).

---

## Output files

| File | Contents |
|------|----------|
| `<name>_orth_1.cif` | Best supercell, space group P1, all atoms listed |
| `<name>_orth_1.POSCAR` | Same structure for VASP / other codes |
| `<name>_orth_2.*` | Runner-up (with `--top 2`) |

The cell is rotated so **a** points along x and **b** lies in the xy plane (the LAMMPS
convention). Atoms are wrapped into the box, so molecules may look cut at the edges.
That is normal for a periodic cell.

---

## Checks the script does for you

- Warns if two atoms are closer than 0.5 Å, which usually means the CIF symmetry did
  not expand correctly.
- Checks that the supercell has exactly `N x` the atoms of the unit cell.
