# Orthorhombify (the supercell search)

The `orthorhombify` command finds a new box made of whole lattice vectors whose angles are as
close to 90 degrees as possible. Because it only recombines lattice vectors, the crystal
is **not strained**: every atom keeps its exact environment. The new box is just bigger.

```bash
orthorhombify my_crystal.cif --max-cells 8 --top 2 --outdir output
```

---

## What it does with each kind of cell

The script reads the cell lengths and angles and picks a method by itself. All seven
crystal systems are covered:

| Cell shape | Looks like | What happens | Typical result |
|------------|------------|--------------|----------------|
| **Cubic, tetragonal, orthorhombic** | all angles 90 | Nothing to do; written unchanged | N = 1 |
| **Monoclinic** | one angle is not 90 | Exact 2D search (below) | within 0.001 to 0.5 degrees, N = 1 to 8 |
| **Hexagonal** (and trigonal in hexagonal axes) | a = b, gamma = 120 (or 60) | Exact 2D search | **exactly 90**, N = 2 (a, a + 2b, c) |
| **Rhombohedral** (trigonal in rhombohedral axes) | a = b = c, alpha = beta = gamma | 3D search | **exactly 90**, N = 6 |
| **Triclinic** | anything else | 3D search (below) | usually 1 to 4 degrees with N up to 8 to 12 |

The name printed (for example `Cell shape is HEXAGONAL`) describes the **shape of the
cell**, not its symmetry. A primitive cell of a centred lattice can look monoclinic or
rhombohedral even if the crystal is orthorhombic or cubic; the search still works.

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

### Hexagonal

A hexagonal cell is a monoclinic-shaped cell with a = b and gamma = 120. The 2D search
finds the classic orthohexagonal cell `a' = a`, `b' = a + 2b`, `c' = c`, which is exactly
rectangular, with twice the atoms.

### Triclinic and rhombohedral

All three angles can be off, so the search is in 3D. Every short lattice vector
`u*a + v*b + w*c` (whole numbers up to `--max-index`, default 5) is tried; the most
perpendicular pairs are completed with the best third vector, keeping the box at most
`--max-cells` unit cells. The three new axes are then ordered so each one is as close
as possible to the old axis of the same name.

For a rhombohedral cell this finds the exact orthorhombic box of 6 cells. For a general
triclinic cell no exact answer exists; bigger boxes get closer:

| Triclinic test cell | N <= 3 | N <= 8 | N <= 12 |
|---------------------|--------|--------|---------|
| 6.5, 11.0, 8.7 Å, 95/103/88 | 4.04 | 3.78 | |
| 7.1, 8.3, 9.2 Å, 101/96/105 | | 2.23 | 1.42 |

(largest angle deviation from 90, degrees)

`--method hnf` runs the original v1 search instead (with `--unimod`, `--fast`, `--gate`,
`--stop-on-tol`). It is slower and usually worse: 6.0 degrees for the second cell above.

---

## Choosing a result

For monoclinic, hexagonal, rhombohedral and triclinic cells the script prints the best
angle it found for each box size:

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
orthorhombify my_crystal.cif --force-orthogonal
```

This squares the box and keeps the fractional coordinates. It introduces a tiny strain,
which is printed (for example `max Green-Lagrange strain = 5.67e-04` for FOX-7).

---

## Output files

| File | Contents |
|------|----------|
| `<name>_orth_1.cif` | Best supercell, space group P1, all atoms listed |
| `<name>_orth_1.POSCAR` | Same structure for VASP / other codes |
| `<name>_orth_1.json` | The transformation, both cells, and how every axis and plane maps (see [Lattice directions](lattice-directions.md)). Keep it: it lets you go back later |
| `<name>_orth_2.*` | Runner-up (with `--top 2`) |

The cell is rotated so **a** points along x and **b** lies in the xy plane (the LAMMPS
convention). Atoms are wrapped into the box, so molecules may look cut at the edges.
That is normal for a periodic cell.

---

## Checks the script does for you

- Warns if two atoms are closer than 0.5 Å, which usually means the CIF symmetry did
  not expand correctly.
- Checks that the supercell has exactly `N x` the atoms of the unit cell.
