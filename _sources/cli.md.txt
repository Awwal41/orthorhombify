# CLI reference

Installing the package adds three commands. All print this list with `--help` and
their version with `--version`.

| Command | What it does | From a clone, without installing |
|---------|--------------|----------------------------------|
| `orthorhombify` | Find near-orthorhombic supercells (also `python -m orthorhombify`) | `python orthorhombify_v2.py` |
| `orthorhombify-lammps` | LAMMPS data + input files with charges | `python make_lammps_inputs.py` |
| `orthorhombify-transform` | Any other cell, or back to the original | - |

---

## orthorhombify

```bash
orthorhombify FILE [FILE ...] [options]
```

`FILE` is one or more CIF files. Wildcards like `"*.cif"` work on every system.
For every result it writes `.cif`, `.POSCAR` and `.json` and prints the
[lattice directions](lattice-directions.md).

### Main options

| Option | Default | What it does |
|--------|---------|--------------|
| `--max-cells N` | 6 | Largest box to try, in unit cells |
| `--top K` | 2 | How many results to write (1 to 5) |
| `--tol T` | off | Pick the **smallest** box within T degrees of 90 |
| `--anis-max R` | off | Skip boxes whose longest edge is more than R times the shortest |
| `--force-orthogonal` | off | Make the final box exactly 90/90/90 (small strain, printed) |
| `--outdir DIR` | `.` | Where to write the files |

### Less common options

| Option | Default | What it does |
|--------|---------|--------------|
| `--mode` | `auto` | `auto` (from the cell shape), `monoclinic` (2D search) or `triclinic` (3D search) |
| `--method` | `vectors` | 3D search: `vectors` (fast, thorough) or `hnf` (the original v1 search) |
| `--angle-tol` | 0.001 | How close to 90 an angle must be to count as 90 |
| `--max-index M` | 12 (2D), 5 (3D) | Largest whole-number coefficient tried for the new axes |

### `--method hnf` only

| Option | Default | What it does |
|--------|---------|--------------|
| `--unimod M` | 2 | How strongly axes may be mixed (1 or 2) |
| `--fast` | off | Only diagonal supercells. Faster, worse angles |
| `--gate X` | off | Skip candidates whose angle errors add up to more than X degrees |
| `--stop-on-tol` | off | Stop as soon as the angle errors add up to at most 3 x `--tol` |

---

## orthorhombify-lammps

```bash
orthorhombify-lammps FILE [FILE ...] [options]
```

`FILE` is any structure ASE can read: CIF, POSCAR, XYZ with a cell, and so on.

### Charges

| Option | Default | What it does |
|--------|---------|--------------|
| `--charge El=q` | O=-0.3668, N=0.3668, C=-0.4320, H=0.2160 | Set an element charge. Repeat for more |
| `--no-neutralize` | off | Keep the charges as given even if the box is not neutral |

### Force field and run

| Option | Default | What it does |
|--------|---------|--------------|
| `--ffield FILE` | `ffield.reax` | ReaxFF force-field file name |
| `--lgvdw` | off | ReaxFF-lg force field |
| `--temp K` | 300 | Temperature (K) |
| `--press P` | 1 | Pressure (atm) |
| `--timestep T` | 0.1 | Timestep (fs) |
| `--nvt-steps N` | 10000 | NVT steps |
| `--npt-steps N` | 50000 | NPT steps |
| `--barostat` | `aniso` | `aniso` (keep angles), `iso`, or `tri` (relax angles too) |
| `--no-minimize` | off | Skip the energy minimization |
| `--seed S` | 4928459 | Random seed for starting velocities |

### Box and output

| Option | Default | What it does |
|--------|---------|--------------|
| `--replicate NX NY NZ` | off | Enlarge the box in LAMMPS |
| `--min-length L` | 20 | Warn if a box edge is shorter than L Å |
| `--outdir DIR` | `lammps` | Where to write the folders |

---

## orthorhombify-transform

```bash
orthorhombify-transform FILE [FILE ...] (--matrix "P" | --inverse JSON) [options]
```

`FILE` is any structure ASE can read; LAMMPS `.data` files (atom_style charge) are read
directly. See [Going back and other cells](transform.md).

| Option | Default | What it does |
|--------|---------|--------------|
| `--matrix "P"` | | New cell = P @ old cell. Rows split by `;`, fractions allowed, 3 numbers = diagonal |
| `--inverse JSON` | | Back to the original cell, using the `.json` written by `orthorhombify` |
| `--replicate NX NY NZ` | worked out | With `--inverse`: the LAMMPS replicate to undo |
| `--reduce` | `auto` | Smaller cell: `average` the copies of each atom, or `dedup` (merge folded atoms) |
| `--tol T` | 0.5 | Merge distance / allowed spread between copies (Å) |
| `--outdir DIR` | `.` | Where to write the files |
| `--suffix S` | `_original` / `_transformed` | Added to the output file name |
