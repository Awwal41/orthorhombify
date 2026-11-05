# CLI reference

Both scripts print the same list with `--help`.

---

## orthorhombify_v2.py

```bash
python orthorhombify_v2.py FILE [FILE ...] [options]
```

`FILE` is one or more CIF files. Wildcards like `"*.cif"` work on every system.

### Main options

| Option | Default | What it does |
|--------|---------|--------------|
| `--max-cells N` | 6 | Largest box to try, in unit cells |
| `--top K` | 2 | How many results to write (1 to 5) |
| `--tol T` | off | Monoclinic: pick the **smallest** box within T degrees of 90 |
| `--anis-max R` | off | Skip boxes whose longest edge is more than R times the shortest |
| `--force-orthogonal` | off | Make the final box exactly 90/90/90 (small strain, printed) |
| `--outdir DIR` | `.` | Where to write the files |

### Less common options

| Option | Default | What it does |
|--------|---------|--------------|
| `--mode` | `auto` | `auto`, `monoclinic` or `triclinic`. Overrides detection |
| `--angle-tol` | 0.001 | How close to 90 an angle must be to count as 90 |
| `--max-index M` | 12 | Monoclinic: largest whole-number coefficient tried |

### Triclinic only

| Option | Default | What it does |
|--------|---------|--------------|
| `--unimod M` | 2 | How strongly axes may be mixed (1 or 2) |
| `--fast` | off | Only diagonal supercells. Faster, worse angles |
| `--gate X` | off | Skip candidates whose angle errors add up to more than X degrees |
| `--stop-on-tol` | off | Stop as soon as the angle errors add up to at most 3 x `--tol` |

---

## make_lammps_inputs.py

```bash
python make_lammps_inputs.py FILE [FILE ...] [options]
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
