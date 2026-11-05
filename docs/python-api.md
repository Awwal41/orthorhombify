# Using it from Python

Everything the commands do is also available as functions, so you can use Orthorhombify
in your own scripts and notebooks.

```python
import orthorhombify as orth
```

---

## The whole workflow in a few lines

```python
from ase.io import read, write
import orthorhombify as orth

atoms = read("my_crystal.cif")

# 1. find the best near-orthorhombic supercells
best = orth.find_supercells(atoms, max_cells=8, top=1)[0]
print(best.n_cells, best.lengths, best.angles, best.max_dev)

# 2. build it and save it
supercell = best.build(atoms)
write("my_crystal_orth.cif", supercell)

# 3. LAMMPS files with neutralized starting charges
files = orth.write_lammps_inputs(supercell, "lammps/my_crystal", "my_crystal",
                                 ffield="ffield.reax", replicate=(3, 2, 1))
print(files.data_file, files.input_file)
print(files.charges.per_element)      # final charge of each element
```

Try it on the bundled example by using `read(orth.example_path())`.

---

## Finding supercells

```python
orth.detect_system(atoms)              # ('monoclinic', 1)  -> unique axis b
orth.best_angle_by_size(atoms, 8)      # {1: 12.06, 2: 6.65, ..., 7: 0.001, 8: 0.001}

results = orth.find_supercells(atoms, max_cells=8, top=3)
small   = orth.find_supercells(atoms, max_cells=8, top=1, tol=1.0)       # smallest box within 1 degree
compact = orth.find_supercells(atoms, max_cells=8, top=1, anis_max=4.0)  # no long thin boxes
```

`find_supercells` also accepts a 3x3 cell matrix or `[a, b, c, alpha, beta, gamma]`
instead of `Atoms`, which is handy for quick checks:

```python
orth.find_supercells([9.936, 12.61, 12.277, 90, 90, 116.66], max_cells=8, top=1)[0].angles
# (90.0, 90.0, 89.90)
```

Each result is a `Candidate`:

| Attribute | Meaning |
|-----------|---------|
| `P` | 3x3 integer matrix: new cell rows = `P @` old cell rows |
| `n_cells` | Number of unit cells in the new box |
| `lengths` | a, b, c (Å) |
| `angles` | alpha, beta, gamma (degrees) |
| `max_dev` | Largest angle deviation from 90 degrees |
| `anisotropy` | Longest edge / shortest edge |
| `build(atoms, force_orthogonal=False)` | Returns the supercell `Atoms` |

With `force_orthogonal=True` the box is made exactly 90/90/90 and the strain is stored in
`supercell.info["orthorhombify_strain"]`.

---

## Charges and LAMMPS

```python
orth.DEFAULT_CHARGES                   # {'O': -0.3668, 'N': 0.3668, 'C': -0.432, 'H': 0.216}

info = orth.assign_charges(supercell, charges={"O": -0.35}, neutralize=True)
info.net_before, info.shift, info.per_element, info.charges

orth.suggest_replicate(supercell)      # (3, 2, 1), or None if the box is big enough

text = orth.make_input_script(supercell, "my_crystal", "my_crystal.data", temp=500)
```

`write_lammps_inputs` takes the same options as the `orthorhombify-lammps` command
(`charges`, `neutralize`, `ffield`, `lgvdw`, `temp`, `press`, `timestep`, `nvt_steps`,
`npt_steps`, `barostat`, `minimize`, `replicate`, `seed`). See the [API reference](api.rst)
for the full signatures.

---

## Lattice directions

```python
m = best.mapping                       # LatticeMapping, exact fractions
print(m.summary(atoms.cell))           # the table printed by the command
m.new_axis_expressions()               # ["a + c", "b", "-4a + 3c"]
m.direction_to_new([1, 1, 0])          # old [1 1 0] in the new cell: [3/7, 1, -1/7]
m.plane_to_old([0, 0, 1])              # new (0 0 1) plane in the old cell: [-1, 0, 1]
m.inverse()                            # the mapping the other way round
```

Any matrix works, not only results of the search: `orth.LatticeMapping([[1, 0, 1], [0, 1, 0], [-4, 0, 3]])`.
More in [Lattice directions](lattice-directions.md).

---

## Other cells and going back

```python
# any cell change: new cell = P @ old cell (fractions allowed)
hexagonal = orth.transform(orthohexagonal, "1 0 0; -1/2 1/2 0; 0 0 1")

# back to the original cell after LAMMPS (replicate is worked out from the box)
eq = read("lammps/my_crystal/my_crystal_equil.data", format="lammps-data",
          atom_style="charge", units="real")
original = orth.back_to_original(eq, "output/my_crystal_orth_1.json")
original.info["orthorhombify_replicate"], original.info["orthorhombify_spread"]
```

More in [Going back and other cells](transform.md).

---

## Many crystals in a loop

```python
import glob, os
from ase.io import read
import orthorhombify as orth

for path in glob.glob("cifs/*.cif"):
    name = os.path.splitext(os.path.basename(path))[0]
    atoms = read(path)
    best = orth.find_supercells(atoms, max_cells=8, top=1, anis_max=5)[0]
    sc = best.build(atoms)
    orth.write_lammps_inputs(sc, f"lammps/{name}", name, replicate=orth.suggest_replicate(sc))
    print(f"{name}: N={best.n_cells}  max dev={best.max_dev:.3f} deg  atoms={len(sc)}")
```
