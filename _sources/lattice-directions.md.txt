# Lattice directions

Every time Orthorhombify makes a new cell it tells you how the **old** axes (a, b, c)
and the **new** axes (a', b', c') are connected, including new axes that are a mix of
two or three old directions. It prints the table and saves it in the `.json` file.

---

## Reading the table

For beta-HMX (`P=[1 0 1; 0 1 0; -4 0 3]`, 7 unit cells):

```text
Lattice directions (old cell a, b, c -> new cell a', b', c'):
  new a' = a + c               along old [1 0 1]         55.05 deg from old a
  new b' = b                   along old [0 1 0]          0.00 deg from old b
  new c' = -4a + 3c            along old [-4 0 3]        42.99 deg from old c
  old a  = (3/7)a' - (1/7)c'   along new [3 0 -1]
  old b  = b'                  along new [0 1 0]
  old c  = (4/7)a' + (1/7)c'   along new [4 0 1]
Planes:
  new (1 0 0) = old (3 0 4)        old (1 0 0) = new (1 0 -4)
  new (0 1 0) = old (0 1 0)        old (0 1 0) = new (0 1 0)
  new (0 0 1) = old (-1 0 1)       old (0 0 1) = new (1 0 3)
```

| Line | Meaning |
|------|---------|
| `new a' = a + c` | The new a axis is the sum of one old a and one old c: a **mix of two old directions** |
| `along old [1 0 1]` | The same thing as a crystal direction of the old cell |
| `55.05 deg from old a` | How far the new axis has turned away from the old one |
| `old a = (3/7)a' - (1/7)c'` | The reverse: the old a axis written in the new cell. Fractions are normal, because the new cell is 7 times bigger |
| `along new [3 0 -1]` | The old a axis as a direction of the new cell (smallest whole numbers) |
| `new (0 0 1) = old (-1 0 1)` | The new cell's c-face is the old (-1 0 1) plane. It contains the new a and b axes |
| `old (0 0 1) = new (1 0 3)` | The old c-face in Miller indices of the new cell |

So for this HMX box, the LAMMPS z axis (new c') runs along the old [-4 0 3] direction,
and the xy plane of the box is the old (-1 0 1) plane. The b axis is untouched.

---

## Any direction or plane

The table only lists the axes. For any other direction or plane use Python:

```python
import orthorhombify as orth

m = orth.LatticeMapping([[1, 0, 1], [0, 1, 0], [-4, 0, 3]])

m.direction_to_new([1, 1, 0])    # old [1 1 0] in the new cell -> [3/7, 1, -1/7]
m.direction_to_old([0, 0, 1])    # new [0 0 1] in the old cell -> [-4, 0, 3]
m.plane_to_new([1, 0, 0])        # old (1 0 0) -> new (1 0 -4)
m.plane_to_old([0, 0, 1])        # new (0 0 1) -> old (-1 0 1)
```

After a search, the same object is on every result: `best.mapping`.

### The rules behind it

With `new cell = P @ old cell` (rows are the cell vectors):

| | old → new | new → old |
|---|-----------|-----------|
| Direction [u v w] | `u' = u @ P⁻¹` | `u = u' @ P` |
| Plane (h k l) | `h' = P @ h` | `h = P⁻¹ @ h'` |

A direction lies in a plane when `h · u = 0`, and that is the same in both cells.

---

## The `.json` file

Next to every output (`<name>_orth_1.json`, `<name>_original.json`, ...):

Shortened (`...` marks left-out entries):

```text
{
  "orthorhombify_version": "2.1.0",
  "source": "beta_HMX.cif",
  "n_cells": 7,
  "max_angle_deviation_deg": 0.001025,
  "old_cell": {"a": 6.5209, "b": 10.761, "c": 7.3062, "alpha": 90.0, "beta": 102.058, "gamma": 90.0, "n_atoms": 56, ...},
  "new_cell": {"a": 8.717642, "b": 10.761, "c": 37.411818, "alpha": 90.0, "beta": 89.998975, "gamma": 90.0, "n_atoms": 392, ...},
  "P": [["1", "0", "1"], ["0", "1", "0"], ["-4", "0", "3"]],
  "P_inverse": [["3/7", "0", "-1/7"], ["0", "1", "0"], ["4/7", "0", "1/7"]],
  "det_P": "7",
  "new_axes": {"a'": {"expression": "a + c", "old_direction": [1, 0, 1]}, ...},
  "old_axes": {"a": {"expression": "(3/7)a' - (1/7)c'", "new_direction": [3, 0, -1]}, ...},
  "new_planes_in_old": {"(1 0 0)": [3, 0, 4], "(0 1 0)": [0, 1, 0], "(0 0 1)": [-1, 0, 1]},
  "old_planes_in_new": {"(1 0 0)": [1, 0, -4], "(0 1 0)": [0, 1, 0], "(0 0 1)": [1, 0, 3]},
  "axis_angles_deg": {"a'": 55.045, "b'": 0.0, "c'": 42.986}
}
```

`orthorhombify-transform --inverse` reads this file to take a structure back to the
original cell. See [Going back and other cells](transform.md).
