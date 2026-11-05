# Examples

Four energetic molecular crystals, all monoclinic, all space group 14.

**Try it now:** beta-HMX comes from the open
[Crystallography Open Database](https://www.crystallography.net/cod/4110146.html)
(COD 4110146) and ships inside the package:

```python
from ase.io import read
import orthorhombify as orth

atoms = read(orth.example_path("beta_HMX"))
best = orth.find_supercells(atoms, max_cells=8, top=1)[0]
print(best.n_cells, best.angles)       # 7 (90.0, 89.999, 90.0)
```

or from the terminal (in a clone it is `src/orthorhombify/data/beta_HMX.cif`):

```bash
orthorhombify path/to/beta_HMX.cif --max-cells 8 --top 1 --outdir output
orthorhombify-lammps output/beta_HMX_orth_1.cif --outdir lammps
```

The other three come from the Cambridge Structural Database (CSD) and are **not** in this
repository. Download them from the [CCDC](https://www.ccdc.cam.ac.uk/structures/) with
the refcodes below, then run the same commands on them.

---

## Summary

| Crystal | Source | Setting | Input angle | Best result | N | Atoms |
|---------|--------|---------|-------------|-------------|---|-------|
| beta-HMX | COD 4110146 | P 21/n | beta = 102.06 | beta = 89.999 | 7 | 392 |
| MDNT | CSD SEGFAL | P 1 1 21/b | gamma = 116.66 | gamma = 89.90 | 7 | 840 |
| TNT (monoclinic) | CSD ZZZMUC06 | P 1 1 21/b | gamma = 111.01 | gamma = 90.13 | 4 | 672 |
| FOX-7 | CSD SEDTUQ03 | P 21/n | beta = 90.065 | already the best cell | 1 | 56 |

N is the number of unit cells in the new box.

---

## beta-HMX

```text
N=1: 12.058  N=2: 6.647  N=3: 5.177  N=4: 0.833  N=5: 0.833  N=6: 0.833  N=7: 0.001
[1] N=7 (392 atoms)  P=[1 0 1; 0 1 0; -4 0 3]
    a,b,c = 8.7176, 10.7610, 37.4118 A   angles = 90.000, 89.999, 90.000 deg
```

At N = 7 the box is orthorhombic to 0.001 degrees. For comparison, the original v1 script
(`--max-cells 12`) gave 90.83 degrees with 448 atoms.

## MDNT

```text
N=1: 20.787  N=2: 3.504  N=3: 1.644  N=5: 1.543  N=7: 0.100
[1] N=7  a,b,c = 9.9360, 78.8857, 12.2770 A   gamma = 89.900   (long and thin)
[2] N=7  a,b,c = 18.1393, 43.2105, 12.2770 A  gamma = 90.110   (more compact)
```

Both results hold 840 atoms. Result 2 has a better shape for MD.

## TNT (monoclinic)

```text
N=1: 20.182  N=2: 0.470  N=4: 0.130
[1] N=4 (672 atoms)  a,b,c = 21.4070, 56.0840, 6.0932 A   gamma = 90.130
```

For a smaller box, `--tol 0.5` picks N = 2: 15.02 x 39.97 x 6.09 Å, 336 atoms,
gamma = 89.53.

The starting charges are far from neutral for TNT (-97.42 e for 672 atoms), so the
neutralizing shift is large (+0.145 e per atom). Consider giving TNT its own charges.
See [LAMMPS files and charges](lammps.md#neutralizing-on-by-default).

## FOX-7

```text
N=1 ... N=8: 0.065
[1] N=1 (56 atoms)  beta = 90.065
```

The unit cell is already the closest to 90 degrees. No supercell does better. If you need
exactly 90, `--force-orthogonal` squares it with a strain of only 5.7e-4. The box is
small (6.9 x 6.6 x 11.3 Å), so use `--replicate 3 4 2` for LAMMPS.
