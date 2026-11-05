# Examples

`beta_HMX.cif` is beta-HMX from the Crystallography Open Database
([COD 4110146](https://www.crystallography.net/cod/4110146.html)), monoclinic `P 21/n`.

```bash
python orthorhombify_v2.py examples/beta_HMX.cif --max-cells 8 --top 1 --outdir output
python make_lammps_inputs.py output/beta_HMX_orth_1.cif --outdir lammps
```

Expected: a 7-cell box (392 atoms), 8.72 x 10.76 x 37.41 Å, beta = 89.999 degrees.

Results for MDNT, TNT and FOX-7 (from the CSD, not included) are in the
[examples page](https://awwal41.github.io/orthorhombify/examples.html) of the docs.
