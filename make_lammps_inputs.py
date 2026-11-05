#!/usr/bin/env python3
"""
Convert structure files (CIF / POSCAR, e.g. the orthorhombify_v2 outputs) into a
LAMMPS data file with per-element initial charges, plus a ready-to-run LAMMPS
input script for a ReaxFF + QEq equilibration (minimize -> NVT -> NPT). QEq
updates the charges every step; the final structure and charges are saved with
write_data.

Each input  name.cif  produces  <outdir>/name/name.data  and  <outdir>/name/in.name
"""
import argparse, glob, os, sys, warnings
from collections import Counter
import numpy as np

try:
    from ase.io import read
    from ase.io.lammpsdata import write_lammps_data
except Exception:
    print("This script needs ASE. Install with: pip install ase", file=sys.stderr)
    raise

warnings.filterwarnings("ignore", message="crystal system .* is not interpreted")

DEFAULT_CHARGES = {"O": -0.3668, "N": 0.3668, "C": -0.4320, "H": 0.2160}
# LAMMPS atom types are assigned in this order (only elements present are used)
TYPE_ORDER = ["C", "H", "N", "O"]


def parse_charge_overrides(items):
    charges = dict(DEFAULT_CHARGES)
    for item in items or []:
        try:
            el, val = item.split("=")
            charges[el.strip()] = float(val)
        except ValueError:
            sys.exit(f"Bad --charge '{item}', expected e.g. O=-0.3668")
    return charges

def species_order(symbols):
    present = set(symbols)
    return [e for e in TYPE_ORDER if e in present] + sorted(present - set(TYPE_ORDER))

def is_orthogonal(atoms, tol=1e-6):
    cell = atoms.cell.array
    return abs(cell[1, 0]) < tol and abs(cell[2, 0]) < tol and abs(cell[2, 1]) < tol \
        and np.allclose(atoms.cell.angles(), 90.0, atol=1e-6)


INPUT_TEMPLATE = """\
# LAMMPS input for {name}  ({natoms} atoms, {formula})
# ReaxFF + QEq equilibration: minimize -> NVT -> NPT. Charges are re-equilibrated
# by QEq every step, starting from the initial charges stored in {data}.
# Run:  lmp -in in.{name}      (put your ReaxFF force-field file next to this script)

# ---------- settings you may want to change ----------
variable        ffield     string  {ffield}
variable        T          equal   {temp}        # K
variable        P          equal   {press}       # atm
variable        dt         equal   {dt}          # fs
variable        nvt_steps  equal   {nvt_steps}
variable        npt_steps  equal   {npt_steps}
variable        seed       equal   {seed}

# ---------- system ----------
units           real
atom_style      charge
boundary        p p p
read_data       {data}
{replicate}
pair_style      reaxff NULL{lgvdw}
pair_coeff      * * ${{ffield}} {elements}
fix             qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

# ---------- charge monitoring (average charge per element) ----------
compute         q all property/atom q
{groups}
thermo_style    custom step temp press pe etotal vol lx ly lz {qcols}
thermo_modify   flush yes
thermo          100

dump            traj all custom 1000 {name}_traj.lammpstrj id type element x y z q
dump_modify     traj element {elements} sort id

# ---------- 1) energy minimization (fixed box) ----------
{minimize}
# ---------- 2) NVT equilibration ----------
reset_timestep  0
timestep        ${{dt}}
velocity        all create ${{T}} ${{seed}} mom yes rot yes dist gaussian
fix             nvt all nvt temp ${{T}} ${{T}} $(100.0*v_dt)
run             ${{nvt_steps}}
unfix           nvt
write_data      {name}_nvt.data

# ---------- 3) NPT equilibration ({baro}) ----------
fix             npt all npt temp ${{T}} ${{T}} $(100.0*v_dt) {baro} ${{P}} ${{P}} $(1000.0*v_dt)
run             ${{npt_steps}}
unfix           npt

# ---------- save equilibrated structure + updated charges ----------
write_data      {name}_equil.data
{print_q}
"""


def make_input(name, atoms, order, data_name, args):
    groups = "\n".join(
        f"group           g{el} type {i}\n"
        f"compute         q{el} g{el} reduce ave c_q"
        for i, el in enumerate(order, 1))
    qcols = " ".join(f"c_q{el}" for el in order)
    print_q = "print           \"Final average charges: " + \
        "  ".join(f"{el} = $(c_q{el}:%.4f)" for el in order) + "\""
    minimize = ("min_style       cg\n"
                "minimize        1.0e-6 1.0e-8 2000 20000\n") if args.minimize else "# (skipped)\n"
    baro = args.barostat
    if baro == "tri" and is_orthogonal(atoms):
        print("  note: 'tri' needs a triclinic box in LAMMPS; using 'aniso' for this orthogonal cell.")
        baro = "aniso"
    return INPUT_TEMPLATE.format(
        name=name, natoms=len(atoms), formula=atoms.get_chemical_formula(),
        data=data_name, ffield=args.ffield, temp=args.temp, press=args.press,
        dt=args.timestep, nvt_steps=args.nvt_steps, npt_steps=args.npt_steps,
        seed=args.seed, lgvdw=" lgvdw yes" if args.lgvdw else "",
        elements=" ".join(order), groups=groups, qcols=qcols,
        minimize=minimize, baro=baro, print_q=print_q,
        replicate="replicate       {} {} {}\n".format(*args.replicate) if args.replicate else "")


def process(path, charges, args):
    atoms = read(path)
    name = os.path.splitext(os.path.basename(path))[0]
    symbols = atoms.get_chemical_symbols()
    missing = sorted(set(symbols) - set(charges))
    if missing:
        print(f"  {path}: no charge given for {missing}; use --charge El=value", file=sys.stderr)
        return False

    q = np.array([charges[s] for s in symbols])
    net = float(q.sum())
    counts = Counter(symbols)
    print("=" * 78)
    print(f"{path}")
    print(f"  {len(atoms)} atoms  " + "  ".join(f"{el}:{counts[el]}" for el in species_order(symbols)))
    print(f"  net charge with initial charges = {net:+.4f} e")
    if abs(net) > 1e-6:
        if args.neutralize:
            q -= net / len(q)
            shifted = {el: charges[el] - net / len(q) for el in counts}
            print(f"  neutralized: shifted every charge by {-net/len(q):+.6f} e -> "
                  + ", ".join(f"{el} = {v:+.4f}" for el, v in sorted(shifted.items())))
        else:
            print(f"  WARNING: --no-neutralize is set, so the system keeps a net charge of {net:+.4f} e.")
    atoms.set_initial_charges(q)

    order = species_order(symbols)
    outdir = os.path.join(args.outdir, name)
    os.makedirs(outdir, exist_ok=True)
    data_name = f"{name}.data"
    with open(os.path.join(outdir, data_name), "w") as fd:
        write_lammps_data(fd, atoms, specorder=order, atom_style="charge",
                          masses=True, units="real")
    with open(os.path.join(outdir, f"in.{name}"), "w") as fd:
        fd.write(make_input(name, atoms, order, data_name, args))

    box = "orthogonal" if is_orthogonal(atoms) else \
        "triclinic (small tilt, angles {:.3f} {:.3f} {:.3f})".format(*atoms.cell.angles())
    print(f"  atom types: " + "  ".join(f"{i}={el}" for i, el in enumerate(order, 1)) + f"   box: {box}")
    reps = args.replicate or (1, 1, 1)
    lengths = atoms.cell.lengths() * np.array(reps)
    if lengths.min() < args.min_length:
        need = [int(np.ceil(args.min_length / L)) for L in atoms.cell.lengths()]
        print(f"  NOTE: simulated box edges {np.round(lengths, 2).tolist()} A; shortest is below "
              f"{args.min_length} A (QEq/ReaxFF cutoff is 10 A). Consider --replicate {need[0]} {need[1]} {need[2]}")
    print(f"  -> {outdir}{os.sep}{data_name}  and  in.{name}")
    return True


def main():
    ap = argparse.ArgumentParser(description="Make LAMMPS data + input files (with initial charges) "
                                             "from CIF/POSCAR structures.")
    ap.add_argument("files", nargs="+", help="Structure files (CIF, POSCAR, ...); wildcards allowed")
    ap.add_argument("--outdir", default="lammps", help="Output directory. Default: lammps")
    ap.add_argument("--charge", action="append", metavar="El=q",
                    help="Override/add an element charge, e.g. --charge O=-0.35 (repeatable). "
                         "Defaults: " + ", ".join(f"{k}={v}" for k, v in DEFAULT_CHARGES.items()))
    ap.add_argument("--no-neutralize", dest="neutralize", action="store_false",
                    help="Keep the element charges exactly as given, even if the total is not zero. "
                         "By default every charge is shifted by the same amount so the total is exactly zero.")
    ap.add_argument("--ffield", default="ffield.reax", help="ReaxFF force-field file name. Default: ffield.reax")
    ap.add_argument("--lgvdw", action="store_true", help="Force field is ReaxFF-lg (adds 'lgvdw yes').")
    ap.add_argument("--temp", type=float, default=300.0, help="Temperature (K). Default: 300")
    ap.add_argument("--press", type=float, default=1.0, help="Pressure (atm). Default: 1")
    ap.add_argument("--timestep", type=float, default=0.1, help="Timestep (fs). Default: 0.1")
    ap.add_argument("--nvt-steps", type=int, default=10000, help="NVT steps. Default: 10000")
    ap.add_argument("--npt-steps", type=int, default=50000, help="NPT steps. Default: 50000")
    ap.add_argument("--barostat", choices=["aniso", "iso", "tri"], default="aniso",
                    help="NPT coupling. aniso (default) keeps box angles fixed; tri also relaxes tilts.")
    ap.add_argument("--no-minimize", dest="minimize", action="store_false",
                    help="Skip the initial energy minimization.")
    ap.add_argument("--replicate", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                    help="Add 'replicate NX NY NZ' after read_data to enlarge the box in LAMMPS.")
    ap.add_argument("--min-length", type=float, default=20.0,
                    help="Warn if any box edge is shorter than this (A). Default: 20 (2x the 10 A cutoff)")
    ap.add_argument("--seed", type=int, default=4928459, help="Velocity seed. Default: 4928459")
    args = ap.parse_args()

    charges = parse_charge_overrides(args.charge)
    paths = []
    for pattern in args.files:
        paths.extend(sorted(glob.glob(pattern)) or [pattern])
    ok = [process(p, charges, args) for p in paths]
    sys.exit(0 if all(ok) else 2)

if __name__ == "__main__":
    main()
