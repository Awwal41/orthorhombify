"""Command-line tools: ``orthorhombify`` and ``orthorhombify-lammps``."""
from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings

import numpy as np

from . import __version__
from .core import (ANGLE_NAMES, AXES, ORTHOGONAL_SYSTEMS, best_angle_by_size, detect_system,
                   find_supercells, min_distance)
from .lammps import DEFAULT_CHARGES, is_orthogonal, suggest_replicate, write_lammps_inputs
from .mapping import LatticeMapping, write_mapping_json
from .transform import back_to_original, parse_matrix, transform


def _read(path):
    from ase.io import read
    with warnings.catch_warnings():
        # ASE warns about the 'monoclinic' cell-setting tag even when the CIF's explicit
        # symmetry operators (which it uses) are present; min_distance() catches real problems.
        warnings.filterwarnings("ignore", message="crystal system .* is not interpreted")
        if path.lower().endswith((".data", ".lmp")):
            return read(path, format="lammps-data", atom_style="charge", units="real")
        return read(path)

def _write_structure(out, atoms):
    from ase.io import write
    write(out + ".cif", atoms)
    write(out + ".POSCAR", atoms, format="vasp", direct=True, sort=True)

def _print_cell(label, atoms):
    print("  {}: a,b,c = {:.4f}, {:.4f}, {:.4f} A   alpha,beta,gamma = {:.3f}, {:.3f}, {:.3f} deg"
          .format(label, *atoms.cell.cellpar()))

def _expand(patterns):
    # Windows shells do not expand wildcards, so do it here
    paths = []
    for pattern in patterns:
        paths.extend(sorted(glob.glob(pattern)) or [pattern])
    return paths

def _fmt_matrix(P):
    return "[" + "; ".join(" ".join(f"{int(x):>3d}" for x in row) for row in P) + "]"


# ---------- orthorhombify ----------
def _orthorhombify_one(path, args):
    atoms = _read(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    cellpar = atoms.cell.cellpar()
    system, uniq = detect_system(atoms, args.angle_tol)
    if args.mode == "monoclinic":
        system = "monoclinic"
        uniq = uniq if uniq is not None else int(np.argmax(np.abs(cellpar[3:] - 90.0)))
    elif args.mode == "triclinic":
        system = "triclinic"

    print("=" * 78)
    print(path)
    print(f"  {len(atoms)} atoms, {atoms.get_chemical_formula()}")
    print("  a,b,c = {:.4f}, {:.4f}, {:.4f} A   alpha,beta,gamma = {:.3f}, {:.3f}, {:.3f} deg"
          .format(*cellpar))
    dmin = min_distance(atoms)
    if dmin < 0.5:
        print(f"  WARNING: shortest interatomic distance is {dmin:.3f} A; the CIF symmetry "
              f"expansion may be wrong (check for duplicate/overlapping atoms).")

    if system in ORTHOGONAL_SYSTEMS:
        print(f"  Cell shape is {system.upper()} (all angles 90 within --angle-tol); "
              "writing the cell unchanged.")
    elif system in ("monoclinic", "hexagonal"):
        print(f"  Cell shape is {system.upper()}, unique axis {AXES[uniq]} "
              f"({ANGLE_NAMES[uniq]} = {cellpar[3 + uniq]:.3f} deg is the only non-90 angle); "
              "using the exact 2D search.")
    else:
        how = "the original HNF search" if args.method == "hnf" else "the 3D lattice-vector search"
        print(f"  Cell shape is {system.upper()}; using {how}.")
    if system not in ORTHOGONAL_SYSTEMS and args.method != "hnf":
        table = best_angle_by_size(atoms, args.max_cells, max_index=args.max_index,
                                   angle_tol=args.angle_tol, mode=args.mode)
        if table:
            print("  Best possible angle deviation by supercell size N (x unit cell):")
            print("    " + "  ".join(f"N={n}: {d:.3f}" for n, d in table.items()))

    results = find_supercells(
        atoms, args.max_cells, args.top, mode=args.mode, method=args.method,
        angle_tol=args.angle_tol, tol=args.tol, anis_max=args.anis_max, max_index=args.max_index,
        unimod=args.unimod, fast=args.fast, gate=args.gate, stop_on_tol=args.stop_on_tol)
    if not results:
        print("  No candidates found. Increase --max-cells / --max-index, or relax --anis-max / --gate.",
              file=sys.stderr)
        return False

    os.makedirs(args.outdir, exist_ok=True)
    for rank, cand in enumerate(results, 1):
        sc = cand.build(atoms, force_orthogonal=args.force_orthogonal)
        out = os.path.join(args.outdir, f"{stem}_orth_{rank}")
        _write_structure(out, sc)
        write_mapping_json(out + ".json", cand.mapping, atoms, sc, source=os.path.basename(path),
                           n_cells=cand.n_cells, max_angle_deviation_deg=round(cand.max_dev, 6),
                           strain=sc.info["orthorhombify_strain"])
        print(f"  [{rank}] N={cand.n_cells} ({len(sc)} atoms)  P={_fmt_matrix(cand.P)}")
        print("       a,b,c = {:.4f}, {:.4f}, {:.4f} A   angles = {:.3f}, {:.3f}, {:.3f} deg   "
              "max dev = {:.3f} deg  anis = {:.2f}".format(*cand.lengths, *cand.angles,
                                                         cand.max_dev, cand.anisotropy))
        if args.force_orthogonal:
            print(f"       forced to exactly 90/90/90; max Green-Lagrange strain = "
                  f"{sc.info['orthorhombify_strain']:.2e}")
        print(cand.mapping.summary(atoms.cell.array, indent="       "))
        print(f"       -> {out}.cif / .POSCAR / .json")
    return True


def orthorhombify_main(argv=None):
    ap = argparse.ArgumentParser(
        prog="orthorhombify",
        description="Find near-orthorhombic integer supercells for monoclinic or triclinic crystals.")
    ap.add_argument("cif", nargs="+", help="One or more input CIF files (wildcards allowed)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    ap.add_argument("--mode", choices=["auto", "monoclinic", "triclinic"], default="auto",
                    help="Force a search type: monoclinic = 2D search around the unique axis, "
                         "triclinic = 3D search. Default: pick from the cell shape.")
    ap.add_argument("--method", choices=["vectors", "hnf"], default="vectors",
                    help="3D search method. vectors (default) is fast and thorough; hnf is the original v1 search.")
    ap.add_argument("--angle-tol", type=float, default=1e-3,
                    help="An angle counts as 90 deg if within this many degrees (detection). Default: 1e-3.")
    ap.add_argument("--max-cells", type=int, default=6,
                    help="Max supercell size N = det(P), i.e. number of unit cells. Default: 6.")
    ap.add_argument("--max-index", type=int, default=None,
                    help="Max |integer coefficient| tried for the new axes. Default: 12 (2D search), 5 (3D search).")
    ap.add_argument("--top", type=int, default=2, help="How many best solutions to export (1..5). Default: 2.")
    ap.add_argument("--tol", type=float, default=None,
                    help="Angle tolerance (deg). [monoclinic] prefer the SMALLEST supercell whose angle "
                         "deviation is <= tol. [triclinic] used with --stop-on-tol.")
    ap.add_argument("--anis-max", type=float, default=None,
                    help="Reject boxes with max/min edge length above this value, e.g. 3.0.")
    ap.add_argument("--force-orthogonal", action="store_true",
                    help="Set the final cell to exactly 90/90/90 (keeps fractional coords; introduces a "
                         "small strain, which is reported).")
    ap.add_argument("--outdir", default=".", help="Output directory. Default: current directory.")
    ap.add_argument("--unimod", type=int, default=2, help="[--method hnf] Unimodular shear coeff limit. Default: 2.")
    ap.add_argument("--fast", action="store_true", help="[--method hnf] Only diagonal H.")
    ap.add_argument("--gate", type=float, default=None, help="[--method hnf] Angle-sum gate (deg), e.g. 12.")
    ap.add_argument("--stop-on-tol", action="store_true", help="[--method hnf] Stop once sum dev <= 3*tol.")
    args = ap.parse_args(argv)

    if not 1 <= args.top <= 5:
        ap.error("--top must be between 1 and 5.")
    ok = [_orthorhombify_one(p, args) for p in _expand(args.cif)]
    return 0 if all(ok) else 2


# ---------- orthorhombify-transform ----------
def _transform_one(path, args):
    import json

    atoms = _read(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    print("=" * 78)
    print(path)
    print(f"  {len(atoms)} atoms, {atoms.get_chemical_formula()}")
    _print_cell("input ", atoms)
    try:
        if args.inverse:
            out = back_to_original(atoms, args.inverse, replicate=args.replicate,
                                   reduce=args.reduce, tol=args.tol)
        else:
            out = transform(atoms, parse_matrix(args.matrix), reduce=args.reduce, tol=args.tol)
    except ValueError as err:
        print(f"  ERROR: {err}", file=sys.stderr)
        return False

    mapping = LatticeMapping(out.info["orthorhombify_P"])
    _print_cell("output", out)
    print(f"  {len(out)} atoms in the new cell (det = {mapping.det}); method: {out.info['orthorhombify_reduce']}"
          + (f", copies of each atom agree within {out.info['orthorhombify_spread']:.3f} A"
             if "orthorhombify_spread" in out.info else ""))
    if "orthorhombify_warning" in out.info:
        print(f"  WARNING: {out.info['orthorhombify_warning']}")
    if args.inverse:
        reps = out.info["orthorhombify_replicate"]
        with open(args.inverse) as fd:
            ref = json.load(fd)["old_cell"]
        if tuple(reps) != (1, 1, 1):
            print(f"  undid LAMMPS replicate {reps[0]} {reps[1]} {reps[2]}")
        names = ("a", "b", "c", "alpha", "beta", "gamma")
        print("  change vs the original cell in " + os.path.basename(args.inverse) + ":")
        print("    " + "   ".join(f"{n} {ref[n]:.4f} -> {v:.4f} ({100 * (v - ref[n]) / ref[n]:+.2f}%)"
                                  for n, v in zip(names, out.cell.cellpar())))
    print(mapping.summary(atoms.cell.array, indent="  "))

    os.makedirs(args.outdir, exist_ok=True)
    out_stem = os.path.join(args.outdir, stem + (args.suffix or ("_original" if args.inverse else "_transformed")))
    _write_structure(out_stem, out)
    write_mapping_json(out_stem + ".json", mapping, atoms, out, source=os.path.basename(path),
                       method=out.info["orthorhombify_reduce"])
    print(f"  -> {out_stem}.cif / .POSCAR / .json")
    return True


def transform_main(argv=None):
    ap = argparse.ArgumentParser(
        prog="orthorhombify-transform",
        description="Apply any cell transformation (bigger or smaller cell), or undo an orthorhombify "
                    "step to get back to the original cell (for example after a LAMMPS run).")
    ap.add_argument("files", nargs="+", help="Structure files (CIF, POSCAR, LAMMPS .data, ...); wildcards allowed")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    how = ap.add_mutually_exclusive_group(required=True)
    how.add_argument("--matrix", metavar="'P'",
                     help="new_cell = P @ old_cell. Rows separated by ';', fractions allowed, e.g. "
                          "'1 0 1; 0 1 0; -4 0 3' or '1/2 -1/2 0; 1/2 1/2 0; 0 0 1'. Three numbers = diagonal.")
    how.add_argument("--inverse", metavar="JSON",
                     help="Go back to the original cell using the .json written by orthorhombify.")
    ap.add_argument("--replicate", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                    help="[--inverse] LAMMPS replicate used on the supercell. Default: worked out from the box.")
    ap.add_argument("--reduce", choices=["auto", "average", "dedup"], default="auto",
                    help="How to fill a smaller cell: average the copies of each atom (needs Orthorhombify "
                         "atom order, e.g. after LAMMPS), or merge folded atoms closer than --tol. Default: auto.")
    ap.add_argument("--tol", type=float, default=0.5, help="Merge / spread tolerance (A). Default: 0.5")
    ap.add_argument("--outdir", default=".", help="Output directory. Default: current directory.")
    ap.add_argument("--suffix", default=None,
                    help="Added to the file name. Default: _original (--inverse) or _transformed.")
    args = ap.parse_args(argv)
    ok = [_transform_one(p, args) for p in _expand(args.files)]
    return 0 if all(ok) else 2


# ---------- orthorhombify-lammps ----------
def _parse_charges(items):
    charges = {}
    for item in items or []:
        try:
            el, val = item.split("=")
            charges[el.strip()] = float(val)
        except ValueError:
            raise SystemExit(f"Bad --charge '{item}', expected e.g. O=-0.3668")
    return charges

def _lammps_one(path, charges, args):
    atoms = _read(path)
    name = os.path.splitext(os.path.basename(path))[0]
    outdir = os.path.join(args.outdir, name)
    try:
        files = write_lammps_inputs(
            atoms, outdir, name, charges=charges, neutralize=args.neutralize, ffield=args.ffield,
            lgvdw=args.lgvdw, temp=args.temp, press=args.press, timestep=args.timestep,
            nvt_steps=args.nvt_steps, npt_steps=args.npt_steps, barostat=args.barostat,
            minimize=args.minimize, replicate=args.replicate, seed=args.seed)
    except KeyError as err:
        print(f"  {path}: {err.args[0]}; use --charge El=value", file=sys.stderr)
        return False

    info = files.charges
    counts = {el: atoms.get_chemical_symbols().count(el) for el in files.types}
    print("=" * 78)
    print(path)
    print(f"  {len(atoms)} atoms  " + "  ".join(f"{el}:{n}" for el, n in counts.items()))
    print(f"  net charge with initial charges = {info.net_before:+.4f} e")
    if info.shift:
        print(f"  neutralized: shifted every charge by {info.shift:+.6f} e -> "
              + ", ".join(f"{el} = {v:+.4f}" for el, v in sorted(info.per_element.items())))
    elif abs(info.net_before) > 1e-6:
        print(f"  WARNING: --no-neutralize is set, so the system keeps a net charge of {info.net_before:+.4f} e.")

    box = "orthogonal" if is_orthogonal(atoms) else \
        "triclinic (small tilt, angles {:.3f} {:.3f} {:.3f})".format(*atoms.cell.angles())
    print("  atom types: " + "  ".join(f"{i}={el}" for i, el in enumerate(files.types, 1)) + f"   box: {box}")
    need = suggest_replicate(atoms, args.min_length, args.replicate)
    if need:
        lengths = atoms.cell.lengths() * np.array(args.replicate or (1, 1, 1))
        print(f"  NOTE: simulated box edges {np.round(lengths, 2).tolist()} A; shortest is below "
              f"{args.min_length} A (QEq/ReaxFF cutoff is 10 A). Consider --replicate {need[0]} {need[1]} {need[2]}")
    print(f"  -> {files.data_file}  and  {os.path.basename(files.input_file)}")
    return True


def lammps_main(argv=None):
    ap = argparse.ArgumentParser(
        prog="orthorhombify-lammps",
        description="Make LAMMPS data + input files (with initial charges) from CIF/POSCAR structures.")
    ap.add_argument("files", nargs="+", help="Structure files (CIF, POSCAR, ...); wildcards allowed")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
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
    args = ap.parse_args(argv)

    charges = _parse_charges(args.charge)
    ok = [_lammps_one(p, charges, args) for p in _expand(args.files)]
    return 0 if all(ok) else 2
