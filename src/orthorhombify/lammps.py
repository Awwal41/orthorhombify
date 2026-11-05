"""Write LAMMPS data files with per-element starting charges and a ReaxFF + QEq
equilibration input script (minimize -> NVT -> NPT)."""
from __future__ import annotations

import os
from collections import Counter
from dataclasses import dataclass

import numpy as np
from ase import Atoms
from ase.io.lammpsdata import write_lammps_data

DEFAULT_CHARGES = {"O": -0.3668, "N": 0.3668, "C": -0.4320, "H": 0.2160}
#: LAMMPS atom types are assigned in this order (only elements present are used).
TYPE_ORDER = ["C", "H", "N", "O"]


@dataclass
class ChargeInfo:
    """Charges assigned by :func:`assign_charges`.

    Attributes:
        charges: Per-atom charges (e), in atom order.
        net_before: Total charge before neutralizing.
        shift: Amount added to every charge (0 if not neutralized).
        per_element: Final charge of each element.
    """

    charges: np.ndarray
    net_before: float
    shift: float
    per_element: dict


@dataclass
class LammpsFiles:
    """Paths and charge details returned by :func:`write_lammps_inputs`."""

    data_file: str
    input_file: str
    types: list
    charges: ChargeInfo


def species_order(symbols) -> list[str]:
    """Element order used for LAMMPS atom types: C, H, N, O first, then the rest."""
    present = set(symbols)
    return [e for e in TYPE_ORDER if e in present] + sorted(present - set(TYPE_ORDER))


def is_orthogonal(atoms: Atoms, tol: float = 1e-6) -> bool:
    """True if the cell is exactly 90/90/90 and aligned with x, y, z."""
    cell = atoms.cell.array
    return bool(abs(cell[1, 0]) < tol and abs(cell[2, 0]) < tol and abs(cell[2, 1]) < tol
                and np.allclose(atoms.cell.angles(), 90.0, atol=1e-6))


def assign_charges(atoms: Atoms, charges: dict | None = None, neutralize: bool = True) -> ChargeInfo:
    """Per-atom charges from per-element values.

    Args:
        atoms: The structure.
        charges: ``{element: charge}``. Missing entries fall back to
            :data:`DEFAULT_CHARGES`.
        neutralize: Shift every charge by the same amount so the total is zero.

    Raises:
        KeyError: An element in ``atoms`` has no charge.
    """
    table = dict(DEFAULT_CHARGES)
    table.update(charges or {})
    symbols = atoms.get_chemical_symbols()
    missing = sorted(set(symbols) - set(table))
    if missing:
        raise KeyError(f"No charge given for {missing}")
    q = np.array([table[s] for s in symbols], dtype=float)
    net = float(q.sum())
    shift = -net / len(q) if (neutralize and abs(net) > 1e-12) else 0.0
    q = q + shift
    per_element = {el: table[el] + shift for el in Counter(symbols)}
    return ChargeInfo(charges=q, net_before=net, shift=shift, per_element=per_element)


def suggest_replicate(atoms: Atoms, min_length: float = 20.0, replicate=None):
    """Replication ``(nx, ny, nz)`` that makes every box edge >= ``min_length``.

    Returns ``None`` if the box (times ``replicate``, if given) is already big enough.
    """
    reps = np.array(replicate or (1, 1, 1))
    lengths = atoms.cell.lengths()
    if (lengths * reps).min() >= min_length:
        return None
    return tuple(int(np.ceil(min_length / L)) for L in lengths)


_INPUT_TEMPLATE = """\
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


def make_input_script(atoms: Atoms, name: str, data_file: str, *, ffield: str = "ffield.reax",
                      lgvdw: bool = False, temp: float = 300.0, press: float = 1.0,
                      timestep: float = 0.1, nvt_steps: int = 10000, npt_steps: int = 50000,
                      barostat: str = "aniso", minimize: bool = True, replicate=None,
                      seed: int = 4928459) -> str:
    """Text of the LAMMPS input script (see :func:`write_lammps_inputs` for arguments)."""
    order = species_order(atoms.get_chemical_symbols())
    if barostat not in ("aniso", "iso", "tri"):
        raise ValueError("barostat must be 'aniso', 'iso' or 'tri'")
    if barostat == "tri" and is_orthogonal(atoms):
        barostat = "aniso"  # LAMMPS 'tri' needs a triclinic box
    groups = "\n".join(f"group           g{el} type {i}\n"
                       f"compute         q{el} g{el} reduce ave c_q"
                       for i, el in enumerate(order, 1))
    return _INPUT_TEMPLATE.format(
        name=name, natoms=len(atoms), formula=atoms.get_chemical_formula(),
        data=os.path.basename(data_file), ffield=ffield, temp=temp, press=press,
        dt=timestep, nvt_steps=nvt_steps, npt_steps=npt_steps, seed=seed,
        lgvdw=" lgvdw yes" if lgvdw else "", elements=" ".join(order), groups=groups,
        qcols=" ".join(f"c_q{el}" for el in order),
        minimize=("min_style       cg\nminimize        1.0e-6 1.0e-8 2000 20000\n"
                  if minimize else "# (skipped)\n"),
        baro=barostat,
        replicate="replicate       {} {} {}\n".format(*replicate) if replicate else "",
        print_q='print           "Final average charges: '
                + "  ".join(f"{el} = $(c_q{el}:%.4f)" for el in order) + '"')


def write_lammps_inputs(atoms: Atoms, outdir: str, name: str, *, charges: dict | None = None,
                        neutralize: bool = True, ffield: str = "ffield.reax", lgvdw: bool = False,
                        temp: float = 300.0, press: float = 1.0, timestep: float = 0.1,
                        nvt_steps: int = 10000, npt_steps: int = 50000, barostat: str = "aniso",
                        minimize: bool = True, replicate=None, seed: int = 4928459) -> LammpsFiles:
    """Write ``<outdir>/<name>.data`` and ``<outdir>/in.<name>``.

    Args:
        atoms: The structure (any cell; LAMMPS tilt factors are written if needed).
        outdir: Folder to write into (created if missing).
        name: Base name of the files.
        charges: ``{element: charge}`` overrides; defaults are :data:`DEFAULT_CHARGES`.
        neutralize: Shift every charge by the same amount so the total is zero.
        ffield: ReaxFF force-field file name used in ``pair_coeff``.
        lgvdw: Force field is ReaxFF-lg (adds ``lgvdw yes``).
        temp: Temperature (K).
        press: Pressure (atm).
        timestep: Timestep (fs).
        nvt_steps: Number of NVT steps.
        npt_steps: Number of NPT steps.
        barostat: ``"aniso"`` (keeps angles), ``"iso"``, or ``"tri"`` (relaxes angles).
        minimize: Run an energy minimization first.
        replicate: ``(nx, ny, nz)`` to add a ``replicate`` command after ``read_data``.
        seed: Random seed for the starting velocities.
    """
    info = assign_charges(atoms, charges, neutralize)
    atoms = atoms.copy()
    atoms.set_initial_charges(info.charges)
    order = species_order(atoms.get_chemical_symbols())

    os.makedirs(outdir, exist_ok=True)
    data_file = os.path.join(outdir, f"{name}.data")
    input_file = os.path.join(outdir, f"in.{name}")
    with open(data_file, "w") as fd:
        write_lammps_data(fd, atoms, specorder=order, atom_style="charge", masses=True, units="real")
    with open(input_file, "w") as fd:
        fd.write(make_input_script(
            atoms, name, data_file, ffield=ffield, lgvdw=lgvdw, temp=temp, press=press,
            timestep=timestep, nvt_steps=nvt_steps, npt_steps=npt_steps, barostat=barostat,
            minimize=minimize, replicate=replicate, seed=seed))
    return LammpsFiles(data_file=data_file, input_file=input_file, types=order, charges=info)
