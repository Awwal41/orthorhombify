import numpy as np
import pytest
from ase.build import bulk
from ase.cell import Cell
from ase.io import read

import orthorhombify as orth
from orthorhombify.cli import lammps_main, orthorhombify_main


@pytest.fixture(scope="module")
def hmx():
    return read(orth.example_path("beta_HMX"))


def test_detect_system():
    assert orth.detect_system([5, 6, 7, 90, 90, 90]) == ("orthorhombic", None)
    assert orth.detect_system([5, 6, 7, 90, 102, 90]) == ("monoclinic", 1)
    assert orth.detect_system([5, 6, 7, 90, 90, 116]) == ("monoclinic", 2)
    assert orth.detect_system([5, 6, 7, 95, 103, 88])[0] == "triclinic"


def test_beta_hmx_monoclinic(hmx):
    best = orth.find_supercells(hmx, max_cells=8, top=1)[0]
    assert best.n_cells == 7
    assert best.max_dev < 0.01
    assert best.angles[0] == pytest.approx(90.0) and best.angles[2] == pytest.approx(90.0)

    sc = best.build(hmx)
    assert len(sc) == 7 * len(hmx)
    assert sc.get_volume() == pytest.approx(7 * hmx.get_volume())
    assert orth.min_distance(sc) == pytest.approx(orth.min_distance(hmx), abs=1e-6)


def test_tol_prefers_smallest_box(hmx):
    best = orth.find_supercells(hmx, max_cells=8, top=1, tol=1.0)[0]
    assert best.n_cells == 4 and best.max_dev <= 1.0


def test_best_angle_by_size_is_non_increasing(hmx):
    table = orth.best_angle_by_size(hmx, max_cells=8)
    values = [table[n] for n in sorted(table)]
    assert values == sorted(values, reverse=True)
    assert table[7] < 0.01


def test_unique_axis_c():
    atoms = bulk("Cu", "sc", a=3.0)
    atoms.set_cell(Cell.fromcellpar([9.936, 12.61, 12.277, 90, 90, 116.66]), scale_atoms=True)
    best = orth.find_supercells(atoms, max_cells=8, top=1)[0]
    assert best.n_cells == 7 and best.max_dev < 0.2


def test_triclinic_runs():
    atoms = bulk("Cu", "sc", a=3.0)
    atoms.set_cell(Cell.fromcellpar([6.5, 11.0, 8.7, 95, 103, 88]), scale_atoms=True)
    res = orth.find_supercells(atoms, max_cells=3, top=1, unimod=1, gate=15)
    assert res and res[0].sum_dev < 15


def test_force_orthogonal(hmx):
    sc = orth.find_supercells(hmx, max_cells=8, top=1)[0].build(hmx, force_orthogonal=True)
    assert np.allclose(sc.cell.angles(), 90.0)
    assert sc.info["orthorhombify_strain"] < 1e-4


def test_neutralize():
    atoms = bulk("Cu", "sc", a=3.0) * (2, 1, 1)
    atoms.set_chemical_symbols(["C", "O"])
    info = orth.assign_charges(atoms)
    assert info.net_before == pytest.approx(-0.4320 - 0.3668)
    assert info.charges.sum() == pytest.approx(0.0)
    raw = orth.assign_charges(atoms, neutralize=False)
    assert raw.shift == 0.0 and raw.charges.sum() == pytest.approx(info.net_before)
    with pytest.raises(KeyError):
        atoms.set_chemical_symbols(["C", "S"])
        orth.assign_charges(atoms)


def test_write_lammps_inputs(tmp_path, hmx):
    files = orth.write_lammps_inputs(hmx, str(tmp_path), "hmx", replicate=(3, 2, 2))
    assert files.types == ["C", "H", "N", "O"]
    data = open(files.data_file).read()
    assert "Atoms # charge" in data and "4 atom types" in data
    script = open(files.input_file).read()
    assert "replicate       3 2 2" in script
    assert "pair_coeff      * * ${ffield} C H N O" in script
    assert "fix             qeq all qeq/reaxff" in script


def test_cli_end_to_end(tmp_path):
    out = tmp_path / "output"
    assert orthorhombify_main([orth.example_path(), "--max-cells", "8", "--top", "1",
                               "--outdir", str(out)]) == 0
    cif = out / "beta_HMX_orth_1.cif"
    assert len(read(cif)) == 392
    assert lammps_main([str(cif), "--outdir", str(tmp_path / "lammps")]) == 0
    assert (tmp_path / "lammps" / "beta_HMX_orth_1" / "in.beta_HMX_orth_1").is_file()
