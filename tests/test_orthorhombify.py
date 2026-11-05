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
    assert orth.detect_system([5, 5, 7, 90, 90, 90]) == ("tetragonal", None)
    assert orth.detect_system([5, 5, 5, 90, 90, 90]) == ("cubic", None)
    assert orth.detect_system([5, 6, 7, 90, 102, 90]) == ("monoclinic", 1)
    assert orth.detect_system([5, 6, 7, 90, 90, 116]) == ("monoclinic", 2)
    assert orth.detect_system([3, 3, 5, 90, 90, 120]) == ("hexagonal", 2)
    assert orth.detect_system([5, 5, 5, 70, 70, 70]) == ("rhombohedral", None)
    assert orth.detect_system([5, 6, 7, 95, 103, 88])[0] == "triclinic"


@pytest.mark.parametrize("cellpar, n_cells", [([3, 3, 5, 90, 90, 120], 2),
                                              ([3, 3, 5, 90, 90, 60], 2),
                                              ([5, 5, 5, 70, 70, 70], 6)])
def test_hexagonal_and_rhombohedral_become_exactly_orthorhombic(cellpar, n_cells):
    best = orth.find_supercells(cellpar, max_cells=8, top=1)[0]
    assert best.n_cells == n_cells and best.max_dev < 1e-6
    assert np.linalg.det(best.P) > 0


def test_vector_search_beats_hnf_on_triclinic():
    cp = [7.1, 8.3, 9.2, 101, 96, 105]
    new = orth.find_supercells(cp, max_cells=8, top=1)[0]
    old = orth.find_supercells(cp, max_cells=8, top=1, method="hnf")[0]
    assert new.max_dev < old.max_dev
    assert np.linalg.det(new.P) > 0


def test_lattice_mapping_hmx():
    m = orth.LatticeMapping([[1, 0, 1], [0, 1, 0], [-4, 0, 3]])
    assert m.det == 7
    assert m.new_axis_expressions() == ["a + c", "b", "-4a + 3c"]
    assert m.old_axis_expressions()[0] == "(3/7)a' - (1/7)c'"
    assert m.plane_to_old([0, 0, 1]) == [-1, 0, 1]
    # the new (001) plane contains the new a and b axes
    for axis in m.new_axes_as_old_directions()[:2]:
        assert np.dot(m.plane_to_old([0, 0, 1]), axis) == 0
    # directions and planes round-trip, and the zone law h.u is preserved
    u, h = [1, 1, 0], [1, 2, 3]
    assert m.direction_to_old(m.direction_to_new(u)) == u
    assert np.dot(m.plane_to_new(h), m.direction_to_new([0, 0, 1])) != 0
    assert m.inverse().inverse().P == m.P


def test_transform_and_back_exact(hmx):
    best = orth.find_supercells(hmx, max_cells=8, top=1)[0]
    sc = best.build(hmx)
    back = orth.transform(sc, best.mapping.P_inv)
    assert back.info["orthorhombify_reduce"] == "average"
    assert len(back) == len(hmx)
    assert np.allclose(back.cell.cellpar(), hmx.cell.cellpar(), atol=1e-6)
    dedup = orth.transform(sc, best.mapping.P_inv, reduce="dedup", tol=0.1)
    assert len(dedup) == len(hmx)


def test_parse_matrix():
    assert orth.parse_matrix("2 2 1") == [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
    m = orth.parse_matrix("1/2 -1/2 0; 1/2 1/2 0; 0 0 1")
    assert str(m[0][0]) == "1/2"


def test_transform_accepts_text():
    from ase.build import bulk
    ortho = orth.transform(bulk("Mg", "hcp", a=3.21, c=5.21), "1 0 0; 1 2 0; 0 0 1")
    hexa = orth.transform(ortho, "1 0 0; -1/2 1/2 0; 0 0 1")
    assert len(ortho) == 4 and len(hexa) == 2
    assert np.allclose(hexa.cell.cellpar(), [3.21, 3.21, 5.21, 90, 90, 120])


def test_back_to_original_after_md(tmp_path, hmx):
    """Simulate LAMMPS: replicate 2x1x1, strain the box, add thermal noise, write data."""
    assert orthorhombify_main([orth.example_path(), "--max-cells", "8", "--top", "1",
                               "--outdir", str(tmp_path)]) == 0
    sc = read(tmp_path / "beta_HMX_orth_1.cif")
    md = sc.repeat((2, 1, 1))
    md.set_cell(md.cell.array * np.array([1.01, 0.99, 1.02])[:, None], scale_atoms=True)
    md.positions += np.random.default_rng(0).normal(0, 0.05, md.positions.shape)
    orth.write_lammps_inputs(md, str(tmp_path), "md")

    md_read = read(tmp_path / "md.data", format="lammps-data", atom_style="charge", units="real")
    back = orth.back_to_original(md_read, str(tmp_path / "beta_HMX_orth_1.json"))
    assert back.info["orthorhombify_replicate"] == (2, 1, 1)
    assert len(back) == len(hmx)
    assert back.get_chemical_formula() == hmx.get_chemical_formula()
    assert back.cell.cellpar()[1] == pytest.approx(hmx.cell.cellpar()[1] * 0.99, rel=1e-6)
    assert back.info["orthorhombify_spread"] < 0.5

    from orthorhombify.cli import transform_main
    assert transform_main([str(tmp_path / "md.data"), "--inverse", str(tmp_path / "beta_HMX_orth_1.json"),
                           "--outdir", str(tmp_path)]) == 0
    assert (tmp_path / "md_original.cif").is_file()


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
    res = orth.find_supercells(atoms, max_cells=3, top=1, method="hnf", unimod=1, gate=15)
    assert res and res[0].sum_dev < 15
    vec = orth.find_supercells(atoms, max_cells=3, top=1)[0]
    assert vec.max_dev <= res[0].max_dev + 1e-9 and vec.n_cells <= 3


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
