"""Find near-orthorhombic integer supercells of monoclinic and triclinic crystals.

Rows convention throughout: cell rows are lattice vectors and a supercell is
``new_rows = P @ old_rows`` with an integer matrix ``P``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import permutations

import numpy as np
from ase import Atoms
from ase.build import make_supercell
from ase.cell import Cell
from ase.neighborlist import neighbor_list

AXES = "abc"
ANGLE_NAMES = ("alpha", "beta", "gamma")


@dataclass(frozen=True, eq=False)
class Candidate:
    """One supercell found by :func:`find_supercells`.

    Attributes:
        P: 3x3 integer matrix, ``new_cell_rows = P @ old_cell_rows``.
        lengths: New cell lengths a, b, c (Å).
        angles: New cell angles alpha, beta, gamma (degrees).
        n_cells: Number of unit cells in the supercell, ``|det(P)|``.
        max_dev: Largest deviation of any angle from 90 degrees.
        sum_dev: Sum of the three angle deviations from 90 degrees.
        anisotropy: Longest edge divided by the shortest edge.
        volume: Supercell volume (Å^3).
    """

    P: np.ndarray
    lengths: tuple
    angles: tuple
    n_cells: int
    max_dev: float
    sum_dev: float
    anisotropy: float
    volume: float

    def build(self, atoms: Atoms, force_orthogonal: bool = False) -> Atoms:
        """Build this supercell from ``atoms``. See :func:`build_supercell`."""
        return build_supercell(atoms, self.P, force_orthogonal=force_orthogonal)

    @property
    def mapping(self):
        """:class:`~orthorhombify.LatticeMapping` between the input cell and this supercell."""
        from .mapping import LatticeMapping
        return LatticeMapping(self.P)


# ---------- geometry on ROWS ----------
def _angle_deg(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return 0.0
    c = max(-1.0, min(1.0, float(np.dot(u, v)) / (nu * nv)))
    return math.degrees(math.acos(c))

def _lengths_angles(B):
    a, b, c = B
    lengths = (float(np.linalg.norm(a)), float(np.linalg.norm(b)), float(np.linalg.norm(c)))
    angles = (_angle_deg(b, c), _angle_deg(a, c), _angle_deg(a, b))
    return lengths, angles

def _candidate(P, A_rows):
    P = np.asarray(P, dtype=int)
    B = P @ A_rows
    lengths, angles = _lengths_angles(B)
    devs = [abs(x - 90.0) for x in angles]
    return Candidate(
        P=P, lengths=lengths, angles=angles,
        n_cells=int(round(abs(np.linalg.det(P)))),
        max_dev=max(devs), sum_dev=sum(devs),
        anisotropy=max(lengths) / max(1e-12, min(lengths)),
        volume=float(abs(np.linalg.det(B))),
    )

def _cell_rows(obj):
    if isinstance(obj, Atoms):
        return obj.cell.array.copy()
    return np.array(Cell.new(obj).array)


# ---------- crystal-system detection ----------
ORTHOGONAL_SYSTEMS = ("cubic", "tetragonal", "orthorhombic")
PLANE_SYSTEMS = ("monoclinic", "hexagonal")      # one angle != 90: exact 2D search
SPACE_SYSTEMS = ("rhombohedral", "triclinic")    # 3D search


def detect_system(cell, angle_tol: float = 1e-3, length_tol: float = 1e-4):
    """Classify a cell by its shape (lengths and angles, not by symmetry).

    Args:
        cell: ASE ``Atoms``, a 3x3 cell matrix, or ``[a, b, c, alpha, beta, gamma]``.
        angle_tol: Two angles are equal if within this many degrees.
        length_tol: Two lengths are equal if within this relative difference.

    Returns:
        ``(system, unique_axis)``. ``system`` is one of ``"cubic"``, ``"tetragonal"``,
        ``"orthorhombic"``, ``"hexagonal"`` (two equal edges at 120 or 60 degrees),
        ``"monoclinic"``, ``"rhombohedral"`` (equal edges, three equal angles) or
        ``"triclinic"``. ``unique_axis`` is 0, 1 or 2 (a, b, c) for hexagonal and
        monoclinic cells (the axis at 90 degrees to the other two), otherwise ``None``.
    """
    a, b, c, al, be, ga = Cell(_cell_rows(cell)).cellpar()
    lengths, angles = (a, b, c), (al, be, ga)
    off = [abs(x - 90.0) > angle_tol for x in angles]

    def same(x, y):
        return abs(x - y) <= length_tol * max(x, y)

    if sum(off) == 0:
        if same(a, b) and same(b, c):
            return "cubic", None
        if same(a, b) or same(b, c) or same(a, c):
            return "tetragonal", None
        return "orthorhombic", None
    if sum(off) == 1:
        # alpha = angle(b,c) -> unique a ; beta -> unique b ; gamma -> unique c
        k = off.index(True)
        i, j = [m for m in range(3) if m != k]
        if same(lengths[i], lengths[j]) and min(abs(angles[k] - 120), abs(angles[k] - 60)) <= angle_tol:
            return "hexagonal", k
        return "monoclinic", k
    if same(a, b) and same(b, c) and abs(al - be) <= angle_tol and abs(be - ga) <= angle_tol:
        return "rhombohedral", None
    return "triclinic", None


# ---------- MONOCLINIC: exact 2D search perpendicular to the unique axis ----------
def _mono_pairs(A_rows, uniq, max_cells, max_index, anis_max):
    """All in-plane vector pairs (v1, v2) with 1 <= |det| <= max_cells."""
    i1, i2 = [k for k in range(3) if k != uniq]
    p, q, u = A_rows[i1], A_rows[i2], A_rows[uniq]
    lu = float(np.linalg.norm(u))

    R = max_index
    coeffs = [(i, j) for i in range(0, R + 1) for j in range(-R, R + 1)
              if (i, j) != (0, 0) and math.gcd(i, abs(j)) == 1 and not (i == 0 and j < 0)]
    C = np.array(coeffs, dtype=int)
    V = C[:, :1] * p + C[:, 1:] * q
    L = np.linalg.norm(V, axis=1)

    det = np.outer(C[:, 0], C[:, 1]) - np.outer(C[:, 1], C[:, 0])
    cosm = (V @ V.T) / np.outer(L, L)
    dev = np.degrees(np.arcsin(np.clip(np.abs(cosm), 0.0, 1.0)))

    mask = (np.abs(det) >= 1) & (np.abs(det) <= max_cells)
    mask &= np.triu(np.ones_like(mask, dtype=bool), k=1)
    if anis_max is not None:
        Lmax = np.maximum(np.maximum.outer(L, L), lu)
        Lmin = np.minimum(np.minimum.outer(L, L), lu)
        mask &= (Lmax / Lmin) <= anis_max

    ii, jj = np.nonzero(mask)
    return C, L, ii, jj, np.abs(det[ii, jj]), dev[ii, jj], (i1, i2)

def _search_monoclinic(A_rows, uniq, max_cells, max_index, top, anis_max, tol):
    C, L, ii, jj, N, d, (i1, i2) = _mono_pairs(A_rows, uniq, max_cells, max_index, anis_max)
    if len(ii) == 0:
        return []
    if tol is not None:
        # smallest supercell that meets the tolerance first, then best angle
        order = np.lexsort((d, N, d > tol))
    else:
        order = np.lexsort((L[ii] + L[jj], N, np.round(d, 6)))

    results, seen = [], set()
    for k in order:
        ca, cb = C[ii[k]].copy(), C[jj[k]].copy()
        # shorter in-plane vector first, then make the transform right-handed
        if L[jj[k]] < L[ii[k]]:
            ca, cb = cb, ca
        if ca[0] * cb[1] - ca[1] * cb[0] < 0:
            cb = -cb
        P = np.zeros((3, 3), dtype=int)
        P[uniq, uniq] = 1
        P[i1, i1], P[i1, i2] = ca
        P[i2, i1], P[i2, i2] = cb
        cand = _candidate(P, A_rows)
        key = (round(float(d[k]), 4), cand.n_cells, tuple(sorted(round(x, 3) for x in cand.lengths)))
        if key in seen:
            continue
        seen.add(key)
        results.append(cand)
        if len(results) >= top:
            break
    return results


# ---------- TRICLINIC / RHOMBOHEDRAL: 3D search over short lattice vectors ----------
def _space_pool(A_rows, max_cells, max_index, anis_max, max_pairs=20000, pair_gate=20.0):
    """Candidate triples (v1, v2, v3) of lattice vectors with 1 <= |det| <= max_cells.

    Every primitive lattice vector with coefficients up to ``max_index`` is tried;
    the most perpendicular pairs are completed with the best third vector.
    Returns (C, pool) with pool rows (max_dev, N, length_sum, i, j, k).
    """
    # A box of three nearly perpendicular vectors has volume >= ~0.5 * L1*L2*L3, and no
    # lattice vector is shorter than the shortest one, so long pairs cannot fit max_cells.
    vmax = max_cells * abs(np.linalg.det(A_rows)) / 0.5
    R = max_index
    C = np.array([c for c in np.ndindex(2 * R + 1, 2 * R + 1, 2 * R + 1)]) - R
    C = C[np.any(C != 0, axis=1)]
    first = C[np.arange(len(C)), np.argmax(C != 0, axis=1)]
    C = C[(first > 0) & (np.gcd.reduce(np.abs(C), axis=1) == 1)]
    V = C @ A_rows
    L = np.linalg.norm(V, axis=1)
    keep = L * L.min() ** 2 <= vmax
    C, V, L = C[keep], V[keep], L[keep]
    cosm = np.abs((V / L[:, None]) @ (V / L[:, None]).T)
    np.fill_diagonal(cosm, 1.0)

    fits = np.outer(L, L) * L.min() <= vmax
    i, j = np.nonzero(np.triu((cosm < math.sin(math.radians(pair_gate))) & fits, k=1))
    order = np.lexsort((L[i] + L[j], np.round(cosm[i, j], 9)))[:max_pairs]
    i, j = i[order], j[order]
    rows = []
    for s in range(0, len(i), 512):
        a, b = i[s:s + 512], j[s:s + 512]
        det = np.abs(np.cross(C[a], C[b]) @ C.T)                    # (pairs, vectors)
        ok = (det >= 1) & (det <= max_cells)
        if anis_max is not None:
            hi = np.maximum(np.maximum(L[a], L[b])[:, None], L[None, :])
            lo = np.minimum(np.minimum(L[a], L[b])[:, None], L[None, :])
            ok &= hi / lo <= anis_max
        worst = np.maximum(np.maximum(cosm[a], cosm[b]), cosm[a, b][:, None])
        dev = np.degrees(np.arcsin(np.clip(worst, 0, 1)))
        dev[~ok] = np.inf
        k = np.argsort(dev, axis=1)[:, :3]                         # best 3 third vectors per pair
        rr = np.repeat(np.arange(len(a)), k.shape[1])
        kk = k.ravel()
        d = dev[rr, kk]
        good = np.isfinite(d)
        rr, kk, d = rr[good], kk[good], d[good]
        rows.append(np.column_stack([d, det[rr, kk], L[a[rr]] + L[b[rr]] + L[kk], a[rr], b[rr], kk]))
    return C, np.vstack(rows) if rows else np.zeros((0, 6))

def _aligned(P_rows, A_rows):
    """Order and sign the three new axes so each is as parallel as possible to the old
    axis of the same name, keeping the transformation right-handed."""
    V = P_rows @ A_rows
    U = V / np.linalg.norm(V, axis=1)[:, None]
    W = A_rows / np.linalg.norm(A_rows, axis=1)[:, None]
    perm = max(permutations(range(3)), key=lambda p: sum(abs(U[p[i]] @ W[i]) for i in range(3)))
    P = P_rows[list(perm)].copy()
    cos = np.array([(P[i] @ A_rows) @ W[i] / np.linalg.norm(P[i] @ A_rows) for i in range(3)])
    P[cos < 0] *= -1
    if np.linalg.det(P) < 0:
        P[int(np.argmin(np.abs(cos)))] *= -1
    return P

def _search_space(A_rows, max_cells, max_index, top, anis_max, tol):
    C, pool = _space_pool(A_rows, max_cells, max_index, anis_max)
    if not len(pool):
        return []
    dev, N, Ls = pool[:, 0], pool[:, 1], pool[:, 2]
    if tol is not None:
        order = np.lexsort((Ls, dev, N, dev > tol))
    else:
        order = np.lexsort((Ls, N, np.round(dev, 6)))
    results, seen = [], set()
    for r in order:
        idx = pool[r, 3:].astype(int)
        cand = _candidate(_aligned(C[idx], A_rows), A_rows)
        key = (round(cand.max_dev, 4), cand.n_cells, tuple(sorted(round(x, 3) for x in cand.lengths)))
        if key in seen:
            continue
        seen.add(key)
        results.append(cand)
        if len(results) >= top:
            break
    return results


# ---------- TRICLINIC: 3D HNF x unimodular search (v1 algorithm) ----------
def _small_unimodular_basis(max_coeff=1):
    mats, seen = [], set()
    def add(M):
        key = tuple(M.flatten())
        if key not in seen:
            seen.add(key)
            mats.append(M)

    add(np.eye(3, dtype=int))
    for p in permutations(range(3)):
        P = np.zeros((3, 3), dtype=int)
        for i, j in enumerate(p):
            P[i, j] = 1
        if round(np.linalg.det(P)) == 1:
            add(P)
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            for k in range(-max_coeff, max_coeff + 1):
                if k:
                    S = np.eye(3, dtype=int)
                    S[i, j] = k
                    add(S)
    base = mats[:]
    for A in base:
        for B in base:
            M = A @ B
            if round(np.linalg.det(M)) == 1:
                add(M)
    return mats

def _stream_hnf(det_max, fast):
    for h11 in range(1, det_max + 1):
        for h22 in range(1, det_max // h11 + 1):
            for h33 in range(1, det_max // (h11 * h22) + 1):
                if fast:
                    yield np.diag([h11, h22, h33]).astype(int)
                    continue
                for h12 in range(h22):
                    for h13 in range(h33):
                        for h23 in range(h33):
                            yield np.array([[h11, h12, h13], [0, h22, h23], [0, 0, h33]], dtype=int)

def _search_triclinic(A_rows, max_cells, top, unimod, fast, tol, gate, anis_max, stop_on_tol):
    U_list = _small_unimodular_basis(unimod)
    best, seen = [], set()

    def score(c):
        return (c.sum_dev, c.max_dev, c.volume, c.anisotropy)

    for H in _stream_hnf(max_cells, fast):
        if np.linalg.det(H @ A_rows) <= 1e-12:
            continue
        current_gate = gate
        if current_gate is None and len(best) == top:
            current_gate = best[-1].sum_dev + 5.0

        local = None
        for U in U_list:
            P = U @ H
            if np.linalg.det(P @ A_rows) <= 1e-12:
                continue
            cand = _candidate(P, A_rows)
            if current_gate is not None and cand.sum_dev > current_gate:
                continue
            if local is None or score(cand) < score(local):
                local = cand

        if local is None or (anis_max is not None and local.anisotropy > anis_max):
            continue
        key = tuple(int(round(x * 1000)) for x in local.angles)
        if key in seen:
            continue
        seen.add(key)
        best.append(local)
        best.sort(key=score)
        del best[top:]
        if stop_on_tol and tol is not None and best[0].sum_dev <= 3.0 * tol:
            break
    return best


# ---------- public API ----------
DEFAULT_MAX_INDEX = {"plane": 12, "space": 5}


def _search_kind(A_rows, mode, angle_tol):
    """('orthogonal' | 'plane' | 'space', unique_axis) for the chosen mode."""
    system, uniq = detect_system(A_rows, angle_tol)
    if mode == "monoclinic":
        if uniq is None:
            uniq = int(np.argmax(np.abs(Cell(A_rows).angles() - 90.0)))
        return "plane", uniq
    if mode == "triclinic":
        return "space", None
    if mode != "auto":
        raise ValueError(f"Unknown mode {mode!r}")
    if system in ORTHOGONAL_SYSTEMS:
        return "orthogonal", None
    return ("plane", uniq) if system in PLANE_SYSTEMS else ("space", None)


def find_supercells(atoms, max_cells: int = 6, top: int = 2, *, mode: str = "auto",
                    method: str = "vectors", angle_tol: float = 1e-3, tol: float | None = None,
                    anis_max: float | None = None, max_index: int | None = None,
                    unimod: int = 2, fast: bool = False, gate: float | None = None,
                    stop_on_tol: bool = False) -> list[Candidate]:
    """Find integer supercells whose angles are as close to 90 degrees as possible.

    Args:
        atoms: ASE ``Atoms``, a 3x3 cell matrix, or ``[a, b, c, alpha, beta, gamma]``.
        max_cells: Largest supercell to try, in unit cells.
        top: How many candidates to return (best first).
        mode: ``"auto"`` (pick the search from the cell shape), ``"monoclinic"`` (2D
            search around the unique axis) or ``"triclinic"`` (3D search).
        method: 3D search only. ``"vectors"`` (default; fast, searches all short
            lattice vectors) or ``"hnf"`` (the original v1 search).
        angle_tol: Detection tolerance: an angle within this many degrees counts as 90.
        tol: Return the *smallest* supercells within ``tol`` degrees first. With
            ``method="hnf"`` and ``stop_on_tol``: stop once the summed deviation <= 3*tol.
        anis_max: Skip boxes whose longest/shortest edge ratio exceeds this.
        max_index: Largest integer coefficient tried for the new axes.
            Default 12 for the 2D search and 5 for the 3D search.
        unimod, fast, gate, stop_on_tol: ``method="hnf"`` controls (see the CLI docs).

    Returns:
        A list of :class:`Candidate`, best first. Cells that are already
        orthogonal (cubic, tetragonal, orthorhombic) return the identity.
    """
    A_rows = _cell_rows(atoms)
    kind, uniq = _search_kind(A_rows, mode, angle_tol)
    if kind == "orthogonal":
        return [_candidate(np.eye(3, dtype=int), A_rows)]
    if kind == "plane":
        R = max_index or DEFAULT_MAX_INDEX["plane"]
        return _search_monoclinic(A_rows, uniq, max_cells, R, top, anis_max, tol)
    if method == "hnf":
        return _search_triclinic(A_rows, max_cells, top, unimod, fast, tol, gate, anis_max, stop_on_tol)
    if method != "vectors":
        raise ValueError(f"Unknown method {method!r}")
    R = max_index or DEFAULT_MAX_INDEX["space"]
    return _search_space(A_rows, max_cells, R, top, anis_max, tol)


def best_angle_by_size(atoms, max_cells: int = 8, *, max_index: int | None = None,
                       angle_tol: float = 1e-3, mode: str = "auto") -> dict[int, float]:
    """The best achievable angle deviation using at most N unit cells.

    Returns:
        ``{N: max_dev_degrees}`` for N = 1..max_cells. Empty for cells that are
        already orthogonal.
    """
    A_rows = _cell_rows(atoms)
    kind, uniq = _search_kind(A_rows, mode, angle_tol)
    if kind == "orthogonal":
        return {}
    if kind == "plane":
        R = max_index or DEFAULT_MAX_INDEX["plane"]
        _, _, _, _, N, d, _ = _mono_pairs(A_rows, uniq, max_cells, R, None)
    else:
        _, pool = _space_pool(A_rows, max_cells, max_index or DEFAULT_MAX_INDEX["space"], None)
        N, d = pool[:, 1].astype(int), pool[:, 0]
    table, best = {}, None
    for n in range(1, max_cells + 1):
        sel = d[N == n]
        if len(sel):
            best = float(sel.min()) if best is None else min(best, float(sel.min()))
        if best is not None:
            table[n] = best
    return table


def build_supercell(atoms: Atoms, P, force_orthogonal: bool = False) -> Atoms:
    """Build the supercell ``P @ cell`` of ``atoms``.

    The result is rigidly rotated to the standard orientation (a along x, b in the
    xy plane) and atoms are wrapped into the box. No strain is applied unless
    ``force_orthogonal`` is set, which makes the box exactly 90/90/90 while keeping
    fractional coordinates; the largest Green-Lagrange strain component is stored
    in ``result.info["orthorhombify_strain"]``.
    """
    P = np.asarray(P, dtype=int)
    n_expected = len(atoms) * int(round(abs(np.linalg.det(P))))
    sc = make_supercell(atoms, P)
    if len(sc) != n_expected:
        raise RuntimeError(f"Supercell has {len(sc)} atoms, expected {n_expected}.")
    sc.set_cell(Cell.fromcellpar(sc.cell.cellpar()), scale_atoms=True)
    strain = 0.0
    if force_orthogonal:
        old = sc.cell.array.copy()
        sc.set_cell(np.diag(sc.cell.lengths()), scale_atoms=True)
        F = np.linalg.solve(old, sc.cell.array)  # new_rows = old_rows @ F
        E = 0.5 * (F @ F.T - np.eye(3))
        strain = float(np.abs(E).max())
    sc.wrap()
    sc.info["orthorhombify_P"] = P.tolist()
    sc.info["orthorhombify_strain"] = strain
    return sc


def min_distance(atoms: Atoms, cutoff: float = 3.0) -> float:
    """Shortest interatomic distance (Å), including periodic images."""
    d = neighbor_list("d", atoms, cutoff)
    return float(d.min()) if len(d) else float("inf")
