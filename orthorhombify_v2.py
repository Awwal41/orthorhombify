#!/usr/bin/env python3
"""
Orthorhombify v2: find integer supercells of a crystal that are as close to
orthorhombic (all angles 90 deg) as possible, without straining the structure.

  * Monoclinic input (exactly one angle != 90): fast, exhaustive 2D search in the
    plane perpendicular to the unique axis. The unique axis is kept as-is, so the
    two angles involving it stay exactly 90 deg and only one angle has to be fixed.
    Works for any unique-axis setting (a, b or c), e.g. P 21/c, P 21/n, P 1 1 21/b.
  * Triclinic input: the original 3D HNF x unimodular search (v1 algorithm).
  * Orthorhombic input: nothing to do, the cell is written out unchanged.

Rows convention throughout: cell rows are lattice vectors, new_rows = P @ old_rows.
"""
import argparse, glob, math, os, sys, warnings
from itertools import permutations
import numpy as np

try:
    from ase.io import read, write
    from ase.build import make_supercell
    from ase.cell import Cell
    from ase.neighborlist import neighbor_list
except Exception:
    print("This script needs ASE. Install with: pip install ase", file=sys.stderr)
    raise

# ASE warns about the 'monoclinic' cell-setting tag even when the CIF's explicit
# symmetry operators (which it uses) are present; the min-distance check below
# catches a genuinely wrong symmetry expansion.
warnings.filterwarnings("ignore", message="crystal system .* is not interpreted")

AXES = "abc"
ANGLE_NAMES = ("alpha", "beta", "gamma")


# ---------- geometry on ROWS ----------
def angle_deg_rows(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return 0.0
    c = float(np.dot(u, v)) / (nu * nv)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))

def lengths_angles_from_rows(Brows):
    a, b, c = Brows[0, :], Brows[1, :], Brows[2, :]
    la, lb, lc = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
    alpha = angle_deg_rows(b, c)
    beta  = angle_deg_rows(a, c)
    gamma = angle_deg_rows(a, b)
    return (la, lb, lc, alpha, beta, gamma)

def positive_volume_rows(Brows):
    return np.linalg.det(Brows) > 1e-12

def score_tuple(angles, lengths, vol):
    alpha, beta, gamma = angles
    dev_sum = abs(alpha - 90) + abs(beta - 90) + abs(gamma - 90)
    maxdev  = max(abs(alpha - 90), abs(beta - 90), abs(gamma - 90))
    a, b, c = lengths
    anis = (max(a, b, c) / max(1e-12, min(a, b, c)))
    return (dev_sum, maxdev, vol, anis)


# ---------- crystal-system detection ----------
def detect_system(cellpar, angle_tol):
    """Return ('orthorhombic'|'monoclinic'|'triclinic', unique_axis_index or None)."""
    off = [abs(ang - 90.0) > angle_tol for ang in cellpar[3:]]
    n_off = sum(off)
    if n_off == 0:
        return "orthorhombic", None
    if n_off == 1:
        # alpha = angle(b,c) -> unique a ; beta -> unique b ; gamma -> unique c
        return "monoclinic", off.index(True)
    return "triclinic", None


# ---------- MONOCLINIC: exact 2D search perpendicular to the unique axis ----------
def search_monoclinic(A_rows, uniq, max_cells, max_index, top_k, anis_max, tol):
    """
    In-plane lattice vectors p, q (the two non-unique axes). Enumerate primitive
    vectors v = i*p + j*q with |i|,|j| <= max_index, and every pair (v1, v2) with
    1 <= |det| <= max_cells. The new cell is (v1, v2, unique axis); the unique axis
    is perpendicular to the plane, so only angle(v1, v2) deviates from 90.
    """
    i1, i2 = [k for k in range(3) if k != uniq]
    p, q, u = A_rows[i1], A_rows[i2], A_rows[uniq]
    lu = float(np.linalg.norm(u))

    R = max_index
    coeffs = []
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            if (i, j) == (0, 0) or math.gcd(abs(i), abs(j)) != 1:
                continue
            if i < 0 or (i == 0 and j < 0):  # one representative per +/- pair
                continue
            coeffs.append((i, j))
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
    n_pairs = len(ii)
    if n_pairs == 0:
        return [], n_pairs

    N = np.abs(det[ii, jj])
    d = dev[ii, jj]
    if tol is not None:
        # smallest supercell that meets the tolerance first, then best angle
        order = np.lexsort((d, N, d > tol))
    else:
        order = np.lexsort((L[ii] + L[jj], N, np.round(d, 6)))

    vol1 = abs(float(np.linalg.det(A_rows)))
    results, seen = [], set()
    for k in order:
        a_idx, b_idx = ii[k], jj[k]
        ca, cb = C[a_idx].copy(), C[b_idx].copy()
        # shorter in-plane vector first, then make the transform right-handed
        if L[b_idx] < L[a_idx]:
            ca, cb = cb, ca
        if ca[0] * cb[1] - ca[1] * cb[0] < 0:
            cb = -cb
        P = np.zeros((3, 3), dtype=int)
        P[uniq, uniq] = 1
        P[i1, i1], P[i1, i2] = ca
        P[i2, i1], P[i2, i2] = cb
        B = P @ A_rows
        la, lb, lc, al, be, ga = lengths_angles_from_rows(B)
        key = (round(float(d[k]), 4), int(N[k]), tuple(sorted(round(x, 3) for x in (la, lb, lc))))
        if key in seen:
            continue
        seen.add(key)
        S = score_tuple((al, be, ga), (la, lb, lc), vol1 * N[k])
        results.append((S, P, (la, lb, lc), (al, be, ga), vol1 * N[k]))
        if len(results) >= top_k:
            break
    return results, n_pairs

def best_per_multiplicity(A_rows, uniq, max_cells, max_index):
    """Best achievable angle deviation for each supercell size N (trade-off table)."""
    out = {}
    for n in range(1, max_cells + 1):
        r, _ = search_monoclinic(A_rows, uniq, n, max_index, top_k=1, anis_max=None, tol=None)
        if r:
            out[n] = r[0][0][1]
    return out


# ---------- TRICLINIC: original v1 search ----------
def small_unimodular_basis(max_coeff=1):
    mats, seen = [], set()
    def add(M):
        key = tuple(M.flatten())
        if key not in seen:
            seen.add(key); mats.append(M)

    add(np.eye(3, dtype=int))
    for p in permutations(range(3)):
        P = np.zeros((3, 3), dtype=int)
        for i, j in enumerate(p): P[i, j] = 1
        if round(np.linalg.det(P)) == 1:
            add(P)
    for i in range(3):
        for j in range(3):
            if i == j: continue
            for k in range(-max_coeff, max_coeff + 1):
                if k == 0: continue
                S = np.eye(3, dtype=int); S[i, j] = k
                add(S)
    base = mats[:]
    for A in base:
        for B in base:
            C = A @ B
            if round(np.linalg.det(C)) == 1:
                add(C)
    return mats

def stream_hnf(det_max, fast):
    for h11 in range(1, det_max + 1):
        for h22 in range(1, det_max // h11 + 1):
            for h33 in range(1, det_max // (h11 * h22) + 1):
                if fast:
                    yield np.diag([h11, h22, h33]).astype(int)
                else:
                    for h12 in range(0, h22):
                        for h13 in range(0, h33):
                            for h23 in range(0, h33):
                                yield np.array([[h11, h12, h13],
                                                [  0, h22, h23],
                                                [  0,   0, h33]], dtype=int)

def search_triclinic(A_rows, max_cells, top_k, unimod_coeff, fast, tol, gate, anis_max, stop_on_tol):
    U_list = small_unimodular_basis(unimod_coeff)
    best, seen_bins = [], set()
    hnf_count = cand_checked = 0

    def consider(cand):
        S, P, lengths, angles, vol = cand
        if anis_max is not None and S[3] > anis_max:
            return
        key = tuple(int(round(x * 1000)) for x in angles)
        if key in seen_bins:
            return
        seen_bins.add(key)
        best.append(cand)
        best.sort(key=lambda t: t[0])
        if len(best) > top_k:
            best.pop()

    for H in stream_hnf(max_cells, fast):
        hnf_count += 1
        B0_rows = H @ A_rows
        if not positive_volume_rows(B0_rows):
            continue
        vol = float(abs(np.linalg.det(B0_rows)))
        current_gate = gate
        if current_gate is None and len(best) == top_k:
            current_gate = best[-1][0][0] + 5.0

        best_local = None
        for U in U_list:
            B_rows = U @ B0_rows
            if not positive_volume_rows(B_rows):
                continue
            la, lb, lc, alpha, beta, gamma = lengths_angles_from_rows(B_rows)
            S = score_tuple((alpha, beta, gamma), (la, lb, lc), vol)
            cand_checked += 1
            if current_gate is not None and S[0] > current_gate:
                continue
            if best_local is None or S < best_local[0]:
                best_local = (S, U @ H, (la, lb, lc), (alpha, beta, gamma), vol)

        if best_local is not None:
            consider(best_local)
            if stop_on_tol and tol is not None and best and best[0][0][0] <= 3.0 * tol:
                break
    return best, hnf_count, cand_checked


# ---------- building / writing ----------
def build_supercell(atoms, P, force_orthogonal):
    n_expected = len(atoms) * int(round(abs(np.linalg.det(P))))
    sc = make_supercell(atoms, P)
    if len(sc) != n_expected:
        raise RuntimeError(f"Supercell has {len(sc)} atoms, expected {n_expected}.")
    # Rigid rotation to the standard orientation (a || x, b in xy-plane);
    # fractional coordinates are kept, so this introduces no strain.
    sc.set_cell(Cell.fromcellpar(sc.cell.cellpar()), scale_atoms=True)
    strain = 0.0
    if force_orthogonal:
        old = sc.cell.array.copy()
        sc.set_cell(np.diag(sc.cell.lengths()), scale_atoms=True)
        F = np.linalg.solve(old, sc.cell.array)  # new_rows = old_rows @ F
        E = 0.5 * (F @ F.T - np.eye(3))          # Green-Lagrange strain
        strain = float(np.abs(E).max())
    sc.wrap()
    return sc, strain

def min_distance(atoms, cutoff=3.0):
    d = neighbor_list("d", atoms, cutoff)
    return float(d.min()) if len(d) else float("inf")

def fmt_matrix(P):
    return "[" + "; ".join(" ".join(f"{int(x):>3d}" for x in row) for row in P) + "]"


# ---------- driver ----------
def process(path, args):
    atoms = read(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    A_rows = atoms.cell.array
    cellpar = atoms.cell.cellpar()
    system, uniq = detect_system(cellpar, args.angle_tol)
    if args.mode != "auto":
        system = args.mode
        if system == "monoclinic":
            uniq = int(np.argmax(np.abs(cellpar[3:] - 90.0)))

    print("=" * 78)
    print(f"{path}")
    print(f"  {len(atoms)} atoms, {atoms.get_chemical_formula()}")
    print("  a,b,c = {:.4f}, {:.4f}, {:.4f} A   alpha,beta,gamma = {:.3f}, {:.3f}, {:.3f} deg"
          .format(*cellpar))
    dmin = min_distance(atoms)
    if dmin < 0.5:
        print(f"  WARNING: shortest interatomic distance is {dmin:.3f} A; the CIF symmetry "
              f"expansion may be wrong (check for duplicate/overlapping atoms).")

    if system == "orthorhombic":
        print("  Already orthorhombic within --angle-tol; writing the cell unchanged.")
        results = [(score_tuple(tuple(cellpar[3:]), tuple(cellpar[:3]), atoms.get_volume()),
                    np.eye(3, dtype=int), tuple(cellpar[:3]), tuple(cellpar[3:]), atoms.get_volume())]
    elif system == "monoclinic":
        print(f"  Detected MONOCLINIC, unique axis {AXES[uniq]} "
              f"({ANGLE_NAMES[uniq]} = {cellpar[3 + uniq]:.3f} deg is the only non-90 angle).")
        results, n_pairs = search_monoclinic(A_rows, uniq, args.max_cells, args.max_index,
                                             args.top, args.anis_max, args.tol)
        print(f"  In-plane vector pairs evaluated: {n_pairs}")
        table = best_per_multiplicity(A_rows, uniq, args.max_cells, args.max_index)
        print("  Best possible angle deviation by supercell size N (x unit cell):")
        print("    " + "  ".join(f"N={n}: {d:.3f}" for n, d in table.items()))
    else:
        print("  Detected TRICLINIC; using the 3D HNF x unimodular search.")
        results, hnf_count, cand_checked = search_triclinic(
            A_rows, args.max_cells, args.top, args.unimod, args.fast, args.tol,
            args.gate, args.anis_max, args.stop_on_tol)
        print(f"  HNF streamed: {hnf_count}  | candidates evaluated: {cand_checked}")

    if not results:
        print("  No candidates found. Increase --max-cells / --max-index, or relax --anis-max / --gate.",
              file=sys.stderr)
        return False

    os.makedirs(args.outdir, exist_ok=True)
    for rank, (S, P, lengths, angles, vol) in enumerate(results, 1):
        sc, strain = build_supercell(atoms, P, args.force_orthogonal)
        out = os.path.join(args.outdir, f"{stem}_orth_{rank}")
        write(out + ".cif", sc)
        write(out + ".POSCAR", sc, format="vasp", direct=True, sort=True)
        N = int(round(abs(np.linalg.det(P))))
        print(f"  [{rank}] N={N} ({len(sc)} atoms)  P={fmt_matrix(P)}")
        print("       a,b,c = {:.4f}, {:.4f}, {:.4f} A   angles = {:.3f}, {:.3f}, {:.3f} deg   "
              "max dev = {:.3f} deg  anis = {:.2f}".format(*lengths, *angles, S[1], S[3]))
        if args.force_orthogonal:
            print(f"       forced to exactly 90/90/90; max Green-Lagrange strain = {strain:.2e}")
        print(f"       -> {out}.cif / {out}.POSCAR")
    return True

def main():
    ap = argparse.ArgumentParser(
        description="Find near-orthorhombic integer supercells for monoclinic or triclinic crystals.")
    ap.add_argument("cif", nargs="+", help="One or more input CIF files")
    ap.add_argument("--mode", choices=["auto", "monoclinic", "triclinic"], default="auto",
                    help="Force a search type. Default: auto-detect from the cell angles.")
    ap.add_argument("--angle-tol", type=float, default=1e-3,
                    help="An angle counts as 90 deg if within this many degrees (detection). Default: 1e-3.")
    ap.add_argument("--max-cells", type=int, default=6,
                    help="Max supercell size N = det(P), i.e. number of unit cells. Default: 6.")
    ap.add_argument("--max-index", type=int, default=12,
                    help="[monoclinic] Max |integer coefficient| for the in-plane vectors. Default: 12.")
    ap.add_argument("--top", type=int, default=2, help="How many best solutions to export (1..5). Default: 2.")
    ap.add_argument("--tol", type=float, default=None,
                    help="Angle tolerance (deg). [monoclinic] prefer the SMALLEST supercell whose angle "
                         "deviation is <= tol. [triclinic] used with --stop-on-tol.")
    ap.add_argument("--anis-max", type=float, default=None,
                    help="Reject boxes with max/min edge length above this value, e.g. 3.0.")
    ap.add_argument("--force-orthogonal", action="store_true",
                    help="Set the final cell to exactly 90/90/90 (keeps fractional coords; introduces a "
                         "small strain, which is reported). Useful for MD codes needing orthogonal boxes.")
    ap.add_argument("--outdir", default=".", help="Output directory. Default: current directory.")
    # triclinic-only options (same as v1)
    ap.add_argument("--unimod", type=int, default=2, help="[triclinic] Unimodular shear coeff limit. Default: 2.")
    ap.add_argument("--fast", action="store_true", help="[triclinic] Only diagonal H.")
    ap.add_argument("--gate", type=float, default=None, help="[triclinic] Angle-sum gate (deg), e.g. 12.")
    ap.add_argument("--stop-on-tol", action="store_true", help="[triclinic] Stop once sum dev <= 3*tol.")
    args = ap.parse_args()

    if not 1 <= args.top <= 5:
        print("--top must be between 1 and 5.", file=sys.stderr)
        sys.exit(2)

    # Windows shells do not expand wildcards, so do it here
    paths = []
    for pattern in args.cif:
        paths.extend(sorted(glob.glob(pattern)) or [pattern])
    ok = [process(path, args) for path in paths]
    sys.exit(0 if all(ok) else 2)

if __name__ == "__main__":
    main()
