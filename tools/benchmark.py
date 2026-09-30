#!/usr/bin/env python3
"""Compare BharatSolve against SciPy's HiGHS backend on Netlib LP instances
and one MILP. Produces bench_results.csv and prints a summary table.
Usage: python3 tools/benchmark.py
"""
import glob, json, subprocess, time, sys, os
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix

BS = os.path.join(os.path.dirname(__file__), "..", "bin", "bharatsolve")
DATA = os.path.join(os.path.dirname(__file__), "..", "data")

def parse_mps_for_scipy(path):
    """Minimal MPS reader for feeding SciPy's linprog (LP only, no ranges/int)."""
    rows, row_type, rhs = [], {}, {}
    cols, col_entries = [], {}
    obj_row = None
    obj_rhs = 0.0
    bounds = {}
    section = None
    int_block = False
    name = ""
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("*"):
                continue
            indented = line[0] in " \t"
            tok = line.split()
            if not indented:
                kw = tok[0]
                if kw == "NAME": name = tok[1] if len(tok) > 1 else ""; continue
                if kw in ("ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "OBJSENSE"):
                    section = kw; continue
                if kw == "ENDATA": break
                section = None; continue
            if section == "ROWS":
                t, nm = tok[0].upper(), tok[1]
                if t == "N":
                    if obj_row is None: obj_row = nm
                else:
                    rows.append(nm); row_type[nm] = t; rhs[nm] = 0.0
            elif section == "COLUMNS":
                if "'MARKER'" in tok:
                    if "'INTORG'" in tok: int_block = True
                    if "'INTEND'" in tok: int_block = False
                    continue
                cn = tok[0]
                if cn not in col_entries:
                    col_entries[cn] = {}; cols.append(cn); bounds[cn] = [0.0, np.inf]
                for k in range(1, len(tok) - 1, 2):
                    rn, v = tok[k], float(tok[k + 1])
                    if rn == obj_row: col_entries[cn]["__obj__"] = col_entries[cn].get("__obj__", 0.0) + v
                    elif rn in row_type: col_entries[cn][rn] = col_entries[cn].get(rn, 0.0) + v
            elif section == "RHS":
                for k in range(1, len(tok) - 1, 2):
                    rn, v = tok[k], float(tok[k + 1])
                    if rn == obj_row: obj_rhs = v
                    elif rn in rhs: rhs[rn] = v
            elif section == "BOUNDS":
                bt = tok[0].upper(); cn = tok[2]
                v = float(tok[3]) if len(tok) > 3 else 0.0
                lo, hi = bounds.get(cn, [0.0, np.inf])
                if bt == "UP": hi = v
                elif bt == "LO": lo = v
                elif bt == "FX": lo = hi = v
                elif bt == "FR": lo, hi = -np.inf, np.inf
                elif bt == "MI": lo = -np.inf
                elif bt == "PL": hi = np.inf
                bounds[cn] = [lo, hi]
    n, m = len(cols), len(rows)
    c = np.zeros(n)
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    col_index = {cn: j for j, cn in enumerate(cols)}
    for j, cn in enumerate(cols):
        c[j] = col_entries[cn].get("__obj__", 0.0)
    for rn in rows:
        row = np.zeros(n)
        for cn in cols:
            v = col_entries[cn].get(rn)
            if v: row[col_index[cn]] = v
        t = row_type[rn]
        if t == "L": A_ub.append(row); b_ub.append(rhs[rn])
        elif t == "G": A_ub.append(-row); b_ub.append(-rhs[rn])
        elif t == "E": A_eq.append(row); b_eq.append(rhs[rn])
    lb = np.array([bounds[cn][0] for cn in cols])
    ub = np.array([bounds[cn][1] for cn in cols])
    bnds = list(zip(lb, ub))
    objoff = -obj_rhs  # MPS convention: true_obj = c^T x - rhs(objective row)
    return name, c, (np.array(A_ub) if A_ub else None), (np.array(b_ub) if b_ub else None), \
           (np.array(A_eq) if A_eq else None), (np.array(b_eq) if b_eq else None), bnds, objoff

def run_bharatsolve(path, time_limit=20):
    t0 = time.time()
    r = subprocess.run([BS, path, "--json", "/tmp/_bs.json", "--time", str(time_limit)],
                        capture_output=True, text=True)
    dt = time.time() - t0
    try:
        j = json.load(open("/tmp/_bs.json"))
        return j["status"], j["objective"], dt, j["verification"]["ok"]
    except Exception:
        return "error", None, dt, False

def run_scipy(path, time_limit=20, integer=False):
    try:
        name, c, A_ub, b_ub, A_eq, b_eq, bnds, objoff = parse_mps_for_scipy(path)
    except Exception as e:
        return "parse_error", None, 0.0
    t0 = time.time()
    if integer:
        from scipy.optimize import milp, LinearConstraint, Bounds
        cons = []
        if A_ub is not None: cons.append(LinearConstraint(A_ub, -np.inf, b_ub))
        if A_eq is not None: cons.append(LinearConstraint(A_eq, b_eq, b_eq))
        lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
        res = milp(c, constraints=cons, bounds=Bounds(lo, hi),
                   integrality=np.ones(len(c)), options={"time_limit": time_limit})
        dt = time.time() - t0
        obj = (res.fun + objoff) if res.success else None
        return ("optimal" if res.success else "not_optimal"), obj, dt
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bnds, method="highs",
                   options={"time_limit": time_limit})
    dt = time.time() - t0
    obj = (res.fun + objoff) if res.success else None
    return ("optimal" if res.success else "not_optimal"), obj, dt

def main():
    files = sorted(glob.glob(os.path.join(DATA, "*.mps")))
    rows = []
    print(f"{'instance':<10} {'BS status':<12} {'BS obj':>14} {'BS t(s)':>9} {'verify':>7} | "
          f"{'HiGHS status':<12} {'HiGHS obj':>14} {'HiGHS t(s)':>10} {'match':>6}")
    for path in files:
        inst = os.path.basename(path).replace(".mps", "")
        is_int = "'INTORG'" in open(path).read()
        bs_status, bs_obj, bs_t, bs_ver = run_bharatsolve(path)
        sp_status, sp_obj, sp_t = run_scipy(path, integer=is_int)
        match = ""
        if bs_obj is not None and sp_obj is not None:
            match = "yes" if abs(bs_obj - sp_obj) < 1e-3 * max(1.0, abs(sp_obj)) else "NO"
        print(f"{inst:<10} {bs_status:<12} {bs_obj if bs_obj is not None else float('nan'):>14.4f} "
              f"{bs_t:>9.4f} {str(bs_ver):>7} | {sp_status:<12} "
              f"{sp_obj if sp_obj is not None else float('nan'):>14.4f} {sp_t:>10.4f} {match:>6}")
        rows.append([inst, bs_status, bs_obj, bs_t, bs_ver, sp_status, sp_obj, sp_t, match])
    import csv
    with open(os.path.join(os.path.dirname(__file__), "..", "bench_results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["instance", "bs_status", "bs_obj", "bs_time_s", "bs_verified",
                    "highs_status", "highs_obj", "highs_time_s", "objective_match"])
        w.writerows(rows)
    print("\nwrote bench_results.csv")

if __name__ == "__main__":
    main()
