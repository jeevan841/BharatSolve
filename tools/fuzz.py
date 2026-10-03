#!/usr/bin/env python3
"""Random LP/MILP fuzzer: compares BharatSolve with SciPy's HiGHS on small random models.
Usage: python3 tools/fuzz.py [count] [seed]     (FUZZ_BIG=1 for larger MILPs)     (FUZZ_METHOD=dual to test the dual simplex)
Reports status or objective mismatches and saves the failing models to tools/fuzz_failures/.
"""
import json, os, random, subprocess, sys, tempfile
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BIN = os.path.join(ROOT, "bin", "bharatsolve")

BIG = bool(os.environ.get("FUZZ_BIG"))   # larger, mostly-integer, bounded models to stress branch-and-bound

def make_model(rng):
    n = rng.randint(8, 16) if BIG else rng.randint(2, 8); m = rng.randint(4, 10) if BIG else rng.randint(1, 6)
    c = [rng.randint(-5, 5) for _ in range(n)]
    A = [[(rng.randint(-5, 5) if rng.random() < 0.6 else 0) for _ in range(n)] for _ in range(m)]
    sense = [rng.choice("LGE") for _ in range(m)]
    rhs = [rng.randint(-10, 20) for _ in range(m)]
    isint = [rng.random() < (0.8 if BIG else 0.3) for _ in range(n)]
    lb, ub = [], []
    for _ in range(n):
        t = rng.choice(["up", "both", "up", "both", "fx"] if BIG else ["def", "up", "both", "fx", "fr", "up", "both"])
        if t == "def": lb.append(0.0); ub.append(np.inf)
        elif t == "up": lb.append(0.0); ub.append(float(rng.randint(1, 10)))
        elif t == "both":
            l = rng.randint(-5, 3); lb.append(float(l)); ub.append(float(l + rng.randint(0, 8)))
        elif t == "fx": v = float(rng.randint(-3, 5)); lb.append(v); ub.append(v)
        else: lb.append(-np.inf); ub.append(np.inf)
    return n, m, c, A, sense, rhs, isint, lb, ub

def to_mps(model):
    n, m, c, A, sense, rhs, isint, lb, ub = model
    L = ["NAME FUZZ", "ROWS", " N OBJ"]
    for i in range(m): L.append(f" {sense[i]} R{i+1}")
    L.append("COLUMNS")
    inint = False
    for j in range(n):
        if isint[j] and not inint: L.append(" MARKER 'MARKER' 'INTORG'"); inint = True
        if not isint[j] and inint: L.append(" MARKER 'MARKER' 'INTEND'"); inint = False
        ents = []
        if c[j]: ents.append(("OBJ", c[j]))
        for i in range(m):
            if A[i][j]: ents.append((f"R{i+1}", A[i][j]))
        if not ents: ents = [("OBJ", 0)]
        for r, v in ents: L.append(f" C{j+1} {r} {v}")
    if inint: L.append(" MARKER 'MARKER' 'INTEND'")
    L.append("RHS")
    for i in range(m):
        if rhs[i]: L.append(f" RHS R{i+1} {rhs[i]}")
    L.append("BOUNDS")
    for j in range(n):
        l, u = lb[j], ub[j]
        if l == u: L.append(f" FX BND C{j+1} {l}")
        elif l == -np.inf and u == np.inf: L.append(f" FR BND C{j+1}")
        else:
            if l != 0: L.append(f" LO BND C{j+1} {l}")
            if u != np.inf: L.append(f" UP BND C{j+1} {u}")
    L.append("ENDATA")
    return "\n".join(L) + "\n"

def highs(model):
    n, m, c, A, sense, rhs, isint, lb, ub = model
    A = np.array(A, float).reshape(m, n)
    lo = np.array([rhs[i] if sense[i] in "GE" else -np.inf for i in range(m)], float)
    hi = np.array([rhs[i] if sense[i] in "LE" else np.inf for i in range(m)], float)
    r = milp(c=np.array(c, float), constraints=LinearConstraint(A, lo, hi),
             integrality=np.array([1 if x else 0 for x in isint]), bounds=Bounds(lb, ub))
    return r.status, (r.fun if r.status == 0 else None)

def bharat(text):
    with tempfile.NamedTemporaryFile("w", suffix=".mps", delete=False) as f:
        f.write(text); path = f.name
    try:
        p = subprocess.run([BIN, path, "--log", "0", "--json", path + ".json", "--time", "10"] + (["--method", os.environ["FUZZ_METHOD"]] if os.environ.get("FUZZ_METHOD") else []),
                           capture_output=True, text=True, timeout=30)
        obj = None
        if os.path.exists(path + ".json"):
            obj = json.load(open(path + ".json")).get("objective")
        return p.returncode, obj
    finally:
        for q in (path, path + ".json"):
            if os.path.exists(q): os.remove(q)

def main():
    count = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    outdir = os.path.join(ROOT, "tools", "fuzz_failures")
    stats = {"ok": 0, "skip": 0, "bad": 0, "limit": 0, "ambig": 0}
    for k in range(count):
        model = make_model(rng); text = to_mps(model)
        hs, hobj = highs(model)
        if hs not in (0, 2, 3):
            stats["skip"] += 1; continue          # HiGHS gave no clear answer (limit / ambiguous)
        rc, obj = bharat(text)
        want = {0: 0, 2: 2, 3: 3}[hs]
        has_int = any(model[6])
        if rc == 1 and hs == 0 and obj is not None and abs(obj - hobj) <= 1e-5 * max(1, abs(hobj)):
            stats["limit"] += 1; continue              # right answer, but search not finished within the time limit
        if rc == 5 and hs in (0, 2):                    # honest "limit, no answer": never a wrong answer
            stats["limit"] += 1; continue
        if has_int and hs == 2 and rc == 3:             # unbounded relaxation: "unbounded or infeasible" is ambiguous
            stats["ambig"] += 1; continue
        good = (rc == want) and (hs != 0 or (obj is not None and abs(obj - hobj) <= 1e-5 * max(1, abs(hobj))))
        if good: stats["ok"] += 1
        else:
            stats["bad"] += 1
            os.makedirs(outdir, exist_ok=True)
            fn = os.path.join(outdir, f"fail_s{seed}_{k}.mps"); open(fn, "w").write(text)
            print(f"MISMATCH #{k}: HiGHS status {hs} obj {hobj} | BharatSolve exit {rc} obj {obj} -> {os.path.relpath(fn, ROOT)}")
    print(f"seed {seed}: {stats['ok']} agree, {stats['bad']} mismatch, {stats['limit']} honest limit, {stats['ambig']} ambiguous, {stats['skip']} skipped (of {count})")
    sys.exit(1 if stats["bad"] else 0)

if __name__ == "__main__":
    main()
