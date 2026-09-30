#!/usr/bin/env python3
"""BharatSolve UI server: wraps the bharatsolve CLI so the browser can upload
an MPS file, choose options, and see the solve log + JSON result. No solver
logic lives here -- this is a thin process wrapper around the C++ binary,
consistent with IF-01 (CLI is the source of truth; UI is optional/outside
the core, per the SRS)."""
import json, os, subprocess, tempfile, time
from flask import Flask, request, jsonify, send_from_directory

HERE = os.path.dirname(os.path.abspath(__file__))
BS = os.path.join(HERE, "..", "bin", "bharatsolve")
SAMPLES = os.path.join(HERE, "..", "data")

app = Flask(__name__, static_folder=None)

@app.route("/")
def index():
    return send_from_directory(HERE, "index.html")

@app.route("/api/samples")
def samples():
    files = sorted(f for f in os.listdir(SAMPLES) if f.endswith(".mps"))
    return jsonify(files)

@app.route("/api/sample/<name>")
def sample_text(name):
    path = os.path.join(SAMPLES, name)
    if not os.path.abspath(path).startswith(os.path.abspath(SAMPLES)) or not os.path.isfile(path):
        return jsonify({"error": "not found"}), 404
    with open(path) as f:
        return jsonify({"name": name, "text": f.read()})

@app.route("/api/solve", methods=["POST"])
def solve():
    if "file" in request.files and request.files["file"].filename:
        text = request.files["file"].read().decode("utf-8", errors="replace")
    elif "text" in request.form:
        text = request.form.get("text", "")
    else:
        body = request.get_json(silent=True) or {}
        text = body.get("text", "")
    if not text.strip():
        return jsonify({"error": "no model provided"}), 400

    time_limit = float(request.form.get("time_limit", request.args.get("time_limit", 20)))
    no_presolve = request.form.get("no_presolve", "false") == "true"
    no_scale = request.form.get("no_scale", "false") == "true"

    with tempfile.TemporaryDirectory() as td:
        mps_path = os.path.join(td, "model.mps")
        json_path = os.path.join(td, "sol.json")
        with open(mps_path, "w") as f:
            f.write(text)
        cmd = [BS, mps_path, "--json", json_path, "--time", str(time_limit), "--log", "2"]
        if no_presolve: cmd.append("--no-presolve")
        if no_scale: cmd.append("--no-scale")
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 10)
        wall = time.time() - t0
        result = {"log": proc.stderr, "wall_seconds": wall, "exit_code": proc.returncode}
        if os.path.exists(json_path):
            with open(json_path) as f:
                result["solution"] = json.load(f)
        else:
            result["error"] = proc.stderr or "solver produced no output"
        return jsonify(result)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5055))
    print(f"BharatSolve UI: http://127.0.0.1:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
