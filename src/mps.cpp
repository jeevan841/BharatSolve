// Free/fixed-format MPS reader. No external dependencies (NFR-01).
#include "bs.hpp"
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace bs {

static std::vector<std::string> tok(const std::string& line) {
  std::vector<std::string> out;
  std::istringstream ss(line);
  std::string t;
  while (ss >> t) out.push_back(t);
  return out;
}

bool readMpsStream(std::istream& in, Model& M, std::string& err) {
  std::string line;
  std::string section;
  std::unordered_map<std::string, int> rowIdx, colIdx;
  std::vector<std::string> objRow;
  std::string objName;
  // temp COO per column, built as we see COLUMNS entries in order
  std::vector<std::vector<std::pair<int, double>>> colEntries;  // per col: (row, val)
  std::vector<char> rowType;  // 'N','L','G','E'
  std::vector<double> rhs, rangeVal;
  std::vector<char> hasRange;
  bool inInteger = false;
  int lineNo = 0;

  auto getCol = [&](const std::string& name) -> int {
    auto it = colIdx.find(name);
    if (it != colIdx.end()) return it->second;
    int j = (int)colEntries.size();
    colIdx[name] = j;
    colEntries.emplace_back();
    M.cname.push_back(name);
    M.lb.push_back(0.0);
    M.ub.push_back(INF);
    M.c.push_back(0.0);
    M.isint.push_back(inInteger ? 1 : 0);
    return j;
  };
  auto getRow = [&](const std::string& name) -> int {
    auto it = rowIdx.find(name);
    if (it != rowIdx.end()) return it->second;
    int i = (int)rowType.size();
    rowIdx[name] = i;
    rowType.push_back('L');
    rhs.push_back(0.0);
    rangeVal.push_back(0.0);
    hasRange.push_back(0);
    M.rname.push_back(name);
    return i;
  };

  while (std::getline(in, line)) {
    lineNo++;
    if (line.empty()) continue;
    if (line[0] == '*') continue;
    if (!isspace((unsigned char)line[0])) {
      auto tt = tok(line);
      if (tt.empty()) continue;
      section = tt[0];
      if (section == "NAME") M.name = tt.size() > 1 ? tt[1] : "UNNAMED";
      continue;
    }
    auto f = tok(line);
    if (f.empty()) continue;
    if (section == "ROWS") {
      if (f.size() < 2) { err = "ROWS: bad line " + std::to_string(lineNo); return false; }
      char t = toupper(f[0][0]);
      if (t == 'N') {
        if (objName.empty()) objName = f[1];  // first N row is objective; others ignored
      } else {
        int i = getRow(f[1]);
        rowType[i] = t;
      }
    } else if (section == "COLUMNS") {
      if (f.size() >= 3 && f[1] == "'MARKER'") {
        if (f.size() >= 3 && (f[2] == "'INTORG'" || (f.size() > 3 && f[3] == "'INTORG'"))) inInteger = true;
        else inInteger = false;
        continue;
      }
      if (f.size() < 3) { err = "COLUMNS: bad line " + std::to_string(lineNo); return false; }
      int j = getCol(f[0]);
      for (size_t k = 1; k + 1 < f.size(); k += 2) {
        const std::string& rn = f[k];
        double v = std::stod(f[k + 1]);
        if (rn == objName) { M.c[j] += v; continue; }
        int i = getRow(rn);
        colEntries[j].emplace_back(i, v);
      }
    } else if (section == "RHS") {
      size_t start = 1;
      for (size_t k = start; k + 1 < f.size(); k += 2) {
        if (f[k] == objName) { M.objoff = -std::stod(f[k + 1]); continue; }
        int i = getRow(f[k]);
        rhs[i] = std::stod(f[k + 1]);
      }
    } else if (section == "RANGES") {
      for (size_t k = 1; k + 1 < f.size(); k += 2) {
        int i = getRow(f[k]);
        rangeVal[i] = std::stod(f[k + 1]);
        hasRange[i] = 1;
      }
    } else if (section == "BOUNDS") {
      if (f.size() < 3) continue;
      std::string bt = f[0];
      for (auto& c : bt) c = toupper(c);
      int j = getCol(f[2]);
      double v = f.size() > 3 ? std::stod(f[3]) : 0.0;
      if (bt == "UP") { M.ub[j] = v; if (v < 0 && M.lb[j] == 0) M.lb[j] = -INF; }
      else if (bt == "LO") M.lb[j] = v;
      else if (bt == "FX") { M.lb[j] = v; M.ub[j] = v; }
      else if (bt == "FR") { M.lb[j] = -INF; M.ub[j] = INF; }
      else if (bt == "MI") M.lb[j] = -INF;
      else if (bt == "PL") M.ub[j] = INF;
      else if (bt == "BV") { M.lb[j] = 0; M.ub[j] = 1; M.isint[j] = 1; }
      else if (bt == "UI") { M.ub[j] = v; M.isint[j] = 1; }
      else if (bt == "LI") { M.lb[j] = v; M.isint[j] = 1; }
    } else if (section == "OBJSENSE") {
      std::string s = f[0];
      for (auto& c : s) c = toupper(c);
      if (s == "MAX" || s == "MAXIMIZE") M.maximize = true;
    } else if (section == "SOS") {
      // SOS constraints: not yet supported (roadmap item), skip lines.
    } else if (section == "ENDATA") {
      break;
    }
  }

  if (objName.empty()) { err = "no objective (N) row found"; return false; }
  M.n = (int)colEntries.size();
  M.m = (int)rowType.size();
  // finalize row bounds from type/rhs/range
  M.rl.assign(M.m, -INF);
  M.ru.assign(M.m, INF);
  for (int i = 0; i < M.m; i++) {
    double b = rhs[i], r = rangeVal[i];
    switch (rowType[i]) {
      case 'L': M.rl[i] = hasRange[i] ? b - std::fabs(r) : -INF; M.ru[i] = b; break;
      case 'G': M.rl[i] = b; M.ru[i] = hasRange[i] ? b + std::fabs(r) : INF; break;
      case 'E': M.rl[i] = b; M.ru[i] = b;
        if (hasRange[i]) { if (r >= 0) M.ru[i] = b + r; else M.rl[i] = b + r; }
        break;
      default: M.rl[i] = -INF; M.ru[i] = INF;
    }
  }
  // CSC assembly
  M.cs.assign(M.n + 1, 0);
  for (int j = 0; j < M.n; j++) M.cs[j + 1] = M.cs[j] + (int)colEntries[j].size();
  M.ri.assign(M.cs[M.n], 0);
  M.av.assign(M.cs[M.n], 0.0);
  for (int j = 0; j < M.n; j++) {
    int p = M.cs[j];
    for (auto& e : colEntries[j]) { M.ri[p] = e.first; M.av[p] = e.second; p++; }
  }
  return true;
}

bool readMps(const std::string& path, Model& out, std::string& err) {
  std::ifstream f(path);
  if (!f) { err = "cannot open file: " + path; return false; }
  bool ok = readMpsStream(f, out, err);
  if (ok && out.name.empty()) out.name = "UNNAMED";
  return ok;
}

void Model::toMinimize() {
  if (maximize) for (auto& v : c) v = -v;
}

const char* statusName(Status s) {
  switch (s) {
    case Status::Optimal: return "optimal";
    case Status::Feasible: return "feasible";
    case Status::Infeasible: return "infeasible";
    case Status::Unbounded: return "unbounded";
    case Status::Unverified: return "unverified";
    case Status::LimitNoSol: return "limit_no_solution";
    case Status::InputError: return "input_error";
    default: return "unsupported";
  }
}
int exitCode(Status s) {
  switch (s) {
    case Status::Optimal: return 0;
    case Status::Feasible: return 1;
    case Status::Infeasible: return 2;
    case Status::Unbounded: return 3;
    case Status::Unverified: return 4;
    case Status::LimitNoSol: return 5;
    case Status::InputError: return 10;
    default: return 20;
  }
}

}  // namespace bs
