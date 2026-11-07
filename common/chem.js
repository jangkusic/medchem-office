/* 의약화학 오피스 — 물성 계산 엔진 (순수 JavaScript, 인터넷 불필요)
 * - SMILES / MOL(V2000) 읽기, 암묵 수소, 고리, 방향족 판정
 * - 분자식, 분자량, 정확 질량
 * - Wildman–Crippen logP·MR (RDKit Data/Crippen.txt 표), TPSA (Ertl, RDKit 규칙)
 * - HBD/HBA (RDKit 정의 + Lipinski NHOH/NO), 회전 결합 수 (RDKit strict), 고리·방향족 고리 수, Fsp3
 * - Lipinski / Veber 판정, Morgan(ECFP4 유사) 지문과 Tanimoto 유사도, 부분 구조 검색
 */
(function (root) {
  "use strict";

  // ---------- 원소 ----------
  var EL = {
    H: [1, 1.008, 1.00782503], He: [2, 4.0026, 4.00260325], Li: [3, 6.94, 7.0160045], Be: [4, 9.0122, 9.0121822],
    B: [5, 10.81, 11.0093054], C: [6, 12.011, 12], N: [7, 14.007, 14.0030740], O: [8, 15.999, 15.9949146],
    F: [9, 18.998, 18.9984032], Ne: [10, 20.18, 19.9924402], Na: [11, 22.99, 22.9897693], Mg: [12, 24.305, 23.9850419],
    Al: [13, 26.982, 26.9815384], Si: [14, 28.085, 27.9769265], P: [15, 30.974, 30.9737615], S: [16, 32.067, 31.9720707],
    Cl: [17, 35.453, 34.9688527], Ar: [18, 39.948, 39.9623831], K: [19, 39.098, 38.9637069], Ca: [20, 40.078, 39.9625912],
    Ti: [22, 47.867, 47.9479463], Cr: [24, 51.996, 51.9405075], Mn: [25, 54.938, 54.9380451], Fe: [26, 55.845, 55.9349375],
    Co: [27, 58.933, 58.933195], Ni: [28, 58.693, 57.9353429], Cu: [29, 63.546, 62.9295975], Zn: [30, 65.38, 63.9291422],
    Ga: [31, 69.723, 68.9255736], Ge: [32, 72.63, 73.9211778], As: [33, 74.922, 74.9215965], Se: [34, 78.971, 79.9165213],
    Br: [35, 79.904, 78.9183371], Kr: [36, 83.798, 83.911507], Rb: [37, 85.468, 84.9117897], Sr: [38, 87.62, 87.9056121],
    Pd: [46, 106.42, 105.903486], Ag: [47, 107.868, 106.905097], Cd: [48, 112.414, 113.9033585], Sn: [50, 118.71, 119.9021947],
    Sb: [51, 121.76, 120.9038157], Te: [52, 127.6, 129.9062244], I: [53, 126.904, 126.904473], Xe: [54, 131.29, 131.9041535],
    Cs: [55, 132.905, 132.905452], Ba: [56, 137.327, 137.9052472], Pt: [78, 195.084, 194.9647911], Au: [79, 196.967, 196.9665687],
    Hg: [80, 200.59, 201.970643], Pb: [82, 207.2, 207.9766521], Bi: [83, 208.98, 208.9803987], Sc: [21, 44.956, 44.9559119],
    V: [23, 50.942, 50.9439595], Y: [39, 88.906, 88.9058483], Zr: [40, 91.224, 89.9047044], Nb: [41, 92.906, 92.9063781],
    Mo: [42, 95.95, 97.9054082], Ru: [44, 101.07, 101.9043493], Rh: [45, 102.906, 102.905504], In: [49, 114.818, 114.903878],
    Hf: [72, 178.49, 179.94655], Ta: [73, 180.948, 180.9479958], W: [74, 183.84, 183.9509312], Re: [75, 186.207, 186.9557531],
    Os: [76, 190.23, 191.9614807], Ir: [77, 192.217, 192.9629264], Tl: [81, 204.38, 204.9744275], Po: [84, 209, 208.9824304]
  };
  var BYNUM = {}; Object.keys(EL).forEach(function (s) { BYNUM[EL[s][0]] = s; });
  var VALENCE = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1, 3, 5], Si: [4], Se: [2, 4, 6], As: [3, 5], H: [1] };

  // ---------- 분자 ----------
  function Mol() { this.atoms = []; this.bonds = []; }
  Mol.prototype.addAtom = function (a) {
    a.idx = this.atoms.length; a.charge = a.charge || 0; a.arom = !!a.arom; a.nbrs = [];
    if (a.hExplicit == null) a.hExplicit = null; this.atoms.push(a); return a;
  };
  Mol.prototype.addBond = function (i, j, order, arom) {
    var b = { idx: this.bonds.length, a: i, b: j, order: order, arom: !!arom };
    this.bonds.push(b); this.atoms[i].nbrs.push({ atom: j, bond: b }); this.atoms[j].nbrs.push({ atom: i, bond: b }); return b;
  };
  function other(b, i) { return b.a === i ? b.b : b.a; }

  // ---------- SMILES 읽기 ----------
  var ORGANIC = ["Cl", "Br", "B", "C", "N", "O", "P", "S", "F", "I", "b", "c", "n", "o", "p", "s"];
  function parseSmiles(smi) {
    smi = String(smi).trim().split(/\s+/)[0];
    var m = new Mol(), i = 0, prev = null, stack = [], ring = {}, bondSym = null;
    function readBond() { var c = smi[i]; if ("-=#:$/\\.".indexOf(c) >= 0) { i++; return c; } return null; }
    function link(a) {
      if (prev != null && bondSym !== ".") {
        var pa = m.atoms[prev], na = m.atoms[a], o = 1, ar = false;
        if (bondSym === "=") o = 2; else if (bondSym === "#") o = 3; else if (bondSym === "$") o = 4;
        else if (bondSym === ":") { ar = true; o = 1.5; }
        else if (!bondSym && pa.arom && na.arom) { ar = true; o = 1.5; }
        m.addBond(prev, a, o, ar);
      }
      prev = a; bondSym = null;
    }
    while (i < smi.length) {
      var c = smi[i];
      if (c === "(") { stack.push(prev); i++; continue; }
      if (c === ")") { prev = stack.pop(); i++; continue; }
      var bs = readBond(); if (bs) { if (bs === "/" || bs === "\\") bs = null; bondSym = bs === null ? bondSym : bs; if (bs === ".") { prev = null; bondSym = null; } continue; }
      if (c === "%" || /\d/.test(c)) {
        var num; if (c === "%") { num = smi.substr(i + 1, 2); i += 3; } else { num = c; i++; }
        if (ring[num]) {
          var r = ring[num], pa = m.atoms[r.atom], na = m.atoms[prev], sym = bondSym || r.bond, o = 1, ar = false;
          if (sym === "=") o = 2; else if (sym === "#") o = 3; else if (sym === ":") { ar = true; o = 1.5; }
          else if (!sym && pa.arom && na.arom) { ar = true; o = 1.5; }
          m.addBond(r.atom, prev, o, ar); delete ring[num];
        } else ring[num] = { atom: prev, bond: bondSym };
        bondSym = null; continue;
      }
      if (c === "[") {
        var end = smi.indexOf("]", i), body = smi.slice(i + 1, end); i = end + 1;
        var mm = body.match(/^(\d+)?([A-Z][a-z]?|[a-z][a-z]?|\*)(@*)(?:TH\d|AL\d|SP\d|TB\d+|OH\d+)?(H\d*)?([+-]+\d*)?(:\d+)?$/);
        if (!mm) throw new Error("SMILES를 읽을 수 없어요: [" + body + "]");
        var sy = mm[2], arom = false;
        if (/^[a-z]/.test(sy)) { arom = true; sy = sy[0].toUpperCase() + sy.slice(1); }
        if (!EL[sy] && sy !== "*") throw new Error("모르는 원소: " + sy);
        var h = mm[4] ? (mm[4].length > 1 ? +mm[4].slice(1) : 1) : 0, ch = 0;
        if (mm[5]) { var sg = mm[5][0] === "+" ? 1 : -1, rest = mm[5].replace(/[+-]/g, ""); ch = rest ? sg * +rest : sg * mm[5].length; }
        var at = m.addAtom({ el: sy === "*" ? "C" : sy, arom: arom, charge: ch, hExplicit: h, iso: mm[1] ? +mm[1] : 0 });
        link(at.idx); continue;
      }
      var got = null;
      for (var k = 0; k < ORGANIC.length; k++) if (smi.substr(i, ORGANIC[k].length) === ORGANIC[k]) { got = ORGANIC[k]; break; }
      if (!got) { if (c === "*") got = "C"; else throw new Error("SMILES를 읽을 수 없어요: '" + c + "' (" + (i + 1) + "번째 글자)"); }
      i += got.length;
      var ar2 = /^[a-z]/.test(got), el = ar2 ? got.toUpperCase() : got;
      var a2 = m.addAtom({ el: el, arom: ar2, charge: 0, hExplicit: null });
      link(a2.idx);
    }
    if (Object.keys(ring).length) throw new Error("닫히지 않은 고리 번호가 있어요");
    return finish(m, true);
  }

  // ---------- MOL(V2000) 읽기 ----------
  function parseMolfile(txt) {
    var L = String(txt).replace(/\r/g, "").split("\n"), m = new Mol();
    var ci = -1; for (var k = 0; k < Math.min(L.length, 10); k++) if (/V2000/.test(L[k])) { ci = k; break; }
    if (ci < 0) { if (/V3000/.test(txt)) return parseMolV3000(L); throw new Error("MOL 파일 형식을 읽을 수 없어요"); }
    var na = +L[ci].substr(0, 3), nb = +L[ci].substr(3, 3);
    for (var a = 0; a < na; a++) {
      var ln = L[ci + 1 + a], sym = ln.substr(31, 3).trim(), cc = +ln.substr(36, 3) || 0;
      var chg = cc ? 4 - cc : 0; if (cc === 4) chg = 0;
      if (sym === "D" || sym === "T") sym = "H";
      if (!EL[sym]) sym = /^R|^\*|^A$|^Q$|^L$/.test(sym) ? "C" : sym;
      m.addAtom({ el: sym, charge: chg, hExplicit: null, x: +ln.substr(0, 10), y: +ln.substr(10, 10) });
    }
    for (var b = 0; b < nb; b++) {
      var lb = L[ci + 1 + na + b], i1 = +lb.substr(0, 3) - 1, i2 = +lb.substr(3, 3) - 1, t = +lb.substr(6, 3);
      var st = +(lb.substr(9, 3) || 0) || 0;               // 쐐기(1)·점선(6) — 입체배치에 필요해요
      var nb2 = t === 4 ? m.addBond(i1, i2, 1.5, true) : m.addBond(i1, i2, t >= 1 && t <= 3 ? t : 1, false);
      if (nb2) { nb2.stereo = st; nb2.from = i1; }         // from: 쐐기의 좁은 쪽
    }
    for (var p = ci + 1 + na + nb; p < L.length; p++) {
      var s = L[p]; if (/^M  END/.test(s)) break;
      if (/^M  (CHG|ISO|RAD)/.test(s)) {
        var n = +s.substr(6, 3);
        for (var q = 0; q < n; q++) {
          var ai = +s.substr(10 + 8 * q, 3) - 1, v = +s.substr(14 + 8 * q, 3);
          if (/CHG/.test(s)) m.atoms[ai].charge = v;
        }
      }
    }
    m.atoms.forEach(function (at) { if (m.bonds.some(function (bb) { return bb.arom && (bb.a === at.idx || bb.b === at.idx); })) at.arom = true; });
    return finish(m, false);
  }
  function parseMolV3000(L) {
    var m = new Mol(), mode = null, map = {};
    L.forEach(function (s) {
      if (/BEGIN ATOM/.test(s)) { mode = "a"; return; } if (/BEGIN BOND/.test(s)) { mode = "b"; return; }
      if (/END (ATOM|BOND)/.test(s)) { mode = null; return; }
      var t = s.replace(/^M  V30 /, "").trim().split(/\s+/);
      if (mode === "a" && t.length > 4) {
        var chg = 0; t.forEach(function (x) { var mm = x.match(/^CHG=(-?\d+)/); if (mm) chg = +mm[1]; });
        var sym = EL[t[1]] ? t[1] : "C"; map[t[0]] = m.addAtom({ el: sym, charge: chg, hExplicit: null, x: +t[2], y: +t[3] }).idx;
      } else if (mode === "b" && t.length > 3) {
        var ty = +t[1]; if (ty === 4) m.addBond(map[t[2]], map[t[3]], 1.5, true); else m.addBond(map[t[2]], map[t[3]], ty, false);
      }
    });
    m.atoms.forEach(function (at) { if (at.nbrs.some(function (n) { return n.bond.arom; })) at.arom = true; });
    return finish(m, false);
  }

  // ---------- 수소·고리·방향족 ----------
  function bondValence(at) { return at.nbrs.reduce(function (s, n) { return s + (n.bond.arom ? 1 : n.bond.order); }, 0); }
  function implicitH(at) {
    if (at.hExplicit != null) return at.hExplicit;
    var vals = VALENCE[at.el]; if (!vals) return 0;
    var v = bondValence(at), adj = at.charge;
    if (at.arom && !at.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2; })) v += 1;
    // 전하: N+, O+, S+ 는 원자가 +1 / C+, B- 등은 -1 / 음전하 N-, O- 는 -1
    var base;
    if (at.el === "C") base = 4 - Math.abs(adj);
    else if (at.el === "B") base = 3 + adj;
    else if (["N", "O", "S", "P", "Se", "As"].indexOf(at.el) >= 0) base = null; else base = null;
    if (base != null) return Math.max(0, base - v);
    if (at.arom) return Math.max(0, vals[0] + adj - v);
    for (var k = 0; k < vals.length; k++) { var tv = vals[k] + adj; if (tv >= v) return Math.max(0, tv - v); }
    return 0;
  }
  function ringsOf(m) { // 최소 고리 집합(SSSR 근사): 각 결합마다 그것을 지나는 가장 작은 고리
    var n = m.atoms.length, rings = [], seen = {};
    m.bonds.forEach(function (b) {
      // b를 빼고 a→b 최단 경로
      var prev = new Array(n).fill(-1), q = [b.a], vis = new Array(n).fill(false); vis[b.a] = true;
      while (q.length) {
        var u = q.shift(); if (u === b.b) break;
        m.atoms[u].nbrs.forEach(function (nb) { if (nb.bond === b || vis[nb.atom]) return; vis[nb.atom] = true; prev[nb.atom] = u; q.push(nb.atom); });
      }
      if (!vis[b.b]) return;
      var path = [], x = b.b; while (x !== -1) { path.push(x); x = prev[x]; }
      var key = path.slice().sort(function (p, q2) { return p - q2; }).join(",");
      if (!seen[key]) { seen[key] = 1; rings.push(path); }
    });
    // 큰 고리 중 작은 고리들의 합으로 설명되는 것 제거 (순환 기저 크기로 자르기)
    rings.sort(function (p, q2) { return p.length - q2.length; });
    var need = m.bonds.length - n + components(m), basis = [], vecs = [];
    rings.forEach(function (r) {
      if (basis.length >= need) return;
      var v = ringBondVec(m, r);
      if (independent(vecs, v)) { vecs.push(v); basis.push(r); }
    });
    return { sssr: basis, all: rings };
  }
  function components(m) {
    var c = 0, vis = new Array(m.atoms.length).fill(false);
    for (var i = 0; i < m.atoms.length; i++) { if (vis[i]) continue; c++; var st = [i]; vis[i] = true; while (st.length) { var u = st.pop(); m.atoms[u].nbrs.forEach(function (nb) { if (!vis[nb.atom]) { vis[nb.atom] = true; st.push(nb.atom); } }); } }
    return c;
  }
  function ringBondVec(m, r) {
    var v = new Uint8Array(m.bonds.length);
    for (var i = 0; i < r.length; i++) { var a = r[i], b = r[(i + 1) % r.length]; m.atoms[a].nbrs.forEach(function (nb) { if (nb.atom === b) v[nb.bond.idx] = 1; }); }
    return v;
  }
  function independent(vecs, v) { // GF(2) 소거
    var w = Uint8Array.from(v), rows = vecs.map(function (x) { return Uint8Array.from(x); });
    var M = rows.concat([w]), r = 0, cols = v.length;
    for (var c = 0; c < cols && r < M.length; c++) {
      var piv = -1; for (var i = r; i < M.length; i++) if (M[i][c]) { piv = i; break; }
      if (piv < 0) continue; var t = M[r]; M[r] = M[piv]; M[piv] = t;
      for (var j = 0; j < M.length; j++) if (j !== r && M[j][c]) for (var k = 0; k < cols; k++) M[j][k] ^= M[r][k];
      r++;
    }
    return r === M.length;
  }
  function perceiveAromatic(m, rings) {
    // 원자별 π 전자 기여 (RDKit 방식 근사)
    function piCount(at) {
      var ringDouble = false, exoDouble = null;
      at.nbrs.forEach(function (nb) {
        if (nb.bond.order === 2 && !nb.bond.arom) { if (nb.bond.inRing) ringDouble = true; else exoDouble = m.atoms[nb.atom]; }
      });
      var deg = at.nbrs.length + at.h, el = at.el, ch = at.charge;
      if (at.arom) return null; // 이미 방향족 표시
      if (el === "C") {
        if (ringDouble) return 1;
        if (exoDouble) return ["O", "N", "S"].indexOf(exoDouble.el) >= 0 ? 0 : -1;
        if (ch === -1 && deg === 3) return 2; if (ch === 1 && deg === 3) return 0; return -1;
      }
      if (el === "N" || el === "P") {
        if (ringDouble) return 1;
        if (exoDouble && el === "N" && ch === 1) return 1;
        if (deg === 3 && ch === 0) return 2; if (deg === 2 && ch === -1) return 2; return -1;
      }
      if (el === "O" || el === "S" || el === "Se") {
        if (ringDouble && ch === 1) return 1;
        if (ringDouble && el === "S" && ch === 0 && deg <= 3) return -1;
        if (!ringDouble && !exoDouble && at.nbrs.length === 2 && ch === 0) return 2; return -1;
      }
      if (el === "B") { if (!ringDouble && deg === 3) return 0; return -1; }
      return -1;
    }
    var changed = true, guard = 0;
    while (changed && guard++ < 6) {
      changed = false;
      var cand = rings.all.filter(function (r) { return r.length >= 4 && r.length <= 8; });
      // 단일 고리 + 두 고리 합친 고리계(나프탈렌 같은 경우는 단일로도 충분)
      cand.forEach(function (r) {
        if (r.every(function (i) { return m.atoms[i].arom; })) return;
        var sum = 0, ok = true;
        r.forEach(function (i) { var at = m.atoms[i]; var p = at.arom ? 1 : piCount(at); if (p == null) p = 1; if (p < 0) ok = false; sum += p; });
        if (!ok) return;
        if ((sum - 2) % 4 === 0 && sum >= 2) {
          r.forEach(function (i) { m.atoms[i].arom = true; });
          for (var k = 0; k < r.length; k++) { var a = r[k], b = r[(k + 1) % r.length]; m.atoms[a].nbrs.forEach(function (nb) { if (nb.atom === b) { nb.bond.arom = true; } }); }
          changed = true;
        }
      });
    }
  }
  function finish(m, fromSmiles) {
    var R = ringsOf(m);
    m.rings = R.sssr; m.allRings = R.all;
    m.bonds.forEach(function (b) { b.inRing = false; });
    m.atoms.forEach(function (a) { a.inRing = false; a.ringSizes = []; });
    R.all.forEach(function (r) {
      r.forEach(function (i) { m.atoms[i].inRing = true; if (m.atoms[i].ringSizes.indexOf(r.length) < 0) m.atoms[i].ringSizes.push(r.length); });
      for (var k = 0; k < r.length; k++) { var a = r[k], b = r[(k + 1) % r.length]; m.atoms[a].nbrs.forEach(function (nb) { if (nb.atom === b) nb.bond.inRing = true; }); }
    });
    // 고리 밖 결합은 방향족이 아님 (예: 비아릴 결합)
    m.bonds.forEach(function (b) { if (b.arom && !b.inRing) { b.arom = false; b.order = 1; } });
    // 수소 수는 입력 그대로(케쿨레 구조 기준) 고정한 뒤 방향족 판정
    m.atoms.forEach(function (a) { a.h = implicitH(a); a.hExplicit = a.h; });
    perceiveAromatic(m, R);
    m.bonds.forEach(function (b) { if (b.arom && !(m.atoms[b.a].arom && m.atoms[b.b].arom)) b.arom = false; });
    return m;
  }

  // ---------- 분자식·질량 ----------
  function formulaCounts(m) {
    var c = {};
    m.atoms.forEach(function (a) { c[a.el] = (c[a.el] || 0) + 1; if (a.h) c.H = (c.H || 0) + a.h; });
    return c;
  }
  function hill(c) {
    var keys = Object.keys(c), out = "";
    function f(k) { return k + (c[k] > 1 ? c[k] : ""); }
    if (c.C) { out += f("C"); if (c.H) out += f("H"); keys = keys.filter(function (k) { return k !== "C" && k !== "H"; }); }
    keys.sort().forEach(function (k) { out += f(k); });
    return out;
  }
  function netCharge(m) { return m.atoms.reduce(function (s, a) { return s + a.charge; }, 0); }

  // ---------- 작은 SMARTS (Crippen 표용: 트리, 재귀/고리닫힘 없음) ----------
  function parseSmarts(s) {
    var i = 0, atoms = [], edges = [], stack = [], prev = -1, bond = null;
    function atomExpr(txt, bracket) { return { txt: txt, test: compileAtom(txt, bracket) }; }
    while (i < s.length) {
      var c = s[i];
      if (c === "(") { stack.push(prev); i++; continue; }
      if (c === ")") { prev = stack.pop(); i++; continue; }
      if ("-=#:~".indexOf(c) >= 0) { bond = c; i++; continue; }
      var a;
      if (c === "[") { var e = s.indexOf("]", i); a = atomExpr(s.slice(i + 1, e), true); i = e + 1; }
      else {
        var t = s.substr(i, 2) === "Cl" || s.substr(i, 2) === "Br" ? s.substr(i, 2) : c;
        a = atomExpr(t, false); i += t.length;
      }
      atoms.push(a);
      if (prev >= 0) edges.push({ from: prev, to: atoms.length - 1, bond: bond });
      prev = atoms.length - 1; bond = null;
    }
    var kids = atoms.map(function () { return []; });
    edges.forEach(function (e) { kids[e.from].push(e); });
    return { atoms: atoms, kids: kids };
  }
  function compileAtom(txt, bracket) {
    if (!bracket) return prim(txt);
    // ';' 낮은 AND, ',' OR, '&' 또는 붙여쓰기 높은 AND, '!' NOT
    var lows = txt.split(";").map(function (part) {
      var ors = part.split(",").map(function (o) { return hiAnd(o); });
      return function (a, g) { return ors.some(function (f) { return f(a, g); }); };
    });
    return function (a, g) { return lows.every(function (f) { return f(a, g); }); };
  }
  function hiAnd(s) {
    var parts = [], i = 0;
    while (i < s.length) {
      var neg = false; if (s[i] === "&") { i++; continue; }
      while (s[i] === "!") { neg = !neg; i++; }
      var m = s.slice(i).match(/^(#\d+|H\d*|X\d+|D\d+|v\d+|R\d*|r\d+|[+-]\d+|[+]+|[-]+|Cl|Br|[A-Z][a-z]?(?![a-z])|[A-Z]|[cnosp]|[aA*])/);
      if (!m) throw new Error("SMARTS: " + s);
      var tok = m[1]; i += tok.length;
      // 'H' 다음 숫자, 'C' 'H' 구분: 원소 기호가 H로 시작하는 경우(Hg 등) 제외
      var f = prim(tok);
      parts.push(neg ? (function (ff) { return function (a, g) { return !ff(a, g); }; })(f) : f);
    }
    return function (a, g) { return parts.every(function (f) { return f(a, g); }); };
  }
  var ALIPH = { C: 1, N: 1, O: 1, S: 1, P: 1, F: 1, Cl: 1, Br: 1, I: 1, B: 1 };
  function prim(t) {
    if (t === "*") return function () { return true; };
    if (t === "a") return function (a) { return a.arom; };
    if (t === "A") return function (a) { return !a.arom; };
    if (/^#\d+$/.test(t)) { var z = +t.slice(1); return function (a) { return a.z === z; }; }
    if (/^H\d*$/.test(t)) { var hn = t.length > 1 ? +t.slice(1) : 1; return function (a) { return a.th === hn; }; }
    if (/^X\d+$/.test(t)) { var x = +t.slice(1); return function (a) { return a.tx === x; }; }
    if (/^D\d+$/.test(t)) { var d = +t.slice(1); return function (a) { return a.deg === d; }; }
    if (/^v\d+$/.test(t)) { var v = +t.slice(1); return function (a) { return a.val === v; }; }
    if (t === "R") return function (a) { return a.inRing; };
    if (/^[+-]/.test(t)) {
      var ch; if (/^[+-]\d+$/.test(t)) ch = (t[0] === "+" ? 1 : -1) * +t.slice(1); else ch = (t[0] === "+" ? 1 : -1) * t.length;
      return function (a) { return a.charge === ch; };
    }
    if (/^[cnosp]$/.test(t)) { var E = t.toUpperCase(); return function (a) { return a.arom && a.el === E; }; }
    if (EL[t]) return function (a) { return a.el === t && !a.arom; };
    throw new Error("SMARTS 원자: " + t);
  }
  function bondOk(sym, b) {
    if (sym === "~") return true;
    if (sym === ":") return b.arom;
    if (sym === "=") return !b.arom && b.order === 2;
    if (sym === "#") return !b.arom && b.order === 3;
    if (sym === "-") return !b.arom && b.order === 1;
    return b.arom || b.order === 1; // 기본: 단일 또는 방향족
  }
  function matchAt(q, g, start) {
    var used = {}; used[start] = 1;
    function rec(qi, ti) {
      var kids = q.kids[qi];
      function tryK(k) {
        if (k === kids.length) return true;
        var e = kids[k], at = g.atoms[ti];
        for (var n = 0; n < at.nbrs.length; n++) {
          var nb = at.nbrs[n]; if (used[nb.atom]) continue;
          if (!bondOk(e.bond, nb.bond)) continue;
          if (!q.atoms[e.to].test(g.atoms[nb.atom], g)) continue;
          used[nb.atom] = 1;
          if (rec(e.to, nb.atom) && tryK(k + 1)) return true;
          delete used[nb.atom];
        }
        return false;
      }
      return tryK(0);
    }
    return q.atoms[0].test(g.atoms[start], g) && rec(0, start);
  }
  // 수소를 명시 원자로 펼친 그래프 (SMARTS용)
  function explicitGraph(m) {
    var g = new Mol();
    m.atoms.forEach(function (a) { g.addAtom({ el: a.el, arom: a.arom, charge: a.charge }); });
    m.bonds.forEach(function (b) { g.addBond(b.a, b.b, b.order, b.arom); g.bonds[g.bonds.length - 1].inRing = b.inRing; });
    m.atoms.forEach(function (a, i) { for (var k = 0; k < a.h; k++) { var h = g.addAtom({ el: "H", charge: 0 }); g.addBond(i, h.idx, 1, false); } });
    g.atoms.forEach(function (a, i) {
      a.z = EL[a.el] ? EL[a.el][0] : 0;
      a.th = a.nbrs.filter(function (n) { return g.atoms[n.atom].el === "H"; }).length;
      a.tx = a.nbrs.length; a.deg = a.nbrs.length;
      a.val = a.nbrs.reduce(function (s, n) { return s + (n.bond.arom ? 1.5 : n.bond.order); }, 0);
      a.inRing = i < m.atoms.length ? m.atoms[i].inRing : false;
    });
    return g;
  }

  // ---------- Crippen (RDKit Data/Crippen.txt) ----------
  var CRIPPEN = [
    ["C1", "[CH4]", 0.1441, 2.503], ["C1", "[CH3]C", 0.1441, 2.503], ["C1", "[CH2](C)C", 0.1441, 2.503],
    ["C2", "[CH](C)(C)C", 0, 2.433], ["C2", "[C](C)(C)(C)C", 0, 2.433],
    ["C3", "[CH3][N,O,P,S,F,Cl,Br,I]", -0.2035, 2.753], ["C3", "[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]", -0.2035, 2.753],
    ["C4", "[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]", -0.2051, 2.731], ["C4", "[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]", -0.2051, 2.731],
    ["C5", "[C]=[!C;A;!#1]", -0.2783, 5.007],
    ["C6", "[CH2]=C", 0.1551, 3.513], ["C6", "[CH1](=C)[A;!#1]", 0.1551, 3.513], ["C6", "[CH0](=C)([A;!#1])[A;!#1]", 0.1551, 3.513], ["C6", "[C](=C)=C", 0.1551, 3.513],
    ["C7", "[CX2]#[A;!#1]", 0.0017, 3.888], ["C8", "[CH3]c", 0.08452, 2.464], ["C9", "[CH3]a", -0.1444, 2.412],
    ["C10", "[CH2X4]a", -0.0516, 2.488], ["C11", "[CHX4]a", 0.1193, 2.582], ["C12", "[CH0X4]a", -0.0967, 2.576],
    ["C13", "[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]", -0.5443, 4.041],
    ["C14", "[c][#9]", 0, 3.257], ["C15", "[c][#17]", 0.245, 3.564], ["C16", "[c][#35]", 0.198, 3.18], ["C17", "[c][#53]", 0, 3.104],
    ["C18", "[cH]", 0.1581, 3.35], ["C19", "[c](:a)(:a):a", 0.2955, 4.346], ["C20", "[c](:a)(:a)-a", 0.2713, 3.904],
    ["C21", "[c](:a)(:a)-C", 0.136, 3.509], ["C22", "[c](:a)(:a)-N", 0.4619, 4.067], ["C23", "[c](:a)(:a)-O", 0.5437, 3.853],
    ["C24", "[c](:a)(:a)-S", 0.1893, 2.673], ["C25", "[c](:a)(:a)=[C,N,O]", -0.8186, 3.135],
    ["C26", "[C](=C)(a)[A;!#1]", 0.264, 4.305], ["C26", "[C](=C)(c)a", 0.264, 4.305], ["C26", "[CH1](=C)a", 0.264, 4.305], ["C26", "[C]=c", 0.264, 4.305],
    ["C27", "[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]", 0.2148, 2.693], ["CS", "[#6]", 0.08129, 3.243],
    ["H1", "[#1][#6,#1]", 0.123, 1.057], ["H2", "[#1]O[CX4,c]", -0.2677, 1.395], ["H2", "[#1]O[!C;!N;!O;!S]", -0.2677, 1.395],
    ["H2", "[#1][!C;!N;!O]", -0.2677, 1.395], ["H3", "[#1][#7]", 0.2142, 0.9627], ["H3", "[#1]O[#7]", 0.2142, 0.9627],
    ["H4", "[#1]OC=[#6,#7,O,S]", 0.298, 1.805], ["H4", "[#1]O[O,S]", 0.298, 1.805], ["HS", "[#1]", 0.1125, 1.112],
    ["N1", "[NH2+0][A;!#1]", -1.019, 2.262], ["N2", "[NH+0]([A;!#1])[A;!#1]", -0.7096, 2.173], ["N3", "[NH2+0]a", -1.027, 2.827],
    ["N4", "[NH1+0]([!#1;A,a])a", -0.5188, 3], ["N5", "[NH+0]=[!#1;A,a]", 0.08387, 1.757], ["N6", "[N+0](=[!#1;A,a])[!#1;A,a]", 0.1836, 2.428],
    ["N7", "[N+0]([A;!#1])([A;!#1])[A;!#1]", -0.3187, 1.839], ["N8", "[N+0](a)([!#1;A,a])[A;!#1]", -0.4458, 2.819], ["N8", "[N+0](a)(a)a", -0.4458, 2.819],
    ["N9", "[N+0]#[A;!#1]", 0.01508, 1.725], ["N10", "[NH3,NH2,NH;+,+2,+3]", -1.95, 0], ["N11", "[n+0]", -0.3239, 2.202], ["N12", "[n;+,+2,+3]", -1.119, 0],
    ["N13", "[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]", -0.3396, 0.2604], ["N13", "[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]", -0.3396, 0.2604],
    ["N13", "[NH0;+,+2,+3](=[#6])=[#7]", -0.3396, 0.2604], ["N14", "[N;+,+2,+3]#[A;!#1]", 0.2887, 3.359], ["N14", "[N;-,-2,-3]", 0.2887, 3.359],
    ["N14", "[N;+,+2,+3](=[N;-,-2,-3])=N", 0.2887, 3.359], ["NS", "[#7]", -0.4806, 2.134],
    ["O1", "[o]", 0.1552, 1.08], ["O2", "[OH,OH2]", -0.2893, 0.8238], ["O3", "[O]([A;!#1])[A;!#1]", -0.0684, 1.085], ["O4", "[O](a)[!#1;A,a]", -0.4195, 1.182],
    ["O5", "[O]=[#7,#8]", 0.0335, 3.367], ["O5", "[OX1;-,-2,-3][#7]", 0.0335, 3.367], ["O6", "[OX1;-,-2,-2][#16]", -0.3339, 0.7774], ["O6", "[O;-0]=[#16;-0]", -0.3339, 0.7774],
    ["O12", "[O-]C(=O)", -1.326, 0], ["O7", "[OX1;-,-2,-3][!#1;!N;!S]", -1.189, 0], ["O8", "[O]=c", 0.1788, 3.135],
    ["O9", "[O]=[CH]C", -0.1526, 0], ["O9", "[O]=C(C)([A;!#1])", -0.1526, 0], ["O9", "[O]=[CH][N,O]", -0.1526, 0], ["O9", "[O]=[CH2]", -0.1526, 0], ["O9", "[O]=[CX2]=O", -0.1526, 0],
    ["O10", "[O]=[CH]c", 0.1129, 0.2215], ["O10", "[O]=C([C,c])[a;!#1]", 0.1129, 0.2215], ["O10", "[O]=C(c)[A;!#1]", 0.1129, 0.2215],
    ["O11", "[O]=C([!#1;!#6])[!#1;!#6]", 0.4833, 0.389], ["OS", "[#8]", -0.1188, 0.6865],
    ["F", "[#9-0]", 0.4202, 1.108], ["Cl", "[#17-0]", 0.6895, 5.853], ["Br", "[#35-0]", 0.8456, 8.927], ["I", "[#53-0]", 0.8857, 14.02],
    ["Hal", "[#9,#17,#35,#53;-]", -2.996, 0], ["Hal", "[#53;+,+2,+3]", -2.996, 0], ["Hal", "[+;#3,#11,#19,#37,#55]", -2.996, 0],
    ["P", "[#15]", 0.8612, 6.92], ["S2", "[S;-,-2,-3,-4,+1,+2,+3,+5,+6]", -0.0024, 7.365], ["S2", "[S-0]=[N,O,P,S]", -0.0024, 7.365],
    ["S1", "[S;A]", 0.6482, 7.591], ["S3", "[s;a]", 0.6237, 6.691],
    ["Me1", "[#3,#11,#19,#37,#55]", -0.3808, 5.754], ["Me1", "[#4,#12,#20,#38,#56]", -0.3808, 5.754], ["Me1", "[#5,#13,#31,#49,#81]", -0.3808, 5.754],
    ["Me1", "[#14,#32,#50,#82]", -0.3808, 5.754], ["Me1", "[#33,#51,#83]", -0.3808, 5.754], ["Me1", "[#34,#52,#84]", -0.3808, 5.754],
    ["Me2", "[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]", -0.0025, 0], ["Me2", "[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]", -0.0025, 0],
    ["Me2", "[#72,#73,#74,#75,#76,#77,#78,#79,#80]", -0.0025, 0]
  ].map(function (r) { return { id: r[0], q: parseSmarts(r[1]), logp: r[2], mr: r[3] }; });
  function crippen(m) {
    var g = explicitGraph(m), lp = 0, mr = 0;
    g.atoms.forEach(function (a, i) {
      for (var k = 0; k < CRIPPEN.length; k++) { if (matchAt(CRIPPEN[k].q, g, i)) { lp += CRIPPEN[k].logp; mr += CRIPPEN[k].mr; return; } }
    });
    return { logp: lp, mr: mr };
  }

  // ---------- TPSA (RDKit MolSurf.cpp, N·O만) ----------
  function tpsa(m) {
    var res = 0;
    m.atoms.forEach(function (a) {
      if (a.el !== "N" && a.el !== "O") return;
      var nN = a.nbrs.length, nH = a.h, ch = a.charge, nS = 0, nD = 0, nT = 0, nA = 0;
      a.nbrs.forEach(function (nb) { var b = nb.bond; if (b.arom) nA++; else if (b.order === 1) nS++; else if (b.order === 2) nD++; else if (b.order === 3) nT++; });
      var in3 = a.ringSizes.indexOf(3) >= 0, t = -1;
      if (a.el === "N") {
        if (nN === 1) {
          if (nH === 0 && ch === 0 && nT === 1) t = 23.79; else if (nH === 1 && ch === 0 && nD === 1) t = 23.85;
          else if (nH === 2 && ch === 0 && nS === 1) t = 26.02; else if (nH === 2 && ch === 1 && nD === 1) t = 25.59;
          else if (nH === 3 && ch === 1 && nS === 1) t = 27.64;
        } else if (nN === 2) {
          if (nH === 0 && ch === 0 && nS === 1 && nD === 1) t = 12.36; else if (nH === 0 && ch === 0 && nT === 1 && nD === 1) t = 13.60;
          else if (nH === 1 && ch === 0 && nS === 2 && in3) t = 21.94; else if (nH === 1 && ch === 0 && nS === 2 && !in3) t = 12.03;
          else if (nH === 0 && ch === 1 && nT === 1 && nS === 1) t = 4.36; else if (nH === 1 && ch === 1 && nD === 1 && nS === 1) t = 13.97;
          else if (nH === 2 && ch === 1 && nS === 2) t = 16.61; else if (nH === 0 && ch === 0 && nA === 2) t = 12.89;
          else if (nH === 1 && ch === 0 && nA === 2) t = 15.79; else if (nH === 1 && ch === 1 && nA === 2) t = 14.14;
        } else if (nN === 3) {
          if (nH === 0 && ch === 0 && nS === 3 && in3) t = 3.01; else if (nH === 0 && ch === 0 && nS === 3 && !in3) t = 3.24;
          else if (nH === 0 && ch === 0 && nS === 1 && nD === 2) t = 11.68; else if (nH === 0 && ch === 1 && nS === 2 && nD === 1) t = 3.01;
          else if (nH === 1 && ch === 1 && nS === 3) t = 4.44; else if (nH === 0 && ch === 0 && nA === 3) t = 4.41;
          else if (nH === 0 && ch === 0 && nS === 1 && nA === 2) t = 4.93; else if (nH === 0 && ch === 0 && nD === 1 && nA === 2) t = 8.39;
          else if (nH === 0 && ch === 1 && nA === 3) t = 4.10; else if (nH === 0 && ch === 1 && nS === 1 && nA === 2) t = 3.88;
        } else if (nN === 4) { if (nH === 0 && nS === 4 && ch === 1) t = 0; }
        if (t < 0) t = Math.max(0, 30.5 - nN * 8.2 + nH * 1.5);
      } else {
        if (nN === 1) {
          if (nH === 0 && ch === 0 && nD === 1) t = 17.07; else if (nH === 1 && ch === 0 && nS === 1) t = 20.23;
          else if (nH === 0 && ch === -1 && nS === 1) t = 23.06;
        } else if (nN === 2) {
          if (nH === 0 && ch === 0 && nS === 2 && in3) t = 12.53; else if (nH === 0 && ch === 0 && nS === 2 && !in3) t = 9.23;
          else if (nH === 0 && ch === 0 && nA === 2) t = 13.14;
        }
        if (t < 0) t = Math.max(0, 28.5 - nN * 8.6 + nH * 1.5);
      }
      res += t;
    });
    return res;
  }

  // ---------- 수소결합 주개/받개, 회전 결합 ----------
  function valence(m, a) { return a.nbrs.reduce(function (s, n) { return s + (n.bond.arom ? 1.5 : n.bond.order); }, 0) + a.h; }
  function hasDoubleTo(m, a, els, acyclicOnly) {
    return a.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2 && els.indexOf(m.atoms[n.atom].el) >= 0 && (!acyclicOnly || !n.bond.inRing); });
  }
  function hbd(m) {
    return m.atoms.filter(function (a) {
      var v = Math.round(valence(m, a));
      if (a.el === "N" && !a.arom) return a.h > 0 && ((a.charge === 0 && v === 3) || (a.charge === 1 && v === 4));
      if ((a.el === "O" || a.el === "S") && !a.arom) return a.h === 1 && a.charge === 0;
      if (a.el === "N" && a.arom) return a.h === 1 && a.charge === 0;
      return false;
    }).length;
  }
  function hba(m) {
    return m.atoms.filter(function (a) {
      var v = Math.round(valence(m, a));
      if ((a.el === "O" || a.el === "S") && !a.arom) {
        if (a.charge < 0) return true;
        if (v !== 2) return false;
        if (a.h === 1) return !a.nbrs.some(function (n) { return hasDoubleTo(m, m.atoms[n.atom], ["O", "N", "P", "S"], false); });
        return a.h === 0;
      }
      if (a.el === "N" && !a.arom) {
        if (v !== 3) return false;
        return !a.nbrs.some(function (n) { return n.bond.order === 1 && !n.bond.arom && hasDoubleTo(m, m.atoms[n.atom], ["O", "N", "P", "S"], true); });
      }
      if (a.arom && a.charge === 0) { if (a.el === "N") return a.h === 0 && a.nbrs.length === 2; return a.el === "O" || a.el === "S"; }
      return false;
    }).length;
  }
  function lipinskiCounts(m) {
    var no = 0, nhoh = 0;
    m.atoms.forEach(function (a) { if (a.el === "N" || a.el === "O") { no++; nhoh += a.h; } });
    return { no: no, nhoh: nhoh };
  }
  function rotatable(m) {
    function heavyDeg(a) { return a.nbrs.length; }
    function inTriple(a) { return a.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 3; }); }
    function isCX3(a, X) { return a.el === "C" && a.nbrs.filter(function (n) { return m.atoms[n.atom].el === X; }).length === 3; }
    function isTBu(a) { return a.el === "C" && a.nbrs.length === 4 && a.nbrs.filter(function (n) { var b = m.atoms[n.atom]; return b.el === "C" && b.h === 3 && b.nbrs.length === 1; }).length === 3; }
    function bad(a) { return heavyDeg(a) === 1 || inTriple(a) || isCX3(a, "F") || isCX3(a, "Cl") || isCX3(a, "Br") || isTBu(a) || (a.el === "C" && a.h === 3); }
    function amideLike(c, x) { // [CD3](=[N,O,S])-!@[#7,O,S!D1]
      if (c.el !== "C" || c.nbrs.length !== 3) return false;
      if (!c.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2 && ["N", "O", "S"].indexOf(m.atoms[n.atom].el) >= 0; })) return false;
      return x.el === "N" || x.el === "O" || (x.el === "S" && x.nbrs.length !== 1);
    }
    function amidinium(c, x) { return c.el === "C" && c.nbrs.length === 3 && c.nbrs.some(function (n) { var t = m.atoms[n.atom]; return !n.bond.arom && n.bond.order === 2 && t.el === "N" && t.charge === 1; }) && x.el === "N" && x.nbrs.length !== 1; }
    return m.bonds.filter(function (b) {
      if (b.inRing || b.order !== 1 || b.arom) return false;
      var a = m.atoms[b.a], c = m.atoms[b.b];
      if (bad(a) || bad(c)) return false;
      if (amideLike(a, c) || amideLike(c, a) || amidinium(a, c) || amidinium(c, a)) return false;
      return true;
    }).length;
  }

  // ---------- 지문(ECFP4 유사)·유사도 ----------
  function hash32(arr) { var h = 2166136261; for (var i = 0; i < arr.length; i++) { h ^= arr[i] & 0xffffffff; h = Math.imul(h, 16777619); h ^= h >>> 13; } return h >>> 0; }
  function morgan(m, radius, nbits) {
    radius = radius || 2; nbits = nbits || 2048;
    var ids = m.atoms.map(function (a) { return hash32([EL[a.el] ? EL[a.el][0] : 0, a.nbrs.length, a.h, a.charge + 8, a.inRing ? 1 : 0, a.arom ? 1 : 0]); });
    var bits = new Set(); ids.forEach(function (x) { bits.add(x % nbits); });
    for (var r = 0; r < radius; r++) {
      ids = m.atoms.map(function (a, i) {
        var env = a.nbrs.map(function (n) { return [n.bond.arom ? 4 : n.bond.order, ids[n.atom]]; }).sort(function (p, q) { return p[0] - q[0] || p[1] - q[1]; });
        var arr = [r + 1, ids[i]]; env.forEach(function (e) { arr.push(e[0], e[1]); }); return hash32(arr);
      });
      ids.forEach(function (x) { bits.add(x % nbits); });
    }
    return bits;
  }
  // 같은 물질인지 가려내는 열쇠값. 원자 환경을 반복해서 다듬은 뒤 정렬해 해시하므로
  // 같은 구조를 다시 그리거나 원자 순서가 달라져도 같은 값이 나와요.
  function structKey(m) {
    if (!m || !m.atoms.length) return "";
    var ids = m.atoms.map(function (a) {
      return hash32([EL[a.el] ? EL[a.el][0] : 0, a.nbrs.length, a.h, a.charge + 8, a.arom ? 1 : 0]);
    });
    for (var r = 0; r < m.atoms.length + 2 && r < 12; r++) {
      ids = m.atoms.map(function (a, i) {
        var env = a.nbrs.map(function (n) { return [n.bond.arom ? 4 : n.bond.order, ids[n.atom]]; })
          .sort(function (p, q) { return p[0] - q[0] || p[1] - q[1]; });
        var arr = [ids[i]]; env.forEach(function (e) { arr.push(e[0], e[1]); });
        return hash32(arr);
      });
    }
    var sorted = ids.slice().sort(function (a, b) { return a - b; });
    var ch = m.atoms.reduce(function (s, a) { return s + a.charge; }, 0);
    return hash32(sorted.concat([m.atoms.length, m.bonds.length, ch + 64])).toString(36) + "-" +
      hash32(sorted.map(function (x, i) { return (x ^ (i * 2654435761)) >>> 0; })).toString(36);
  }
  function tanimoto(a, b) { var inter = 0; a.forEach(function (x) { if (b.has(x)) inter++; }); var u = a.size + b.size - inter; return u ? inter / u : 0; }

  // ---------- 부분 구조 (질의 분자 ⊂ 대상) ----------
  function substructure(query, target) {
    var qn = query.atoms.length, tn = target.atoms.length; if (qn > tn) return false;
    var map = new Array(qn).fill(-1), used = new Array(tn).fill(false);
    // BFS 순서
    var order = [], seen = new Array(qn).fill(false);
    for (var s = 0; s < qn; s++) { if (seen[s]) continue; var q = [s]; seen[s] = true; while (q.length) { var u = q.shift(); order.push(u); query.atoms[u].nbrs.forEach(function (n) { if (!seen[n.atom]) { seen[n.atom] = true; q.push(n.atom); } }); } }
    function atomOk(qa, ta) { return qa.el === ta.el && (qa.charge === 0 || qa.charge === ta.charge) && (!qa.arom || ta.arom) && ta.nbrs.length >= qa.nbrs.length; }
    function bondOk2(qb, tb) { if (qb.arom || tb.arom) return qb.arom === tb.arom; return qb.order === tb.order; }
    function rec(k) {
      if (k === order.length) return true;
      var qi = order[k], qa = query.atoms[qi];
      for (var ti = 0; ti < tn; ti++) {
        if (used[ti] || !atomOk(qa, target.atoms[ti])) continue;
        var ok = qa.nbrs.every(function (n) {
          var mj = map[n.atom]; if (mj < 0) return true;
          var tb = null; target.atoms[ti].nbrs.forEach(function (tn2) { if (tn2.atom === mj) tb = tn2.bond; });
          return tb && bondOk2(n.bond, tb);
        });
        if (!ok) continue;
        map[qi] = ti; used[ti] = true;
        if (rec(k + 1)) return true;
        map[qi] = -1; used[ti] = false;
      }
      return false;
    }
    return rec(0);
  }

  // ---------- 종합 ----------
  function props(m) {
    var c = formulaCounts(m), mw = 0, ex = 0;
    Object.keys(c).forEach(function (k) { mw += (EL[k] ? EL[k][1] : 0) * c[k]; ex += (EL[k] ? EL[k][2] : 0) * c[k]; });
    var ch = netCharge(m); ex -= ch * 0.000548579909;
    var cr = crippen(m), lc = lipinskiCounts(m);
    var heavy = m.atoms.length, sp3C = 0, nC = 0;
    m.atoms.forEach(function (a) {
      if (a.el !== "C") return; nC++;
      if (!a.arom && a.nbrs.every(function (n) { return !n.bond.arom && n.bond.order === 1; })) sp3C++;
    });
    var aromRings = m.rings.filter(function (r) { return r.every(function (i) { return m.atoms[i].arom; }); }).length;
    var p = {
      formula: hill(c) + (ch ? (ch > 0 ? (ch > 1 ? ch : "") + "+" : (ch < -1 ? -ch : "") + "-") : ""),
      mw: mw, exact: ex, charge: ch, heavy: heavy,
      clogp: cr.logp, mr: cr.mr, tpsa: tpsa(m), hbd: hbd(m), hba: hba(m), nhoh: lc.nhoh, no: lc.no,
      rotb: rotatable(m), rings: m.rings.length, aromRings: aromRings, fsp3: nC ? sp3C / nC : 0,
      hetero: m.atoms.filter(function (a) { return a.el !== "C" && a.el !== "H"; }).length
    };
    var v = 0; if (p.mw > 500) v++; if (p.clogp > 5) v++; if (p.nhoh > 5) v++; if (p.no > 10) v++;
    p.lipinskiViolations = v; p.lipinski = v <= 1; p.veber = p.rotb <= 10 && p.tpsa <= 140;
    return p;
  }
  function fromAny(txt) {
    txt = String(txt || "").trim(); if (!txt) throw new Error("구조가 비어 있어요");
    if (/V2000|V3000|M  END/.test(txt)) return parseMolfile(txt);
    return parseSmiles(txt);
  }

  // ---------- 이온화 작용기 찾기 (CNS MPO의 pKa·logD 추정용) ----------
  function ionizable(m) {
    var bases = [], acids = [];
    function nb(a) { return a.nbrs.map(function (n) { return { at: m.atoms[n.atom], bond: n.bond }; }); }
    function dbl(a, els) { return a.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2 && els.indexOf(m.atoms[n.atom].el) >= 0; }); }
    m.atoms.forEach(function (a) {
      // 염기: 지방족 아민 (아닐린·아마이드·설폰아마이드·카바메이트·하이드라진·N-옥시 제외)
      if (a.el === "N" && !a.arom && (a.charge === 0 || (a.charge === 1 && a.h > 0))) {
        var ns = nb(a);
        if (ns.every(function (x) { return !x.bond.arom && x.bond.order === 1; })) {
          var weak = ns.some(function (x) {
            var t = x.at;
            if (t.arom) return true;                                  // 아닐린
            if (t.el === "N" || t.el === "O" || t.el === "S" || t.el === "P") return true;
            if (t.el === "C" && (dbl(t, ["O", "S", "N"]) || t.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 3; }))) return true; // 아마이드·아미딘 질소는 아래에서
            if (t.el === "C" && t.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2; })) return true; // 엔아민
            return false;
          });
          if (!weak) bases.push({ atom: a.idx, kind: a.charge === 1 ? "양성자화 아민" : (a.h === 2 ? "1차 아민" : a.h === 1 ? "2차 아민" : "3차 아민") });
        }
        if (a.charge === 1 && a.h === 0 && a.nbrs.length === 4) bases.push({ atom: a.idx, kind: "4차 암모늄", quat: true });
      }
      // 염기: 아미딘·구아니딘 (질소에 카보닐이 붙지 않은 경우)
      if (a.el === "C" && !a.arom) {
        var dN = a.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 2 && m.atoms[n.atom].el === "N"; });
        var sN = a.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 1 && m.atoms[n.atom].el === "N"; });
        if (dN.length === 1 && sN.length >= 1) {
          var acyl = dN.concat(sN).some(function (n) { var t = m.atoms[n.atom]; return t.nbrs.some(function (q) { var u = m.atoms[q.atom]; return u !== a && (dbl(u, ["O", "S"]) || u.el === "S"); }); });
          var dNat = m.atoms[dN[0].atom];
          if (!acyl && !dNat.nbrs.some(function (q) { return m.atoms[q.atom].el === "O"; })) bases.push({ atom: a.idx, kind: sN.length >= 2 ? "구아니딘" : "아미딘" });
        }
        // 산: 카복실산
        if (dbl(a, ["O"]) && a.nbrs.some(function (n) { var t = m.atoms[n.atom]; return t.el === "O" && !n.bond.arom && n.bond.order === 1 && (t.h === 1 || t.charge === -1) && t.nbrs.length === 1; })) acids.push({ atom: a.idx, kind: "카복실산" });
      }
      // 산: 설폰산·포스폰산
      if ((a.el === "S" || a.el === "P") && dbl(a, ["O"]) && a.nbrs.some(function (n) { var t = m.atoms[n.atom]; return t.el === "O" && n.bond.order === 1 && (t.h === 1 || t.charge === -1); })) acids.push({ atom: a.idx, kind: a.el === "S" ? "설폰산" : "포스폰산" });
      // 산: 아실설폰아마이드 NH
      if (a.el === "N" && !a.arom && a.h === 1) {
        var ns2 = nb(a);
        var hasSO2 = ns2.some(function (x) { return x.at.el === "S" && x.at.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 2 && m.atoms[n.atom].el === "O"; }).length === 2; });
        var hasCO = ns2.some(function (x) { return x.at.el === "C" && dbl(x.at, ["O"]); });
        if (hasSO2 && hasCO) acids.push({ atom: a.idx, kind: "아실설폰아마이드" });
      }
    });
    // 산: 테트라졸 NH
    (m.rings || []).forEach(function (r) {
      if (r.length !== 5) return;
      var nN = r.filter(function (i) { return m.atoms[i].el === "N"; }).length;
      if (nN === 4 && r.some(function (i) { return m.atoms[i].el === "N" && m.atoms[i].h === 1; })) acids.push({ atom: r[0], kind: "테트라졸" });
    });
    return { bases: bases, acids: acids };
  }

  // ---------- CNS MPO (Wager et al., ACS Chem. Neurosci. 2010, 1, 435) ----------
  function down(x, lo, hi) { return x <= lo ? 1 : x >= hi ? 0 : (hi - x) / (hi - lo); }
  function cnsMpo(o) {
    var parts = {
      clogp: o.clogp == null ? null : down(o.clogp, 3, 5),
      clogd: o.clogd == null ? null : down(o.clogd, 2, 4),
      mw: o.mw == null ? null : down(o.mw, 360, 500),
      tpsa: o.tpsa == null ? null : (o.tpsa < 20 ? 0 : o.tpsa < 40 ? (o.tpsa - 20) / 20 : o.tpsa <= 90 ? 1 : o.tpsa < 120 ? (120 - o.tpsa) / 30 : 0),
      hbd: o.hbd == null ? null : down(o.hbd, 0.5, 3.5),
      pka: o.pka == null ? null : down(o.pka, 8, 10)
    };
    var known = 0, sum = 0, missing = [];
    Object.keys(parts).forEach(function (k) { if (parts[k] == null) missing.push(k); else { sum += parts[k]; known++; } });
    return { parts: parts, min: sum, max: sum + missing.length, missing: missing, total: missing.length ? null : sum };
  }

  // ---------- MOL 파일을 떨어져 있는 덩어리별로 나누기 ----------
  function pad(n, w) { var t = String(n); while (t.length < w) t = " " + t; return t; }
  function molParts(txt) { // MOL 본문을 원자·결합·속성 줄로 쪼개기
    var L = String(txt).replace(/\r/g, "").split("\n"), ci = -1;
    for (var k = 0; k < Math.min(L.length, 10); k++) if (/V2000/.test(L[k])) { ci = k; break; }
    if (ci < 0) return null;
    var na = +L[ci].substr(0, 3), nb = +L[ci].substr(3, 3);
    var props = [];
    for (var q = ci + 1 + na + nb; q < L.length; q++) { if (/^M  END/.test(L[q])) break; props.push(L[q]); }
    return { head: L.slice(0, ci), na: na, nb: nb, atoms: L.slice(ci + 1, ci + 1 + na), bonds: L.slice(ci + 1 + na, ci + 1 + na + nb), props: props };
  }
  function molfileComponents(txt) { // 떨어져 있는 덩어리마다 원자 번호 목록
    var P = molParts(txt); if (!P || !P.na) return [];
    var par = []; for (var i = 0; i < P.na; i++) par[i] = i;
    function find(x) { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; }
    P.bonds.forEach(function (t) { var a = +t.substr(0, 3) - 1, b = +t.substr(3, 3) - 1, x = find(a), y = find(b); if (x !== y) par[y] = x; });
    var g = {}, order = [];
    for (var j = 0; j < P.na; j++) { var r = find(j); if (!g[r]) { g[r] = []; order.push(r); } g[r].push(j); }
    return order.map(function (k) { return g[k]; });
  }
  function molfileSubset(txt, idx) { // 고른 원자만으로 새 MOL 만들기
    var P = molParts(txt); if (!P) return txt;
    var map = {}; idx.forEach(function (a, n) { map[a] = n + 1; });
    var bs = P.bonds.filter(function (t) { return map[+t.substr(0, 3) - 1] && map[+t.substr(3, 3) - 1]; });
    var out = P.head.slice();
    out.push(pad(idx.length, 3) + pad(bs.length, 3) + "  0  0  0  0  0  0  0  0999 V2000");
    idx.forEach(function (a) { out.push(P.atoms[a]); });
    bs.forEach(function (t) { out.push(pad(map[+t.substr(0, 3) - 1], 3) + pad(map[+t.substr(3, 3) - 1], 3) + t.substr(6)); });
    P.props.forEach(function (t) {
      var mm = t.match(/^M  (CHG|ISO|RAD)/); if (!mm) return;
      var n = +t.substr(6, 3), items = [];
      for (var z = 0; z < n; z++) { var ai = +t.substr(10 + 8 * z, 3) - 1, v = +t.substr(14 + 8 * z, 3); if (map[ai]) items.push([map[ai], v]); }
      while (items.length) {
        var chunk = items.splice(0, 8);
        out.push("M  " + mm[1] + pad(chunk.length, 3) + chunk.map(function (c) { return pad(c[0], 4) + pad(c[1], 4); }).join(""));
      }
    });
    out.push("M  END");
    return out.join("\n");
  }
  function splitMolfile(txt) {
    var cs = molfileComponents(txt);
    if (cs.length <= 1) return [txt];
    return cs.map(function (c) { return molfileSubset(txt, c); });
  }

  // 흔한 염·카운터이온 (분자식으로 판정)
  var SALT_FORMULA = { HCl: 1, HBr: 1, HI: 1, HF: 1, "H2SO4": 1, "H3PO4": 1, "HNO3": 1, "H2O": 1, "CH4O3S": 1, "C2HF3O2": 1, "C2H4O2": 1, "C4H4O4": 1, "C4H6O4": 1, "C4H6O6": 1, "C6H8O7": 1, "C7H6O3": 1, "C7H8O3S": 1, "C2H2O4": 1, "CHNO3": 1, "Na": 1, "K": 1, "Li": 1, "Ca": 1, "Mg": 1, "Zn": 1, "Cl": 1, "Br": 1, "I": 1, "NH3": 1, "C2H6O": 1, "C3H8O": 1, "CH4O": 1, "C6H14O3S": 1 };
  function isSalt(mol) {
    var heavy = mol.atoms.length;
    if (heavy > 12) return false;
    var p = props(mol), f = p.formula.replace(/[+-]\d*$|\d*[+-]$/, "");
    if (SALT_FORMULA[f]) return true;
    if (heavy <= 2) return true;                       // 단일 이온·작은 분자
    if (p.charge !== 0 && heavy <= 8) return true;     // 작은 이온성 조각
    return false;
  }

  root.Chem = { splitMolfile: splitMolfile, molfileComponents: molfileComponents, molfileSubset: molfileSubset, isSalt: isSalt, cnsMpo: cnsMpo, ionizable: ionizable, parseSmiles: parseSmiles, parseMolfile: parseMolfile, fromAny: fromAny, props: props, morgan: morgan, structKey: structKey, tanimoto: tanimoto, substructure: substructure, crippen: crippen, tpsa: tpsa };
  if (typeof module !== "undefined") module.exports = root.Chem;
})(typeof window !== "undefined" ? window : globalThis);
