/* 의약화학 오피스 — 3D 좌표 만들기
 *
 * 편집기의 2D 좌표를 출발점으로 삼고, 쐐기(wedge)·점선(hash) 결합에서 z를 뽑아
 * 입체배치를 살린 뒤, 간단한 힘장으로 풀어서 3D 모양을 만듭니다.
 * 결합 길이·결합각·sp2 평면성·이중결합 회전 억제·원자끼리 밀어내기만 넣은
 * 가벼운 계산이라, 보기용으로는 충분하지만 도킹에 바로 쓸 품질은 아니에요.
 */
(function (root) {
  "use strict";
  var Chem = root.Chem;

  // 공유 결합 반지름 (Å)
  var RCOV = { H: 0.31, B: 0.84, C: 0.76, N: 0.71, O: 0.66, F: 0.57, Si: 1.11, P: 1.07, S: 1.05,
    Cl: 1.02, Br: 1.20, I: 1.39, Se: 1.20, Li: 1.28, Na: 1.66, K: 2.03, Mg: 1.41, Ca: 1.76, Zn: 1.22, Fe: 1.32 };
  // 반데르발스 반지름 (밀어내기용)
  var RVDW = { H: 1.10, C: 1.70, N: 1.55, O: 1.52, F: 1.47, P: 1.80, S: 1.80, Cl: 1.75, Br: 1.85, I: 1.98 };
  function rcov(el) { return RCOV[el] || 0.95; }
  function rvdw(el) { return RVDW[el] || 1.7; }

  function bondLen(m, b) {
    var f = b.arom ? 0.93 : b.order === 3 ? 0.84 : b.order === 2 ? 0.87 : 1.0;
    return (rcov(m.atoms[b.a].el) + rcov(m.atoms[b.b].el)) * f;
  }
  // 원자 주변이 몇 도로 벌어지는지
  function idealAngle(m, i) {
    var a = m.atoms[i], n = a.nbrs.length + (a.h || 0);
    var hasTriple = a.nbrs.some(function (x) { return !x.bond.arom && x.bond.order === 3; });
    var hasDouble = a.nbrs.some(function (x) { return !x.bond.arom && x.bond.order === 2; });
    var dbl = a.nbrs.filter(function (x) { return !x.bond.arom && x.bond.order === 2; }).length;
    if (hasTriple || dbl >= 2) return Math.PI;                    // sp
    if (a.arom || hasDouble) return 120 * Math.PI / 180;          // sp2
    if (n <= 2 && (a.el === "O" || a.el === "S")) return 104.5 * Math.PI / 180;
    if (a.el === "N" && n === 3 && !a.arom && !hasDouble) return 107.5 * Math.PI / 180;
    return 109.47 * Math.PI / 180;                                // sp3
  }
  function isPlanar(m, i) {
    var a = m.atoms[i];
    if (a.arom) return true;
    return a.nbrs.some(function (x) { return !x.bond.arom && x.bond.order === 2; });
  }

  // ---------- 수소 붙이기 ----------
  function withHydrogens(m) {
    var out = { atoms: [], bonds: [] };                 // 원본은 건드리지 않고 복제해서 씁니다
    m.atoms.forEach(function (a) {
      out.atoms.push({ el: a.el, charge: a.charge || 0, arom: !!a.arom, h: a.h || 0,
        x: a.x || 0, y: a.y || 0, z: 0, nbrs: [], heavy: true, src: a.idx });
    });
    m.bonds.forEach(function (b) {
      out.bonds.push({ a: b.a, b: b.b, order: b.order, arom: !!b.arom, stereo: b.stereo || 0, from: b.from == null ? b.a : b.from });
    });
    // 암시적 수소를 실제 원자로
    var n0 = out.atoms.length;
    for (var i = 0; i < n0; i++) {
      var h = out.atoms[i].h | 0;
      for (var k = 0; k < h; k++) {
        out.atoms.push({ el: "H", charge: 0, arom: false, h: 0, x: 0, y: 0, z: 0, nbrs: [], heavy: false, src: -1 });
        out.bonds.push({ a: i, b: out.atoms.length - 1, order: 1, arom: false, stereo: 0 });
      }
    }
    out.bonds.forEach(function (b, bi) {
      out.atoms[b.a].nbrs.push({ atom: b.b, bi: bi });
      out.atoms[b.b].nbrs.push({ atom: b.a, bi: bi });
    });
    return out;
  }

  // ---------- 시작 좌표 ----------
  function seed(g, m) {
    var rnd0 = 7;
    function r0() { rnd0 = (rnd0 * 1103515245 + 12345) & 0x7fffffff; return rnd0 / 0x7fffffff - 0.5; }

    var xs = [], ys = [];
    g.atoms.forEach(function (a) { if (a.heavy) { xs.push(a.x); ys.push(a.y); } });
    var span = xs.length ? Math.max(Math.max.apply(null, xs) - Math.min.apply(null, xs),
                                    Math.max.apply(null, ys) - Math.min.apply(null, ys)) : 0;
    if (span < 0.01) {
      // 2D 좌표가 없는 경우 (SMILES 에서 바로 온 구조) — 이웃을 따라가며 퍼뜨립니다
      var seen = {}, order = [], q = [0];
      seen[0] = 1;
      while (q.length) { var u = q.shift(); order.push(u); g.atoms[u].nbrs.forEach(function (nb) { if (!seen[nb.atom]) { seen[nb.atom] = 1; q.push(nb.atom); } }); }
      g.atoms.forEach(function (a, i) { if (!seen[i]) order.push(i); });
      order.forEach(function (i, k) {
        var a = g.atoms[i];
        if (k === 0) { a.x = a.y = a.z = 0; return; }
        var p = null;
        a.nbrs.forEach(function (nb) { if (p == null && seen[nb.atom] && order.indexOf(nb.atom) < k) p = g.atoms[nb.atom]; });
        if (!p) p = g.atoms[order[k - 1]];
        a.x = p.x + r0() * 2.6; a.y = p.y + r0() * 2.6; a.z = p.z + r0() * 2.6;
      });
      return;
    }
    var scale = 1.5;                                   // 편집기 단위 → Å 비슷하게
    g.atoms.forEach(function (a) { a.x *= scale; a.y *= scale; });

    // 쐐기·점선에서 z를 뽑아 입체배치를 살립니다
    g.pin = {};
    g.bonds.forEach(function (b) {
      if (b.stereo !== 1 && b.stereo !== 6) return;
      var ctr = b.from === b.a ? b.a : b.b;              // 좁은 끝 = 입체중심
      var wide = ctr === b.a ? b.b : b.a;                // 넓은 끝 = 떠오르는 쪽
      var dz = b.stereo === 1 ? 0.95 : -0.95;            // 1 = 쐐기(앞), 6 = 점선(뒤)
      g.atoms[wide].z += dz;
      g.pin[wide] = dz > 0 ? 1 : -1;                     // 계산 도중 반대쪽으로 넘어가지 않게
      // 입체중심에 수소가 하나 있으면 반대편으로 보냅니다.
      // 이걸 안 하면 수소가 쐐기와 같은 쪽으로 가버려서 배치가 뒤집혀요.
      g.atoms[ctr].nbrs.forEach(function (nb) {
        if (g.atoms[nb.atom].heavy || g.pin[nb.atom] != null) return;
        g.atoms[nb.atom].z = -dz;
        g.pin[nb.atom] = dz > 0 ? -1 : 1;
      });
    });
    // 평면에서 벗어나게 살짝 흔들어 줍니다 (완전 평면은 풀리지 않아요)
    var rnd = 1;
    function rand() { rnd = (rnd * 1103515245 + 12345) & 0x7fffffff; return rnd / 0x7fffffff - 0.5; }
    g.atoms.forEach(function (a) {
      a.z += rand() * 0.35;
      if (!a.heavy) { a.x += rand() * 0.6; a.y += rand() * 0.6; a.z += rand() * 0.6; }
    });
    // 수소는 붙어 있는 원자 근처에서 출발
    g.atoms.forEach(function (a, i) {
      if (a.heavy || !a.nbrs.length) return;
      var p = g.atoms[a.nbrs[0].atom];
      a.x += p.x; a.y += p.y;
      if (g.pin[i] == null) a.z += p.z;                  // 핀이 걸린 수소는 정해둔 쪽에 그대로
    });
  }

  // ---------- 힘장 ----------
  function build(g, m) {
    var bonds = [], angles = [], planes = [], torsions = [];
    g.bonds.forEach(function (b, bi) {
      var fake = { a: b.a, b: b.b, order: b.order, arom: b.arom };
      var ea = g.atoms[b.a].el, eb = g.atoms[b.b].el;
      var f = b.arom ? 0.93 : b.order === 3 ? 0.84 : b.order === 2 ? 0.87 : 1.0;
      bonds.push({ i: b.a, j: b.b, r0: (rcov(ea) + rcov(eb)) * f, bi: bi });
    });
    g.atoms.forEach(function (a, i) {
      if (a.nbrs.length < 2) return;
      var t0 = idealAngleG(g, i);
      for (var p = 0; p < a.nbrs.length; p++)
        for (var q = p + 1; q < a.nbrs.length; q++)
          angles.push({ i: a.nbrs[p].atom, j: i, k: a.nbrs[q].atom, t0: a.nbrs.length >= 4 ? 109.47 * Math.PI / 180 : t0 });
      if (a.nbrs.length === 3 && planarG(g, i))
        planes.push({ c: i, a: a.nbrs[0].atom, b: a.nbrs[1].atom, d: a.nbrs[2].atom });
    });
    // 이중결합 회전 막기 — 양 끝 치환기가 한 평면에 있도록
    g.bonds.forEach(function (b) {
      if (b.arom || b.order !== 2) return;
      var A = g.atoms[b.a], B = g.atoms[b.b];
      var na = A.nbrs.filter(function (x) { return x.atom !== b.b; });
      var nb = B.nbrs.filter(function (x) { return x.atom !== b.a; });
      na.forEach(function (x) { nb.forEach(function (y) { torsions.push({ i: x.atom, j: b.a, k: b.b, l: y.atom }); }); });
    });
    // 방향족 고리도 평면 유지
    (m.rings || []).forEach(function (r) {
      if (!r.every(function (i) { return m.atoms[i].arom; })) return;
      for (var i = 0; i + 3 < r.length + 3; i++) {
        var p0 = r[i % r.length], p1 = r[(i + 1) % r.length], p2 = r[(i + 2) % r.length], p3 = r[(i + 3) % r.length];
        torsions.push({ i: p0, j: p1, k: p2, l: p3 });
      }
    });
    return { bonds: bonds, angles: angles, planes: planes, torsions: torsions };
  }
  function idealAngleG(g, i) {
    var a = g.atoms[i], n = a.nbrs.length;
    var dbl = 0, tri = 0;
    a.nbrs.forEach(function (x) { var b = g.bonds[x.bi]; if (b.arom) return; if (b.order === 2) dbl++; if (b.order === 3) tri++; });
    if (tri || dbl >= 2) return Math.PI;
    if (a.arom || dbl === 1) return 120 * Math.PI / 180;
    if (n <= 2 && (a.el === "O" || a.el === "S")) return 104.5 * Math.PI / 180;
    if (a.el === "N" && n === 3) return 107.5 * Math.PI / 180;
    return 109.47 * Math.PI / 180;
  }
  function planarG(g, i) {
    var a = g.atoms[i];
    if (a.arom) return true;
    return a.nbrs.some(function (x) { var b = g.bonds[x.bi]; return !b.arom && b.order === 2; });
  }

  // 1-2, 1-3 쌍은 밀어내기에서 빼요 (결합·각 항이 이미 잡고 있어서)
  function exclusions(g) {
    var ex = {};
    function key(i, j) { return i < j ? i + "," + j : j + "," + i; }
    g.bonds.forEach(function (b) { ex[key(b.a, b.b)] = 1; });
    g.atoms.forEach(function (a, i) {
      for (var p = 0; p < a.nbrs.length; p++)
        for (var q = p + 1; q < a.nbrs.length; q++) ex[key(a.nbrs[p].atom, a.nbrs[q].atom)] = 1;
    });
    return ex;
  }

  function relax(g, ff, ex, steps, pin) {
    var n = g.atoms.length, X = new Float64Array(n * 3), G = new Float64Array(n * 3), V = new Float64Array(n * 3);
    for (var i = 0; i < n; i++) { X[i * 3] = g.atoms[i].x; X[i * 3 + 1] = g.atoms[i].y; X[i * 3 + 2] = g.atoms[i].z; }
    var KB = 12, KA = 3.0, KP = 14, KT = 0, KR = 2.5;
    var step = 0.012;

    function sub(a, b) { return [X[a * 3] - X[b * 3], X[a * 3 + 1] - X[b * 3 + 1], X[a * 3 + 2] - X[b * 3 + 2]]; }
    function norm(v) { return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) || 1e-9; }
    function cross(u, v) { return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]; }
    function dot(u, v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; }
    function add(gi, v, s) { G[gi * 3] += v[0] * s; G[gi * 3 + 1] += v[1] * s; G[gi * 3 + 2] += v[2] * s; }

    for (var it = 0; it < steps; it++) {
      G.fill(0);
      // 결합 길이
      ff.bonds.forEach(function (b) {
        var d = sub(b.i, b.j), r = norm(d), f = 2 * KB * (r - b.r0) / r;
        add(b.i, d, f); add(b.j, d, -f);
      });
      // 결합각 — 세 원자 사이 거리로 바꿔서 간단히
      ff.angles.forEach(function (a) {
        var ri = sub(a.i, a.j), rk = sub(a.k, a.j);
        var li = norm(ri), lk = norm(rk);
        var target = Math.sqrt(li * li + lk * lk - 2 * li * lk * Math.cos(a.t0));
        var d = sub(a.i, a.k), r = norm(d), f = 2 * KA * (r - target) / r;
        add(a.i, d, f); add(a.k, d, -f);
      });
      // sp2 평면성 — 가운데 원자를 세 이웃이 만드는 평면으로
      ff.planes.forEach(function (p) {
        var u = sub(p.b, p.a), v = sub(p.d, p.a), nvec = cross(u, v), ln = norm(nvec);
        nvec = [nvec[0] / ln, nvec[1] / ln, nvec[2] / ln];
        var w = sub(p.c, p.a), h = dot(w, nvec);
        add(p.c, nvec, 2 * KP * h);
        add(p.a, nvec, -2 * KP * h / 3); add(p.b, nvec, -2 * KP * h / 3); add(p.d, nvec, -2 * KP * h / 3);
      });
      // 이중결합·방향족 평면 유지 — 네 원자가 한 평면에 오도록
      ff.torsions.forEach(function (t) {
        var u = sub(t.j, t.i), v = sub(t.k, t.j), w = sub(t.l, t.k);
        var n1 = cross(u, v), n2 = cross(v, w), l1 = norm(n1), l2 = norm(n2);
        var c = dot(n1, n2) / (l1 * l2);
        var s2 = 1 - c * c;                       // sin²φ — 0 또는 180도에서 0
        if (s2 < 1e-6) return;
        var g1 = 2 * KT * (c >= 0 ? -1 : 1) * Math.sqrt(s2);
        var d1 = [n1[0] / l1, n1[1] / l1, n1[2] / l1], d2 = [n2[0] / l2, n2[1] / l2, n2[2] / l2];
        add(t.i, d1, g1); add(t.l, d2, -g1);
      });
      // 떨어진 원자끼리 밀어내기
      for (var a1 = 0; a1 < n; a1++) {
        for (var b1 = a1 + 1; b1 < n; b1++) {
          var k = a1 + "," + b1; if (ex[k]) continue;
          var d2v = sub(a1, b1), r2 = norm(d2v);
          var r0 = (rvdw(g.atoms[a1].el) + rvdw(g.atoms[b1].el)) * 0.80;
          if (r2 >= r0) continue;
          var f2 = 2 * KR * (r2 - r0) / r2;
          add(a1, d2v, f2); add(b1, d2v, -f2);
        }
      }
      // 내리막으로 한 걸음
      var mx = 0;
      for (var q = 0; q < n * 3; q++) { var v2 = Math.abs(G[q]); if (v2 > mx) mx = v2; }
      var s = mx > 1e-9 ? Math.min(step, 0.35 / mx) : step;
      var mom = it < steps * 0.85 ? 0.88 : 0;            // 끝에서는 관성을 꺼서 흔들리지 않게
      for (var q2 = 0; q2 < n * 3; q2++) { V[q2] = mom * V[q2] - G[q2] * s; X[q2] += V[q2]; }
      if (pin) Object.keys(pin).forEach(function (k) {   // 쐐기로 정해진 앞뒤는 지켜요
        var i = +k, z = X[i * 3 + 2];
        if (z * pin[i] < 0.08) { X[i * 3 + 2] = pin[i] * 0.3; V[i * 3 + 2] = 0; }
      });
    }
    for (var i2 = 0; i2 < n; i2++) { g.atoms[i2].x = X[i2 * 3]; g.atoms[i2].y = X[i2 * 3 + 1]; g.atoms[i2].z = X[i2 * 3 + 2]; }
  }

  // 붙어 있는 방향족 고리는 한 평면에 — 눈에 가장 거슬리는 부분이라 마지막에 맞춰줍니다
  function flattenAromatics(g, m) {
    var rings = (m.allRings || m.rings || []).filter(function (r) { return r.every(function (i) { return m.atoms[i].arom; }); });
    if (!rings.length) return;
    var sys = [];                                        // 공유 원자가 있으면 한 덩어리로
    rings.forEach(function (r) {
      var hit = sys.filter(function (S) { return r.some(function (i) { return S.indexOf(i) >= 0; }); });
      if (!hit.length) { sys.push(r.slice()); return; }
      var merged = hit.reduce(function (acc, S) { S.forEach(function (i) { if (acc.indexOf(i) < 0) acc.push(i); }); return acc; }, r.slice());
      sys = sys.filter(function (S) { return hit.indexOf(S) < 0; });
      sys.push(merged);
    });
    sys.forEach(function (S) {
      var pts = S.map(function (i) { return g.atoms[i]; });
      // 고리에 바로 붙은 원자도 같은 평면에 둡니다 (sp2 치환기)
      var ext = S.slice();
      S.forEach(function (i) {
        g.atoms[i].nbrs.forEach(function (nb) {
          if (ext.indexOf(nb.atom) < 0 && g.atoms[nb.atom].el === "H") ext.push(nb.atom);
        });
      });
      var c = [0, 0, 0];
      pts.forEach(function (a) { c[0] += a.x; c[1] += a.y; c[2] += a.z; });
      c = c.map(function (v) { return v / pts.length; });
      var Sm = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      pts.forEach(function (a) {
        var d = [a.x - c[0], a.y - c[1], a.z - c[2]];
        for (var i = 0; i < 3; i++) for (var j = 0; j < 3; j++) Sm[i][j] += d[i] * d[j];
      });
      var tr = Sm[0][0] + Sm[1][1] + Sm[2][2], v = [0.3, 0.5, 0.81];
      for (var it = 0; it < 80; it++) {
        var w = [0, 0, 0];
        for (var i2 = 0; i2 < 3; i2++) { w[i2] = tr * v[i2]; for (var j2 = 0; j2 < 3; j2++) w[i2] -= Sm[i2][j2] * v[j2]; }
        var l = Math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) || 1;
        v = [w[0] / l, w[1] / l, w[2] / l];
      }
      ext.forEach(function (i) {
        var a = g.atoms[i], d = [a.x - c[0], a.y - c[1], a.z - c[2]];
        var h = d[0] * v[0] + d[1] * v[1] + d[2] * v[2];
        a.x -= v[0] * h; a.y -= v[1] * h; a.z -= v[2] * h;
      });
    });
  }

  // 무게중심을 원점으로, 가장 퍼진 방향을 화면에 맞추기
  function center(g) {
    var n = g.atoms.length, cx = 0, cy = 0, cz = 0;
    g.atoms.forEach(function (a) { cx += a.x; cy += a.y; cz += a.z; });
    cx /= n; cy /= n; cz /= n;
    g.atoms.forEach(function (a) { a.x -= cx; a.y -= cy; a.z -= cz; });
  }

  // 평면에서 얼마나 벗어났는지 (0에 가까우면 아직 납작한 것)
  function flatness(g) {
    var pts = g.atoms.filter(function (a) { return a.heavy; });
    if (pts.length < 4) return 1;
    var n = pts.length, c = [0, 0, 0];
    pts.forEach(function (a) { c[0] += a.x; c[1] += a.y; c[2] += a.z; });
    c = c.map(function (v) { return v / n; });
    var S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    pts.forEach(function (a) {
      var d = [a.x - c[0], a.y - c[1], a.z - c[2]];
      for (var i = 0; i < 3; i++) for (var j = 0; j < 3; j++) S[i][j] += d[i] * d[j];
    });
    // 멱반복으로 가장 작은 분산 방향의 크기를 어림
    var tr = S[0][0] + S[1][1] + S[2][2];
    var v = [0.3, 0.5, 0.81], it;
    for (it = 0; it < 60; it++) {                      // (tr·I − S) 의 최대 고유벡터 = S 의 최소
      var w = [0, 0, 0];
      for (var i = 0; i < 3; i++) { w[i] = tr * v[i]; for (var j = 0; j < 3; j++) w[i] -= S[i][j] * v[j]; }
      var l = Math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) || 1;
      v = [w[0] / l, w[1] / l, w[2] / l];
    }
    var out = 0;
    for (var i2 = 0; i2 < 3; i2++) for (var j2 = 0; j2 < 3; j2++) out += v[i2] * S[i2][j2] * v[j2];
    return Math.sqrt(Math.max(0, out) / n);            // 평면에서 벗어난 평균 거리 (Å)
  }

  function embed(m, opts) {
    opts = opts || {};
    var g = withHydrogens(m);
    seed(g, m);
    var ff = build(g, m), ex = exclusions(g);
    relax(g, ff, ex, opts.steps || Math.max(700, Math.min(2600, g.atoms.length * 22)), g.pin);
    flattenAromatics(g, m);
    relax(g, ff, ex, 120, g.pin);                        // 눌러준 뒤 짧게 다시 풀기
    center(g);
    g.flat = flatness(g);
    return g;
  }

  Chem.embed3D = embed;
  Chem.mol3dFlatness = flatness;
})(typeof window !== "undefined" ? window : globalThis);
