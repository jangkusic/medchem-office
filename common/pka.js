/* 의약화학 오피스 — pKa 예상값과 cLogD(pH 7.4) 계산
 * 작용기별 기준 pKa 표에 주변 원자의 전자 끄는 효과를 더하는 규칙 기반 추정이에요.
 * 예측 모델이 아니라 의약화학에서 쓰는 경험 규칙이라, 보통 ±1 정도 차이가 나고
 * 고리가 많이 붙은 방향족 질소는 더 벌어질 수 있어요.
 */
(function (root) {
  "use strict";
  var Chem = root.Chem;

  // 전자 끄는 정도 (거리 3에서의 기준값) — 값이 클수록 염기성을 많이 떨어뜨림
  var W = { Ocarbonyl: 1.6, O: 1.25, Nnitrile: 1.6, N: 0.55, F: 1.7, Cl: 1.1, Br: 1.0, I: 0.7, S: 0.4, P: 0.3 };
  var DIST = { 2: 1.35, 3: 1.0, 4: 0.33, 5: 0.12 };

  function bfs(m, from, max) { // 원자별 결합 거리와 최단 경로 개수 (고리에서는 양쪽으로 가므로 2개가 될 수 있어요)
    var d = {}, c = {}, q = [from]; d[from] = 0; c[from] = 1;
    while (q.length) {
      var u = q.shift(); if (d[u] >= max) continue;
      m.atoms[u].nbrs.forEach(function (n) {
        if (d[n.atom] == null) { d[n.atom] = d[u] + 1; c[n.atom] = c[u]; q.push(n.atom); }
        else if (d[n.atom] === d[u] + 1) c[n.atom] += c[u];
      });
    }
    return { d: d, n: c };
  }
  function isCarbonylO(m, i) {
    var a = m.atoms[i]; if (a.el !== "O" || a.nbrs.length !== 1) return false;
    var b = a.nbrs[0].bond; return !b.arom && b.order === 2;
  }
  function sulfonylO(m, i) {
    var a = m.atoms[i]; if (a.el !== "O") return false;
    return a.nbrs.some(function (n) { return m.atoms[n.atom].el === "S" || m.atoms[n.atom].el === "N"; }) && a.nbrs.length === 1;
  }
  function isAmideN(m, i) {
    var a = m.atoms[i];
    return a.nbrs.some(function (n) {
      var t = m.atoms[n.atom]; if (n.bond.arom || n.bond.order !== 1) return false;
      if (t.el === "C" && t.nbrs.some(function (q) { return !q.bond.arom && q.bond.order === 2 && ["O", "S"].indexOf(m.atoms[q.atom].el) >= 0; })) return true;
      if (t.el === "S" && t.nbrs.filter(function (q) { return !q.bond.arom && q.bond.order === 2 && m.atoms[q.atom].el === "O"; }).length >= 1) return true;
      if (t.el === "N" && t.nbrs.some(function (q) { return m.atoms[q.atom].charge === -1 && m.atoms[q.atom].el === "O"; })) return true; // 나이트로
      return false;
    });
  }
  function arylNbr(m, i) { return m.atoms[i].nbrs.some(function (n) { return m.atoms[n.atom].arom; }); }

  // 지방족 아민 주변의 전자 끄는 효과 합
  function pull(m, from, noAryl) {
    var B = bfs(m, from, 5), d = B.d, total = 0, notes = {};
    Object.keys(d).forEach(function (k) {
      var i = +k, dist = d[i], a = m.atoms[i];
      if (!dist || dist < 2) return;
      var fac = DIST[dist]; if (!fac) return;
      fac *= Math.min(2, B.n[i] || 1); // 고리 양쪽으로 영향이 오면 두 배
      var w = 0, nm = "";
      if (a.el === "O") { if (isCarbonylO(m, i) || sulfonylO(m, i)) { w = W.Ocarbonyl; nm = "카보닐·설포닐"; } else { w = W.O; nm = "산소"; } }
      else if (a.el === "N") {
        var trip = a.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 3; });
        w = trip ? W.Nnitrile : W.N; nm = trip ? "나이트릴" : "질소";
      } else if (W[a.el]) { w = W[a.el]; nm = a.el; }
      else if (a.el === "C" && a.arom && dist === 2) { if (noAryl) return; w = 1.2 / DIST[2]; nm = "방향족 고리"; }
      if (!w) return;
      var v = w * fac; total += v;
      notes[nm] = (notes[nm] || 0) + v;
    });
    total = Math.min(total, 7);
    var txt = Object.keys(notes).sort(function (a, b) { return notes[b] - notes[a]; }).slice(0, 3)
      .map(function (k) { return k + " −" + notes[k].toFixed(1); }).join(", ");
    return { v: total, txt: txt };
  }
  // 고리계(붙어 있는 고리까지)에 있는 전자 끄는 기 개수
  function fused(m, ring) {
    var rs = (m.allRings || m.rings || []), set = ring.slice(), grew = true;
    while (grew) {
      grew = false;
      rs.forEach(function (r) {
        if (r.some(function (i) { return set.indexOf(i) >= 0; }) && r.some(function (i) { return set.indexOf(i) < 0; })) {
          r.forEach(function (i) { if (set.indexOf(i) < 0) set.push(i); }); grew = true;
        }
      });
    }
    return set;
  }
  function ewgLike(m, t) {
    if (t.el === "C" && t.nbrs.some(function (q) { return !q.bond.arom && (q.bond.order === 2 && ["O", "S"].indexOf(m.atoms[q.atom].el) >= 0 || q.bond.order === 3); })) return true;
    if (t.el === "N" && t.charge === 1) return true;
    if (t.el === "S" && t.nbrs.filter(function (q) { return m.atoms[q.atom].el === "O"; }).length >= 2) return true;
    if (t.el === "C" && t.nbrs.filter(function (q) { return m.atoms[q.atom].el === "F"; }).length === 3) return true;
    return false;
  }
  function ringEwg(m, ring) {
    var sys = fused(m, ring), n = 0;
    sys.forEach(function (i) {
      if (ewgLike(m, m.atoms[i])) { n++; return; }                      // 고리 안에 있는 카보닐 (예: 카페인)
      m.atoms[i].nbrs.forEach(function (nb) {
        if (sys.indexOf(nb.atom) >= 0) return;
        if (ewgLike(m, m.atoms[nb.atom])) n++;                          // 고리에 붙은 전자 끄는 기
      });
    });
    return n;
  }
  function ringOf(m, i) {
    var rs = (m.allRings || m.rings || []).filter(function (r) { return r.indexOf(i) >= 0; });
    rs.sort(function (a, b) { return a.length - b.length; });
    return rs[0] || null;
  }

  function estimate(mol) {
    var m = mol, bases = [], acids = [];
    m.atoms.forEach(function (a, i) {
      // ---------- 염기 ----------
      if (a.el === "N" && a.charge <= 1 && !isAmideN(m, i)) {
        if (!a.arom) {
          if (a.nbrs.length === 4 && a.charge === 1) return;                 // 4차 암모늄: pKa 없음
          if (a.nbrs.some(function (n) { return m.atoms[n.atom].el === "O" && m.atoms[n.atom].charge === -1; })) return; // N-옥사이드·나이트로
          // 아미딘·구아니딘
          var amid = null;
          a.nbrs.forEach(function (n) {
            var c = m.atoms[n.atom];
            if (c.el !== "C" || c.arom) return;
            var dN = c.nbrs.filter(function (q) { return !q.bond.arom && q.bond.order === 2 && m.atoms[q.atom].el === "N"; });
            var sN = c.nbrs.filter(function (q) { return !q.bond.arom && q.bond.order === 1 && m.atoms[q.atom].el === "N"; });
            if (dN.length === 1 && sN.length >= 1) amid = sN.length >= 2 ? "guan" : "amid";
          });
          if (amid) { bases.push({ atom: i, pka: amid === "guan" ? 13.5 : 11.5, label: amid === "guan" ? "구아니딘" : "아미딘", basis: (amid === "guan" ? "구아니딘 13.5" : "아미딘 11.5") }); return; }
          if (arylNbr(m, i)) { // 아닐린
            var p0 = a.h >= 1 ? 4.6 : 5.1, c0 = pull(m, i, true);
            bases.push({ atom: i, pka: Math.max(-2, p0 - c0.v), label: "아닐린", basis: "아닐린 " + p0 + (c0.txt ? ", " + c0.txt : "") });
            return;
          }
          var ring3 = a.ringSizes.indexOf(3) >= 0, ring4 = a.ringSizes.indexOf(4) >= 0;
          var base = a.h >= 2 ? 10.6 : a.h === 1 ? 10.9 : 9.9;
          if (ring3) base = 8.0; else if (ring4) base = 11.3;
          var c = pull(m, i);
          bases.push({ atom: i, pka: Math.max(-2, base - c.v), label: a.h >= 2 ? "1차 아민" : a.h === 1 ? "2차 아민" : "3차 아민",
            basis: (a.h >= 2 ? "1차 아민 " : a.h === 1 ? "2차 아민 " : "3차 아민 ") + base.toFixed(1) + (c.txt ? ", " + c.txt : "") });
        } else if (a.h === 0 && a.nbrs.length === 2 && a.charge === 0) {
          // 방향족 질소 (피리딘·이미다졸형)
          var r = ringOf(m, i); if (!r) return;
          var ns = r.filter(function (x) { return m.atoms[x].el === "N"; }).length;
          var hasSO = r.some(function (x) { return m.atoms[x].el === "S" || m.atoms[x].el === "O"; });
          var adjN = a.nbrs.some(function (n) { return m.atoms[n.atom].el === "N" && r.indexOf(n.atom) >= 0; });
          var b2, nm2;
          if (r.length === 6) { b2 = ns >= 3 ? 0.6 : ns === 2 ? 1.3 : 5.2; nm2 = ns >= 2 ? "다이아진" : "피리딘"; }
          else if (r.length === 5) {
            if (hasSO) { b2 = r.some(function (x) { return m.atoms[x].el === "O"; }) ? 0.8 : 2.5; nm2 = r.some(function (x) { return m.atoms[x].el === "O"; }) ? "옥사졸" : "티아졸"; }
            else if (ns >= 3) { b2 = 2.3; nm2 = "트라이아졸"; }
            else if (adjN) { b2 = 2.5; nm2 = "피라졸"; }
            else { b2 = 7.0; nm2 = "이미다졸"; }
          } else return;
          var ew = ringEwg(m, r), pk = b2 - ew * 1.8;
          bases.push({ atom: i, pka: Math.max(-3, pk), label: nm2, basis: nm2 + " " + b2.toFixed(1) + (ew ? ", 전자 끄는 기 " + ew + "개 −" + (ew * 1.8).toFixed(1) : ""), weak: true });
        }
      }
      // ---------- 산 ----------
      if (a.el === "C" && !a.arom) { // 카복실산
        var dO = a.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 2 && m.atoms[n.atom].el === "O"; });
        var oh = a.nbrs.filter(function (n) { var t = m.atoms[n.atom]; return t.el === "O" && n.bond.order === 1 && t.nbrs.length === 1 && (t.h === 1 || t.charge === -1); });
        if (dO.length && oh.length) {
          var ar = arylNbr(m, i), p1 = ar ? 4.2 : 4.6, cc = pull(m, i, true);
          acids.push({ atom: i, pka: Math.max(-1, p1 - cc.v * 0.45), label: "카복실산",
            basis: "카복실산 " + p1 + (cc.txt ? ", " + cc.txt : "") });
        }
      }
      if (a.el === "O" && a.h === 1 && arylNbr(m, i)) { // 페놀
        var r2 = ringOf(m, a.nbrs[0].atom), ew2 = r2 ? ringEwg(m, r2) : 0;
        acids.push({ atom: i, pka: Math.max(0, 10.0 - ew2 * 2.2), label: "페놀", basis: "페놀 10.0" + (ew2 ? ", 전자 끄는 기 " + ew2 + "개" : "") });
      }
      if (a.el === "N" && !a.arom && a.h >= 1) { // 설폰아마이드·아실설폰아마이드
        var so2 = a.nbrs.some(function (n) { var t = m.atoms[n.atom]; return t.el === "S" && t.nbrs.filter(function (q) { return !q.bond.arom && q.bond.order === 2 && m.atoms[q.atom].el === "O"; }).length === 2; });
        var co = a.nbrs.some(function (n) { var t = m.atoms[n.atom]; return t.el === "C" && t.nbrs.some(function (q) { return !q.bond.arom && q.bond.order === 2 && m.atoms[q.atom].el === "O"; }); });
        if (so2 && co) acids.push({ atom: i, pka: 4.5, label: "아실설폰아마이드", basis: "아실설폰아마이드 4.5" });
        else if (so2) acids.push({ atom: i, pka: 10.1, label: "설폰아마이드", basis: "설폰아마이드 10.1" });
      }
      if (a.el === "S" && a.h === 1) acids.push({ atom: i, pka: 8.5, label: "싸이올", basis: "싸이올 8.5" });
    });
    // 테트라졸
    (m.rings || []).forEach(function (r) {
      if (r.length !== 5) return;
      if (r.filter(function (i) { return m.atoms[i].el === "N"; }).length !== 4) return;
      if (!r.some(function (i) { return m.atoms[i].el === "N" && m.atoms[i].h === 1; })) return;
      acids.push({ atom: r[0], pka: 4.9, label: "테트라졸", basis: "테트라졸 4.9" });
    });

    bases.sort(function (x, y) { return y.pka - x.pka; });
    acids.sort(function (x, y) { return x.pka - y.pka; });
    return { bases: bases, acids: acids, base: bases[0] || null, acid: acids[0] || null };
  }

  // pH 7.4에서의 분배계수 (이온화된 만큼 친유성이 떨어짐)
  function logD(logp, pkaBase, pkaAcid, pH) {
    pH = pH == null ? 7.4 : pH;
    var d = logp;
    if (pkaBase != null) d -= Math.log10(1 + Math.pow(10, pkaBase - pH));
    if (pkaAcid != null) d -= Math.log10(1 + Math.pow(10, pH - pkaAcid));
    return d;
  }

  Chem.estimatePka = estimate;
  Chem.logD = logD;
})(typeof window !== "undefined" ? window : globalThis);
