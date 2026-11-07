/* 의약화학 오피스 — 입체이성질체 만들기
 *
 * 그린 구조에서 비칭 탄소와 E/Z 이중결합을 찾아, 가능한 조합마다 분자파일을 하나씩 만듭니다.
 * 비칭 탄소는 쐐기(위)/점선(아래)을 바꿔 끼우고, 이중결합은 한쪽 치환기를
 * 이중결합 축을 기준으로 뒤집어서 E↔Z 를 바꿔요.
 * 만들어진 분자파일을 편집기(Indigo)에 넘기면 SMILES 와 R/S·E/Z 표시가 나옵니다.
 */
(function (root) {
  "use strict";
  var Chem = root.Chem;

  // ---------- V2000 분자파일 쓰기 ----------
  function writeMolfile(mol, title) {
    var L = [title || "", "  MedChemOffice", "", ""];
    L[3] = pad(mol.atoms.length, 3) + pad(mol.bonds.length, 3) + "  0  0  1  0  0  0  0  0999 V2000";   // 5번째 = chiral 플래그 (절대 입체배치)
    mol.atoms.forEach(function (a) {
      L.push(f10(a.x || 0) + f10(a.y || 0) + f10(0) + " " + (a.el + "  ").slice(0, 3) +
        " 0" + pad(chgCode(a.charge || 0), 3) + "  0  0  0  0  0  0  0  0  0  0");
    });
    mol.bonds.forEach(function (b) {
      var i1 = b.from === b.b ? b.b : b.a, i2 = i1 === b.a ? b.b : b.a;   // 쐐기의 좁은 쪽이 먼저
      L.push(pad(i1 + 1, 3) + pad(i2 + 1, 3) + pad(b.arom ? 4 : b.order, 3) + pad(b.stereo || 0, 3) + "  0  0  0");
    });
    var chg = [];
    mol.atoms.forEach(function (a, i) { if (a.charge) chg.push([i + 1, a.charge]); });
    while (chg.length) {
      var part = chg.splice(0, 8);
      L.push("M  CHG" + pad(part.length, 3) + part.map(function (c) { return pad(c[0], 4) + pad(c[1], 4); }).join(""));
    }
    L.push("M  END");
    return L.join("\n");
  }
  function pad(v, n) { return String(v).padStart(n); }
  function f10(v) { return v.toFixed(4).padStart(10); }
  function chgCode(c) { return { 3: 1, 2: 2, 1: 3, "-1": 5, "-2": 6, "-3": 7 }[c] || 0; }

  // ---------- 얕은 복사 (좌표·결합 그대로) ----------
  function clone(mol) {
    var m = { atoms: [], bonds: [] };
    mol.atoms.forEach(function (a) {
      m.atoms.push({ el: a.el, charge: a.charge || 0, arom: !!a.arom, h: a.h, x: a.x || 0, y: a.y || 0, nbrs: [] });
    });
    mol.bonds.forEach(function (b) {
      m.bonds.push({ a: b.a, b: b.b, order: b.order, arom: !!b.arom, stereo: 0, from: b.a, inRing: b.inRing });
    });
    m.bonds.forEach(function (b) {
      m.atoms[b.a].nbrs.push({ atom: b.b, bond: b });
      m.atoms[b.b].nbrs.push({ atom: b.a, bond: b });
    });
    return m;
  }

  // ---------- 가지 모으기 ----------
  function branch(m, start, blocked) {
    var out = [start], seen = {}; seen[start] = 1; seen[blocked] = 1;
    var q = [start];
    while (q.length) {
      var u = q.shift();
      m.atoms[u].nbrs.forEach(function (n) {
        if (seen[n.atom]) return;
        seen[n.atom] = 1; out.push(n.atom); q.push(n.atom);
      });
    }
    return out;
  }
  // 이중결합 축을 기준으로 뒤집기 — E ↔ Z
  function flipAcross(m, pa, pb, atoms) {
    var ax = m.atoms[pa].x, ay = m.atoms[pa].y;
    var dx = m.atoms[pb].x - ax, dy = m.atoms[pb].y - ay;
    var len = Math.hypot(dx, dy) || 1; dx /= len; dy /= len;
    atoms.forEach(function (i) {
      var a = m.atoms[i], vx = a.x - ax, vy = a.y - ay;
      var t = vx * dx + vy * dy;                      // 축 방향 성분은 그대로
      var px = t * dx, py = t * dy;                   // 수직 성분만 부호를 뒤집어요
      a.x = ax + px - (vx - px);
      a.y = ay + py - (vy - py);
    });
  }

  /* 가능한 조합을 전부 만들어요.
   * pot = { centers: [원자번호…], bonds: [결합번호…] }
   * 반환: [{ molfile, bits }] — bits 는 각 자리의 선택(0/1)
   */
  function variants(mol, pot, limit) {
    var slots = pot.centers.length + pot.bonds.length;
    var total = Math.pow(2, slots);
    limit = limit || 32;
    var n = Math.min(total, limit), out = [];
    for (var k = 0; k < n; k++) {
      var m = clone(mol), bits = [];
      pot.centers.forEach(function (ai, s) {
        var up = (k >> s) & 1; bits.push(up);
        var a = m.atoms[ai];
        // 이 중심에 걸린 기존 쐐기는 지우고, 첫 이웃에 새로 하나만 겁니다
        a.nbrs.forEach(function (nb) { if (nb.bond.a === ai || nb.bond.b === ai) { nb.bond.stereo = 0; } });
        var pick = a.nbrs.filter(function (nb) { return !nb.bond.arom && nb.bond.order === 1; })[0] || a.nbrs[0];
        if (!pick) return;
        pick.bond.from = ai;                           // 좁은 쪽이 입체중심
        pick.bond.stereo = up ? 1 : 6;                 // 1 = 쐐기(위), 6 = 점선(아래)
      });
      pot.bonds.forEach(function (bi, s2) {
        var flip = (k >> (pot.centers.length + s2)) & 1; bits.push(flip);
        if (!flip) return;
        var b = m.bonds[bi];
        // 한쪽 끝에 달린 치환기를 전부 축 기준으로 뒤집어요
        var subs = m.atoms[b.b].nbrs.filter(function (nb) { return nb.atom !== b.a; });
        var move = [];
        subs.forEach(function (nb) { branch(m, nb.atom, b.b).forEach(function (i) { if (move.indexOf(i) < 0) move.push(i); }); });
        if (move.length) flipAcross(m, b.a, b.b, move);
      });
      out.push({ molfile: writeMolfile(m, ""), bits: bits });
    }
    return { list: out, total: total, shown: n };
  }

  Chem.writeMolfile = writeMolfile;
  Chem.stereoVariants = variants;
})(typeof window !== "undefined" ? window : globalThis);
