/* 의약화학 도구 — 화합물 목록 공유 저장소
 * 여섯 도구가 모두 이 파일을 통해 같은 화합물 목록을 읽고 씁니다.
 * list 배열은 한 번 만들어지면 교체되지 않아요. 그래야 각 도구가 var list = Store.list 로
 * 붙잡아 둔 참조가 계속 살아 있습니다. 지우기도 splice·length=0 으로 제자리에서 합니다.
 */
(function (root) {
  "use strict";
  var KEY = "medchem-office-compounds-v1";
  var DEF_TGT = "활성";
  var list = [];

  function read() {
    var v = [];
    try { v = JSON.parse(localStorage.getItem(KEY) || "[]"); } catch (e) { v = []; }
    return Array.isArray(v) ? v : [];
  }
  function fill(src) { list.length = 0; src.forEach(function (c) { list.push(c); }); }
  fill(read());

  var onError = null;                                  // 저장 실패를 알릴 곳 (도구에서 지정)
  var onExternal = null;                               // 다른 도구가 목록을 바꿨을 때

  function save() {
    try {
      localStorage.setItem(KEY, JSON.stringify(list.map(function (c) {
        return { id: c.id, code: c.code, smiles: c.smiles, act: c.act, unit: c.unit, acts: actsOf(c),
          memo: c.memo, created: c.created, updated: c.updated || null,
          pka: c.pka == null ? null : c.pka, logd: c.logd == null ? null : c.logd,
          nmr: c.nmr || null, dock: c.dock || null, props: c.props || null, tox: c.tox || null };
      })));
      return true;
    } catch (e) {
      if (onError) onError("저장 공간이 부족해요. CSV로 내보내 두세요.");
      return false;
    }
  }
  function remove(c) { var i = list.indexOf(c); if (i >= 0) list.splice(i, 1); }
  function clear() { list.length = 0; }
  function reload() { fill(read()); keyCache = {}; }

  // 다른 도구(다른 탭)에서 목록이 바뀌면 따라 갱신해요
  root.addEventListener("storage", function (e) {
    if (e.key !== KEY) return;
    reload();
    if (onExternal) onExternal();
  });

  // ---------- 같은 물질 가려내기 ----------
  var keyCache = {};
  function keyOf(smi) {
    if (!smi) return "";
    if (keyCache[smi] == null) { var m = root.Chem && Chem.fromAny(smi); keyCache[smi] = m ? Chem.structKey(m) : ""; }
    return keyCache[smi];
  }
  function findSame(smi, code) {           // 코드를 적었으면 구조+코드, 비웠으면 구조만
    var k = keyOf(smi); if (!k) return -1;
    var hit = -1;
    list.forEach(function (c, i) {
      if (keyOf(c.smiles) !== k) return;
      if (code && (c.code || "").trim() !== code) return;
      if (hit < 0 || (c.created || 0) >= (list[hit].created || 0)) hit = i;
    });
    return hit;
  }
  function codeClash(smi, code) {          // 같은 코드인데 구조가 다른 행
    if (!code) return -1;
    var k = keyOf(smi);
    for (var i = 0; i < list.length; i++) if ((list[i].code || "").trim() === code && keyOf(list[i].smiles) !== k) return i;
    return -1;
  }
  function byId(id) { for (var i = 0; i < list.length; i++) if (list[i].id === id) return list[i]; return null; }

  // ---------- 타깃별 활성값 ----------
  function actsOf(c) {
    if (!c.acts || typeof c.acts !== "object") {
      c.acts = {};
      if (c.act != null && !isNaN(+c.act)) c.acts[DEF_TGT] = { v: +c.act, u: c.unit || "nM" };
    }
    return c.acts;
  }
  function targets() {
    var seen = [];
    list.slice().sort(function (a, b) { return (a.created || 0) - (b.created || 0); })
      .forEach(function (c) { Object.keys(actsOf(c)).forEach(function (t) { if (seen.indexOf(t) < 0) seen.push(t); }); });
    return seen;
  }
  function tgtUnit(t) {                    // 열 전체가 같은 단위면 머리글에 씁니다
    var u = null;
    for (var i = 0; i < list.length; i++) { var a = actsOf(list[i])[t]; if (!a) continue; if (u == null) u = a.u; else if (u !== a.u) return null; }
    return u;
  }
  function renameTarget(from, to) {        // 모든 화합물에서 타깃 이름 바꾸기
    var kept = 0;
    list.forEach(function (c) {
      var A = actsOf(c), v = A[from];
      if (v == null) return;
      if (A[to] != null) { kept++; delete A[from]; return; }   // 받는 열에 값이 있으면 그쪽을 남겨요
      A[to] = v; delete A[from];
    });
    return kept;
  }
  function mergeInto(c, v) {
    var changed = [], tn = (v.tgt || "").trim(), A = actsOf(c);
    if (v.act != null && v.act !== "") {
      var t = tn || DEF_TGT;
      A[t] = { v: +v.act, u: v.unit || "nM" };
      c.act = null; c.unit = null;
      changed.push((t === DEF_TGT ? "활성값" : t + " 활성값") + " 넣었어요");
    } else if (tn && !A[tn] && A[DEF_TGT]) {
      A[tn] = A[DEF_TGT]; delete A[DEF_TGT];               // 값 먼저, 이름 나중에 넣은 경우
      changed.push("이름 없던 활성값을 " + tn + " 열로 옮겼어요");
    }
    if (v.memo) {
      var old = (c.memo || "").trim();
      if (!old) { c.memo = v.memo; changed.push("메모 넣었어요"); }
      else if (old.indexOf(v.memo) < 0) { c.memo = old + " · " + v.memo; changed.push("메모 덧붙였어요"); }
    }
    if (v.pka != null) c.pka = v.pka;
    if (v.logd != null) c.logd = v.logd;
    c.updated = Date.now();
    return changed;
  }
  function autoName(k) {                   // 이름을 비우면 화합물-1, 화합물-2 …
    var used = {}; list.forEach(function (c) { var m = /^화합물-(\d+)$/.exec(c.code || ""); if (m) used[+m[1]] = 1; });
    var n = 1, made = [];
    for (var i = 0; i <= (k || 0); i++) { while (used[n]) n++; used[n] = 1; made.push("화합물-" + n); }
    return k == null ? made[0] : made[k];
  }
  function fresh(code, smiles) {
    return { id: "c" + Date.now().toString(36) + Math.floor(Math.random() * 1e4).toString(36),
      code: code, smiles: smiles, act: null, unit: null, acts: {}, memo: "", created: Date.now(),
      pka: null, logd: null };
  }

  root.Store = {
    KEY: KEY, DEF_TGT: DEF_TGT, list: list,
    save: save, remove: remove, clear: clear, reload: reload, fresh: fresh, byId: byId,
    keyOf: keyOf, findSame: findSame, codeClash: codeClash,
    actsOf: actsOf, targets: targets, tgtUnit: tgtUnit, renameTarget: renameTarget,
    mergeInto: mergeInto, autoName: autoName,
    set onError(f) { onError = f; }, set onExternal(f) { onExternal = f; }
  };
})(window);
