/* 의약화학 오피스 — 위험 신호(구조 경고)와 구조 그림
 * - alerts(mol, props): 대사 취약 자리 / CYP / hERG / DDI 구조 경고를 찾아 원자 번호와 함께 돌려줌
 * - depict(mol, opts): MOL 좌표로 2D 구조를 그리고, 지정한 원자를 색으로 강조
 * 주의: 모두 의약화학에서 널리 쓰는 "구조 경고(structural alert)"이지 예측 모델이 아님.
 */
(function (root) {
  "use strict";
  var Chem = root.Chem;

  var CAT = {
    metab: { name: "대사 취약", color: "#c0392b" },
    cyp: { name: "CYP 억제", color: "#c97a09" },
    herg: { name: "hERG", color: "#7b3fa0" },
    ddi: { name: "DDI 요인", color: "#2a6fb5" }
  };

  // ---------- 그래프 도우미 ----------
  function mk(m) {
    var A = m.atoms;
    return {
      A: A,
      nb: function (a) { return a.nbrs.map(function (n) { return { at: A[n.atom], bond: n.bond }; }); },
      heavy: function (a) { return a.nbrs.length; },
      dbl: function (a, els) { return a.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 2 && (!els || els.indexOf(A[n.atom].el) >= 0); }).map(function (n) { return A[n.atom]; }); },
      sgl: function (a, els) { return a.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 1 && (!els || els.indexOf(A[n.atom].el) >= 0); }).map(function (n) { return A[n.atom]; }); },
      trp: function (a) { return a.nbrs.filter(function (n) { return !n.bond.arom && n.bond.order === 3; }).map(function (n) { return A[n.atom]; }); },
      isMe: function (a) { return a.el === "C" && !a.arom && a.h === 3 && a.nbrs.length === 1; },
      isCarbonyl: function (a) { return a.el === "C" && !a.arom && a.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2 && A[n.atom].el === "O"; }); },
      arylNbr: function (a) { return a.nbrs.some(function (n) { return A[n.atom].arom; }); },
      amide: function (a) { // 질소가 카보닐/설포닐에 붙었는지
        return a.nbrs.some(function (n) {
          var t = A[n.atom]; if (n.bond.arom || n.bond.order !== 1) return false;
          if (t.el === "C" && t.nbrs.some(function (q) { return !q.bond.arom && q.bond.order === 2 && ["O", "S", "N"].indexOf(A[q.atom].el) >= 0; })) return true;
          if (t.el === "S" && t.nbrs.filter(function (q) { return !q.bond.arom && q.bond.order === 2 && A[q.atom].el === "O"; }).length >= 1) return true;
          return false;
        });
      }
    };
  }
  // 고리 안의 원자 집합이 조건에 맞는지
  function ringsOfSize(m, n) { return (m.allRings || m.rings || []).filter(function (r) { return r.length === n; }); }
  function aromRing(m, r) { return r.every(function (i) { return m.atoms[i].arom; }); }

  // ---------- 경고 규칙 ----------
  function alerts(m, p) {
    var G = mk(m), A = m.atoms, out = [], seen = {};
    function add(cat, sev, name, why, fix, atoms) {
      atoms = atoms.filter(function (i, k, s) { return s.indexOf(i) === k; }).sort(function (a, b) { return a - b; });
      var key = cat + "|" + name + "|" + atoms.join(",");
      if (seen[key]) return; seen[key] = 1;
      out.push({ cat: cat, sev: sev, name: name, why: why, fix: fix, atoms: atoms });
    }

    // ===== 대사 취약 =====
    A.forEach(function (a, i) {
      // 방향족 메톡시/에톡시 (O-탈알킬화)
      if (a.el === "O" && !a.arom && a.nbrs.length === 2 && a.h === 0) {
        var ns = G.nb(a), ar = ns.filter(function (x) { return x.at.arom; }), al = ns.filter(function (x) { return !x.at.arom && x.at.el === "C" && !G.isCarbonyl(x.at); });
        if (ar.length === 1 && al.length === 1 && al[0].at.h >= 2) {
          var alk = al[0].at, isMe = alk.h === 3;
          add("metab", "high", (isMe ? "방향족 메톡시" : "방향족 알콕시") + " (Ar–O–" + (isMe ? "CH₃" : "CH₂–") + ")",
            "CYP가 가장 흔하게 공격하는 자리예요. O-탈알킬화로 페놀이 생기고, 그 페놀은 2상 대사(글루쿠론산 포합)로 빠르게 빠져나가요.",
            "할로젠(F, Cl)이나 CF₃·OCF₃로 바꾸기 / 메톡시 두 개를 묶어 메틸렌다이옥시 대신 다이옥솔란·옥세탄 고리 검토 / 메틸을 중수소(CD₃)로 치환",
            [i, alk.idx].concat(ar.map(function (x) { return x.at.idx; })));
        }
      }
      // N-알킬 아민 (N-탈알킬화) — 질소 하나당 한 줄로 묶기
      if (a.el === "N" && !a.arom && a.charge <= 0 && !G.amide(a)) {
        var alphas = G.nb(a).filter(function (x) {
          if (x.bond.arom || x.bond.order !== 1) return false;
          var t = x.at; return t.el === "C" && !t.arom && t.h >= 2 && !G.isCarbonyl(t);
        }).map(function (x) { return x.at; });
        if (alphas.length) {
          var me = alphas.filter(function (t) { return t.h === 3 && t.nbrs.length === 1; });
          add("metab", me.length ? "high" : "mid",
            me.length ? "N-메틸 아민" : "N-알킬 아민 (α-C " + alphas.length + "곳)",
            me.length
              ? "아민의 N-메틸은 CYP3A4·2D6가 가장 먼저 떼어내는 자리예요. 1차 대사물이 빨리 생기고 활성이 달라질 수 있어요."
              : "아민 옆 C–H(α 위치)는 N-탈알킬화가 일어나는 자리예요. 고리 아민이면 고리가 열리는 대사도 함께 봐야 해요.",
            "아민을 고리로 묶기(피롤리딘·피페리딘) / 질소 옆에 입체 장애나 F 넣기 / 염기성을 낮춰 아마이드·설폰아마이드로 바꾸기",
            [i].concat(alphas.map(function (t) { return t.idx; })));
        }
      }
      // 벤질 위치
      if (a.el === "C" && !a.arom && a.h >= 1 && G.arylNbr(a) && !G.isCarbonyl(a)) {
        var arn = a.nbrs.filter(function (n) { return A[n.atom].arom; }).map(function (n) { return n.atom; });
        add("metab", a.h >= 2 ? "mid" : "high", "벤질 위치 C–H",
          "방향족 고리 옆 C–H는 결합이 약해서 CYP 수산화가 잘 일어나요. 3차 탄소일수록 더 취약해요.",
          "불소로 치환(CHF, CF₂) / 이 자리를 고리에 포함시키기 / 헤테로 원자를 넣어 전자 밀도 낮추기",
          [i].concat(arn));
      }
      // 알릴 위치
      if (a.el === "C" && !a.arom && a.h >= 1 && a.nbrs.some(function (n) { var t = A[n.atom]; return !t.arom && !n.bond.arom && n.bond.order === 1 && G.dbl(t, ["C"]).length; })) {
        add("metab", "mid", "알릴 위치 C–H", "이중 결합 옆 C–H도 결합이 약해 수산화가 쉬워요. 에폭사이드 생성 경로도 함께 고려해야 해요.", "이중 결합을 고리로 가두거나 포화시키기 / 불소 치환", [i]);
      }
      // t-부틸
      if (a.el === "C" && !a.arom && a.nbrs.length >= 3 && a.nbrs.filter(function (n) { return G.isMe(A[n.atom]); }).length >= 3) {
        add("metab", "high", "t-부틸 (말단 gem-다이메틸)",
          "말단 메틸이 차례로 산화돼 알코올·카복실산까지 가는 대표적인 대사 자리예요.",
          "한쪽 메틸을 F나 OH로 미리 바꾸기 / 사이클로프로필·옥세탄·바이사이클로[1.1.1]펜탄으로 교체",
          [i].concat(a.nbrs.map(function (n) { return n.atom; })));
      }
      // 에스터
      if (G.isCarbonyl(a) && a.el === "C") {
        var os = G.sgl(a, ["O"]).filter(function (t) { return t.nbrs.length === 2 && t.h === 0; });
        if (os.length === 1) {
          var o = os[0], rest = G.nb(o).filter(function (x) { return x.at.idx !== i; })[0];
          if (rest && !G.isCarbonyl(rest.at)) {
            add("metab", "high", "에스터",
              "혈장·간의 에스터레이스가 빠르게 가수분해해요. 경구 약물에서는 보통 몇 분~몇십 분 반감기예요.",
              "아마이드·옥사다이아졸·케톤 등 생물학적 등가체로 교체 / 입체 장애를 주어 가수분해 늦추기 / 의도한 전구약물이면 그대로 둬도 돼요",
              [i, o.idx, rest.at.idx].concat(G.dbl(a, ["O"]).map(function (t) { return t.idx; })));
          }
        }
        // 카바메이트
        var ns2 = G.sgl(a, ["N"]), os2 = G.sgl(a, ["O"]);
        if (ns2.length && os2.length) {
          add("metab", "mid", "카바메이트", "에스터보다는 안정하지만 가수분해와 카바모일화 반응이 일어날 수 있어요.", "아마이드·우레아로 교체 / 입체 장애 추가", [i, ns2[0].idx, os2[0].idx]);
        }
      }
      // 알데하이드
      if (a.el === "C" && !a.arom && a.h === 1 && G.dbl(a, ["O"]).length) {
        add("metab", "high", "알데하이드", "알데하이드 산화효소·ALDH로 빠르게 산화되고, 단백질 라이신과 결합해 독성을 낼 수 있어요.", "알코올·카복실산·옥심 등으로 바꾸거나 고리에 포함시키기", [i].concat(G.dbl(a, ["O"]).map(function (t) { return t.idx; })));
      }
      // 아닐린
      if (a.el === "N" && !a.arom && a.h >= 1 && G.arylNbr(a) && !G.amide(a)) {
        add("metab", "high", "아닐린 (방향족 1·2차 아민)",
          "N-수산화 → 나이트로소·퀴논이민 같은 반응성 대사체로 이어질 수 있어요. 특이체질 독성 경고로 널리 알려져 있어요.",
          "아마이드·설폰아마이드로 막기 / 오쏘 위치에 치환기를 넣어 입체 장애 주기 / 전자 끄는 기로 산화를 어렵게 하기",
          [i].concat(a.nbrs.map(function (n) { return n.atom; })));
      }
      // 나이트로 방향족
      if (a.el === "N" && a.charge === 1 && G.dbl(a, ["O"]).length && G.sgl(a, ["O"]).some(function (t) { return t.charge === -1; }) && G.arylNbr(a)) {
        add("metab", "high", "나이트로 방향족",
          "장내 세균과 나이트로환원효소가 환원시켜 나이트로소·하이드록실아민을 만들고, 유전독성 위험이 보고돼요.",
          "아마이드·나이트릴·설폰·CF₃ 등 전자 끄는 다른 기로 교체",
          [i].concat(a.nbrs.map(function (n) { return n.atom; })));
      }
      // 말단 알킨
      if (a.el === "C" && a.h === 1 && G.trp(a).length) {
        add("metab", "high", "말단 알킨", "CYP가 산화하면 케텐이 생겨 heme에 비가역적으로 결합해요(기전 기반 억제). 대사와 CYP 억제 모두 문제예요.", "내부 알킨으로 바꾸거나 나이트릴·사이클로프로필로 교체", [i].concat(G.trp(a).map(function (t) { return t.idx; })));
      }
      // 마이클 수용체
      if (a.el === "C" && !a.arom && G.dbl(a, ["C"]).length && a.nbrs.some(function (n) { var t = A[n.atom]; return n.bond.order === 1 && !n.bond.arom && (G.isCarbonyl(t) || (t.el === "C" && G.trp(t).length)); })) {
        add("metab", "mid", "마이클 수용체 (α,β-불포화 카보닐)",
          "글루타티온·단백질 시스테인과 공유 결합해요. 의도한 공유 결합 억제제가 아니라면 비특이 결합 위험이 있어요.",
          "포화시키기 / 반응성을 낮추는 치환기 넣기 / 공유 결합 억제제라면 선택성 확인 필요",
          [i].concat(G.dbl(a, ["C"]).map(function (t) { return t.idx; })));
      }
      // 에폭사이드
      if (a.el === "O" && a.nbrs.length === 2 && a.ringSizes.indexOf(3) >= 0) {
        add("metab", "high", "에폭사이드", "친핵체(DNA·단백질)와 바로 반응하는 알킬화제예요.", "다이올·옥세탄·사이클로프로판으로 교체", [i].concat(a.nbrs.map(function (n) { return n.atom; })));
      }
      // 알킬 할라이드
      if (["Cl", "Br", "I"].indexOf(a.el) >= 0 && a.nbrs.length === 1 && !A[a.nbrs[0].atom].arom) {
        add("metab", "mid", "지방족 할로젠 (알킬화 가능)", "SN2 알킬화와 글루타티온 포합이 일어날 수 있어요.", "F로 바꾸거나 방향족 고리에 붙이기", [i, a.nbrs[0].atom]);
      }
      // 하이드라진·하이드라존
      if (a.el === "N" && a.nbrs.some(function (n) { return A[n.atom].el === "N" && !n.bond.arom; }) && !a.arom) {
        var nn = a.nbrs.filter(function (n) { return A[n.atom].el === "N" && !n.bond.arom; })[0];
        if (nn && !G.amide(a)) add("metab", "mid", "하이드라진·하이드라존 (N–N)", "반응성 대사체와 간독성 경고가 알려진 구조예요.", "아마이드·옥사다이아졸·트라이아졸 등으로 교체", [i, nn.atom]);
      }
      // 사이클로프로필아민
      if (a.el === "N" && !a.arom && a.nbrs.some(function (n) { var t = A[n.atom]; return t.el === "C" && t.ringSizes.indexOf(3) >= 0; })) {
        var cp = a.nbrs.filter(function (n) { return A[n.atom].ringSizes.indexOf(3) >= 0; })[0];
        add("metab", "mid", "사이클로프로필아민", "한 전자 산화 후 고리가 열려 라디칼이 생기고, 기전 기반 억제로 이어질 수 있어요.", "사이클로부틸·이소프로필로 교체하거나 질소를 아마이드로 막기", [i, cp.atom]);
      }
      // 페놀
      if (a.el === "O" && a.h === 1 && G.arylNbr(a)) {
        add("metab", "low", "페놀 (2상 대사)", "글루쿠론산·황산 포합으로 빠르게 배설돼 노출량이 낮아질 수 있어요. 1상 대사는 아니지만 경구 생체이용률에 영향이 커요.",
          "메틸·전구약물로 가리기 / 오쏘 치환으로 포합 늦추기 / 등가체(아마이드, 설폰아마이드, 테트라졸) 검토", [i].concat(a.nbrs.map(function (n) { return n.atom; })));
      }
    });
    // 티오펜·퓨란·피롤
    ringsOfSize(m, 5).forEach(function (r) {
      if (!aromRing(m, r)) return;
      var het = r.filter(function (i) { return A[i].el !== "C"; });
      if (het.length !== 1) return;
      var e = A[het[0]].el;
      if (e === "S" || e === "O") {
        add("metab", "high", e === "S" ? "티오펜" : "퓨란",
          "CYP가 에폭사이드·S-옥사이드를 만들고, 이것이 단백질과 공유 결합해요. 간독성·기전 기반 억제 경고로 유명해요.",
          "피라졸·이속사졸·티아졸·페닐 등으로 교체 / 전자 끄는 기를 붙여 산화를 늦추기", r.slice());
      }
    });
    // para 비어 있는 단일 치환 페닐
    ringsOfSize(m, 6).forEach(function (r) {
      if (!aromRing(m, r) || r.some(function (i) { return A[i].el !== "C"; })) return;
      var sub = r.filter(function (i) { return A[i].nbrs.some(function (n) { return r.indexOf(n.atom) < 0; }); });
      if (sub.length !== 1) return;
      var start = r.indexOf(sub[0]), para = r[(start + 3) % 6];
      if (A[para].h >= 1) {
        add("metab", "mid", "단일 치환 페닐 (para 비어 있음)",
          "전자가 풍부한 벤젠 고리의 para 자리는 CYP 수산화가 가장 잘 일어나는 자리예요.",
          "para에 F·Cl·CF₃를 넣어 막기 / 피리딘 같은 전자 부족 고리로 교체",
          [para, sub[0]]);
      }
    });

    // ===== CYP =====
    // 메틸렌다이옥시페닐
    A.forEach(function (a, i) {
      if (a.el !== "C" || a.arom || a.h !== 2) return;
      var os = G.sgl(a, ["O"]); if (os.length !== 2) return;
      if (os.every(function (o) { return o.nbrs.some(function (n) { return A[n.atom].arom; }); })) {
        add("cyp", "high", "메틸렌다이옥시페닐 (MDP)",
          "CYP가 카빈 중간체를 만들어 heme 철과 단단히 결합해요. 시간 의존적 억제(TDI)의 교과서 사례예요.",
          "다이메톡시·에틸렌다이옥시·다이플루오로메틸렌다이옥시로 교체",
          [i].concat(os.map(function (o) { return o.idx; })));
      }
    });
    // 이미다졸·트리아졸 (heme 배위)
    ringsOfSize(m, 5).forEach(function (r) {
      if (!aromRing(m, r)) return;
      var ns = r.filter(function (i) { return A[i].el === "N"; });
      var free = ns.filter(function (i) { return A[i].h === 0 && A[i].nbrs.length === 2 && A[i].charge === 0; });
      if (ns.length >= 2 && free.length >= 1) {
        add("cyp", "high", ns.length === 2 ? "이미다졸·피라졸형 질소" : "트리아졸형 질소",
          "고리 질소의 비공유 전자쌍이 CYP heme 철과 직접 배위해요. 아졸 항진균제가 강한 CYP3A4 억제를 내는 이유예요.",
          "질소에 입체 장애를 주거나(2-치환) 피라졸·옥사졸·티아졸처럼 배위력이 약한 고리로 교체",
          r.slice());
      }
    });
    // 피리딘형 질소
    ringsOfSize(m, 6).forEach(function (r) {
      if (!aromRing(m, r)) return;
      r.forEach(function (i) {
        var a = A[i];
        if (a.el !== "N" || a.h !== 0 || a.charge !== 0 || a.nbrs.length !== 2) return;
        var ortho = a.nbrs.filter(function (n) { return r.indexOf(n.atom) >= 0; }).map(function (n) { return A[n.atom]; });
        var hindered = ortho.some(function (t) { return t.nbrs.some(function (n) { return r.indexOf(n.atom) < 0; }); });
        add("cyp", hindered ? "low" : "mid", "피리딘형 질소" + (hindered ? " (옆자리 치환 있음)" : ""),
          "방향족 질소도 heme 철과 배위할 수 있어요. 이미다졸보다는 약하지만 CYP 억제 경고로 봐요." + (hindered ? " 옆자리에 치환기가 있어 배위가 어려워 위험은 낮은 편이에요." : ""),
          hindered ? "지금처럼 옆자리 치환을 유지하세요." : "질소 옆(2 또는 6 위치)에 치환기를 넣어 접근을 막기 / 페닐·피리미딘으로 교체 검토",
          r.slice());
      });
    });
    // 친유성 기반
    if (p && p.clogp >= 3.5) {
      add("cyp", p.clogp >= 4.5 ? "mid" : "low", "높은 친유성 (cLogP " + p.clogp.toFixed(1) + ")",
        "cLogP가 높을수록 CYP(특히 3A4·2C9)의 넓은 소수성 포켓에 비특이적으로 잘 끼어들어요. 활성과 무관한 억제가 자주 나와요.",
        "극성기 추가 / 방향족 고리 수 줄이기 / Fsp3 높이기 (지금 " + (p.fsp3 != null ? p.fsp3.toFixed(2) : "–") + ")", []);
    }

    // ===== hERG =====
    var io = Chem.ionizable(m), strongBase = io.bases.filter(function (b) { return !b.quat; });
    var hergPts = [], hergAtoms = [];
    if (strongBase.length) {
      hergPts.push("염기성 질소 " + strongBase.length + "개 (" + strongBase.map(function (b) { return b.kind; }).join(", ") + ")");
      strongBase.forEach(function (b) { hergAtoms.push(b.atom); });
    }
    if (p && p.clogp > 3.7) hergPts.push("cLogP " + p.clogp.toFixed(1) + " (3.7 초과)");
    if (p && p.aromRings >= 2) hergPts.push("방향족 고리 " + p.aromRings + "개");
    if (p && p.rotb >= 6) hergPts.push("회전 결합 " + p.rotb + "개 (유연한 구조)");
    if (p && p.tpsa < 75) hergPts.push("TPSA " + p.tpsa.toFixed(0) + " (75 미만, 극성 낮음)");
    var hergN = hergPts.length;
    if (strongBase.length && hergN >= 3) {
      add("herg", hergN >= 4 ? "high" : "mid", "hERG 결합 요소 " + hergN + "가지",
        "hERG 통로는 염기성 질소(Tyr652·Phe656와 상호작용)와 친유성 방향족 고리를 좋아해요. 이 조합이 많을수록 위험이 커져요.\n해당: " + hergPts.join(" · "),
        "아민 염기성 낮추기(옆에 F 넣기, 고리 안에 산소 넣어 모폴린화) / cLogP 낮추기 / 카복실산·테트라졸 같은 산성기 추가 / 방향족 고리를 포화 고리로 바꾸기",
        hergAtoms);
    } else if (strongBase.length && hergN >= 2) {
      add("herg", "low", "hERG 결합 요소 " + hergN + "가지",
        "염기성 질소가 있지만 다른 위험 요소는 적은 편이에요.\n해당: " + hergPts.join(" · "), "지금 수준을 유지하면서 친유성이 올라가지 않게 주의하세요.", hergAtoms);
    } else if (!strongBase.length && p && p.clogp > 4.5 && p.aromRings >= 3) {
      add("herg", "low", "염기성 질소 없음 · 친유성 높음",
        "염기성 질소가 없어 전형적인 hERG 결합 양상은 아니지만, 친유성이 매우 높으면 중성 분자도 결합하는 사례가 있어요.\n해당: cLogP " + p.clogp.toFixed(1) + " · 방향족 고리 " + p.aromRings + "개", "친유성 낮추기", []);
    }

    // ===== DDI 요인 =====
    var cypHits = out.filter(function (x) { return x.cat === "cyp" && x.sev !== "low"; });
    var pgp = [];
    if (p) {
      if (p.mw > 400) pgp.push("MW " + p.mw.toFixed(0));
      if (p.no > 8) pgp.push("N+O " + p.no + "개");
      if (strongBase.length) pgp.push("염기성 질소");
      if (p.tpsa > 90) pgp.push("TPSA " + p.tpsa.toFixed(0));
    }
    if (cypHits.length) {
      add("ddi", cypHits.some(function (x) { return x.sev === "high"; }) ? "high" : "mid", "CYP 억제를 통한 상호작용 가능",
        "위의 CYP 경고(" + cypHits.map(function (x) { return x.name; }).join(", ") + ") 때문에 같은 효소로 대사되는 약과 함께 쓰면 그 약의 혈중 농도를 올릴 수 있어요.",
        "CYP 억제 시험(IC50, 시간 의존적 억제 포함)으로 확인 / 억제 경고 구조를 바꾸기", []);
    }
    if (pgp.length >= 3) {
      add("ddi", "low", "P-gp 기질 가능성",
        "P-gp 기질에 흔한 특징이 " + pgp.length + "가지 보여요(" + pgp.join(", ") + "). P-gp 억제제와 함께 쓰면 노출량이 달라지고, 뇌 투과도 영향을 받아요.",
        "Caco-2 양방향 투과도 시험으로 확인", []);
    }
    io.acids.forEach(function (ac) {
      add("ddi", "low", "산성기 (" + ac.kind + ")",
        "산성 약물은 혈장 알부민에 강하게 결합하고 OATP·OAT 수송체와 CYP2C9의 기질이 되는 경우가 많아, 와파린 같은 약과 상호작용 사례가 많아요.",
        "필요한 산성기라면 단백 결합률과 2C9 억제를 확인하세요.", [ac.atom]);
    });

    var order = { high: 0, mid: 1, low: 2 };
    out.sort(function (a, b) { return order[a.sev] - order[b.sev] || a.cat.localeCompare(b.cat); });
    return out;
  }

  // ---------- 방향족 결합을 번갈아 이중 결합으로 (그림용 케쿨레 배정) ----------
  function kekulize(m) {
    var A = m.atoms, dbl = {};
    var aromB = m.bonds.filter(function (b) { return b.arom; });
    if (!aromB.length) return dbl;
    var need = {}, adj = {};
    A.forEach(function (a, i) {
      if (!a.arom) return;
      // 이미 고리 밖 이중 결합이 있거나, 전자쌍을 내놓는 헤테로 원자는 이중 결합이 필요 없음
      if (a.nbrs.some(function (n) { return !n.bond.arom && n.bond.order === 2; })) return;
      if (a.el === "O" || a.el === "S" || a.el === "Se") return;
      if ((a.el === "N" || a.el === "P") && (a.h > 0 || a.nbrs.length === 3)) return;
      if (a.el === "C" && a.charge !== 0) return;
      if (a.el === "B") return;
      need[i] = true; adj[i] = [];
    });
    aromB.forEach(function (b) { if (need[b.a] && need[b.b]) { adj[b.a].push({ to: b.b, bond: b }); adj[b.b].push({ to: b.a, bond: b }); } });
    var ids = Object.keys(need).map(Number).sort(function (x, y) { return adj[x].length - adj[y].length; });
    var used = {};
    function solve(k) {
      while (k < ids.length && used[ids[k]]) k++;
      if (k >= ids.length) return true;
      var i = ids[k], opts = adj[i];
      for (var q = 0; q < opts.length; q++) {
        var o = opts[q]; if (used[o.to]) continue;
        used[i] = used[o.to] = 1; dbl[o.bond.idx] = 1;
        if (solve(k + 1)) return true;
        delete used[i]; delete used[o.to]; delete dbl[o.bond.idx];
      }
      return false;
    }
    if (!solve(0)) return null; // 완전한 배정 실패 → 점선으로 표시
    return dbl;
  }

  // ---------- 2D 구조 그림 ----------
  function depict(m, opt) {
    opt = opt || {};
    var hi = opt.highlight || {}, W = opt.width || 600, H = opt.height || 420;
    var A = m.atoms;
    if (!A.length) return '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '"></svg>';
    var has = A.some(function (a) { return (a.x || a.y) && isFinite(a.x); });
    if (!has) return null;
    var xs = A.map(function (a) { return a.x || 0; }), ys = A.map(function (a) { return -(a.y || 0); });
    var bl = [];
    m.bonds.forEach(function (b) { bl.push(Math.hypot(xs[b.a] - xs[b.b], ys[b.a] - ys[b.b])); });
    bl.sort(function (x, y) { return x - y; });
    var unit = bl.length ? bl[Math.floor(bl.length / 2)] || 1 : 1;
    var S = 34 / unit, pad = 26;
    var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs), minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
    var cw = (maxX - minX) * S, ch = (maxY - minY) * S;
    // 작은 분자는 키우고, 큰 분자는 줄여서 그림 칸을 알맞게 채우기
    var sc = Math.min((W - pad * 2) / Math.max(cw, 1), (H - pad * 2) / Math.max(ch, 1));
    sc = Math.max(0.25, Math.min(sc, 1.25));
    S *= sc; cw *= sc; ch *= sc;
    W = Math.min(W, Math.max(240, cw + pad * 2)); H = Math.min(H, Math.max(170, ch + pad * 2));
    var ox = (W - cw) / 2 - minX * S, oy = (H - ch) / 2 - minY * S;
    var P = A.map(function (a, i) { return { x: xs[i] * S + ox, y: ys[i] * S + oy }; });

    var COLOR = { N: "#1b53c4", O: "#c0392b", S: "#b8860b", F: "#17876e", Cl: "#17876e", Br: "#8c4a1f", I: "#7b3fa0", P: "#c97a09", B: "#8a6d3b", Si: "#777" };
    var ink = opt.ink || "#17222a";
    var kek = null; try { kek = kekulize(m); } catch (e) { kek = null; }
    function label(a) {
      if (a.el === "C" && a.charge === 0 && a.nbrs.length > 0) return null;
      var t = a.el;
      if (a.h === 1) t += "H"; else if (a.h > 1) t += "H" + a.h;
      if (a.charge) t += a.charge > 0 ? (a.charge > 1 ? "+" + a.charge : "+") : (a.charge < -1 ? a.charge : "−");
      return t;
    }
    var lab = A.map(label);
    var ringCenter = {};
    (m.rings || []).forEach(function (r) {
      var cx = 0, cy = 0; r.forEach(function (i) { cx += P[i].x; cy += P[i].y; }); cx /= r.length; cy /= r.length;
      for (var k = 0; k < r.length; k++) { var a = r[k], b = r[(k + 1) % r.length]; ringCenter[Math.min(a, b) + "-" + Math.max(a, b)] = { x: cx, y: cy }; }
    });
    function shrink(i, j) {
      var a = P[i], b = P[j], dx = b.x - a.x, dy = b.y - a.y, L = Math.hypot(dx, dy) || 1;
      var s1 = lab[i] ? 10 + (lab[i].length - 1) * 2.5 : 0, s2 = lab[j] ? 10 + (lab[j].length - 1) * 2.5 : 0;
      return { x1: a.x + dx / L * s1, y1: a.y + dy / L * s1, x2: b.x - dx / L * s2, y2: b.y - dy / L * s2, dx: dx / L, dy: dy / L };
    }
    return build();

    function build() {
      var g1 = [], g2 = [], g3 = [];
      Object.keys(hi).forEach(function (k) {
        var i = +k, cols = [].concat(hi[k]);
        cols.forEach(function (c, n) {
          g1.push('<circle cx="' + P[i].x.toFixed(1) + '" cy="' + P[i].y.toFixed(1) + '" r="' + (13 + n * 4) + '" fill="' + (n ? "none" : c) + '" fill-opacity="0.2" stroke="' + c + '" stroke-opacity="0.85" stroke-width="' + (n ? 2 : 1.6) + '"/>');
        });
      });
      m.bonds.forEach(function (b) {
        var s = shrink(b.a, b.b), px = -s.dy, py = s.dx, key = Math.min(b.a, b.b) + "-" + Math.max(b.a, b.b);
        var L = function (o, extra) { return '<line x1="' + (s.x1 + px * o).toFixed(1) + '" y1="' + (s.y1 + py * o).toFixed(1) + '" x2="' + (s.x2 + px * o).toFixed(1) + '" y2="' + (s.y2 + py * o).toFixed(1) + '"' + (extra || "") + "/>"; };
        var ord = b.order;
        if (b.arom) ord = kek ? (kek[b.idx] ? 2 : 1) : 1.5;
        if (ord === 2 || ord === 1.5) {
          var rc = ringCenter[key], side = 0;
          if (rc) { var mx = (s.x1 + s.x2) / 2, my = (s.y1 + s.y2) / 2; side = ((rc.x - mx) * px + (rc.y - my) * py) > 0 ? 1 : -1; }
          if (side) {
            g2.push(L(0));
            g2.push('<line x1="' + (s.x1 + px * side * 5 + s.dx * 6).toFixed(1) + '" y1="' + (s.y1 + py * side * 5 + s.dy * 6).toFixed(1) + '" x2="' + (s.x2 + px * side * 5 - s.dx * 6).toFixed(1) + '" y2="' + (s.y2 + py * side * 5 - s.dy * 6).toFixed(1) + '"' + (ord === 1.5 ? ' stroke-dasharray="5 3"' : "") + "/>");
          } else { g2.push(L(2.6)); g2.push(L(-2.6)); }
        } else if (ord === 3) { g2.push(L(0)); g2.push(L(3.4)); g2.push(L(-3.4)); }
        else g2.push(L(0));
      });
      A.forEach(function (a, i) {
        var t = lab[i]; if (!t) return;
        var col = COLOR[a.el] || ink, fs = 15;
        // H 위치: 이웃 반대쪽
        var dx = 0; a.nbrs.forEach(function (n) { dx += P[n.atom].x - P[i].x; });
        var txt = t;
        if (a.h > 0 && dx > 2) { var hpart = t.slice(a.el.length).replace(/[+−\d]+$/, ""), rest = t.slice(a.el.length + hpart.length); txt = hpart + a.el + rest; }
        g3.push('<text x="' + P[i].x.toFixed(1) + '" y="' + (P[i].y + fs * 0.35).toFixed(1) + '" font-size="' + fs + '" font-family="Arial, Helvetica, sans-serif" font-weight="600" text-anchor="middle" fill="' + col + '" stroke="#fff" stroke-width="3.5" paint-order="stroke" style="paint-order:stroke">' + txt + "</text>");
      });
      return '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + " " + H + '">' +
        "<g>" + g1.join("") + "</g>" +
        '<g stroke="' + ink + '" stroke-width="1.7" stroke-linecap="round">' + g2.join("") + "</g>" +
        g3.join("") + "</svg>";
    }
  }

  Chem.alerts = alerts;
  Chem.depict = depict;
  Chem.ALERT_CAT = CAT;
})(typeof window !== "undefined" ? window : globalThis);
