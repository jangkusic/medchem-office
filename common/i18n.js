/* 의약화학 오피스 — 한국어 / English
 *
 * 한국어 문구를 그대로 열쇠로 씁니다. T("분자량") 은 한국어일 때 "분자량" 을,
 * 영어일 때 사전의 "Molecular weight" 를 돌려줘요. 사전에 없으면 원문 그대로 둡니다.
 * 화면에 이미 그려진 글자는 apply() 가 돌면서 바꿉니다. 원문은 따로 보관해 두므로
 * 한국어로 되돌릴 때 손실이 없어요.
 */
(function (root) {
  "use strict";

  var EN = {
    // ---- 머리말·탭 ----
    "구조 그리기 · 물성": "Draw · Properties",
    "화합물 목록 · SAR": "Compounds · SAR",
    "위험 신호": "Alerts",
    "3D · 도킹": "3D · Docking",
    "준비 중": "coming soon",
    "도구": "Tools",
    "물성": "Properties",
    "도구 목록으로": "Back to tool list",
    "인터넷 없이 이 PC에서만 계산 · 저장": "Everything runs and stays on this computer",
    "도구상자 옮기기": "Move toolbars",
    "단축키 설정": "Shortcuts",

    // ---- 도구상자 옮기기 ----
    "도구상자를 끌어서 옮기세요": "Drag the toolbars where you want them",
    "점선 상자를 원하는 곳으로 끌면 돼요. 위치는 다음에 열 때도 그대로예요.":
      "Drag the dashed boxes anywhere. Positions are remembered next time.",
    "처음 위치로": "Reset positions",
    "완료": "Done",

    // ---- 왼쪽 빈 화면 ----
    "구조 편집기를 여는 중이에요. 처음 열 때는 10초쯤 걸려요.":
      "Opening the structure editor. The first load takes about 10 seconds.",
    "왼쪽에 구조를 그리면 물성이 바로 계산돼요.": "Draw a structure on the left and properties appear here.",
    "아래 칸에 SMILES나 MOL을 붙여넣어도 돼요.": "You can also paste SMILES or MOL below.",

    // ---- 조각(fragment) ----
    "캔버스에 있는 구조": "Structures on the canvas",
    "구조를 누르면 캔버스에서도 그 구조가 선택돼요.": "Click a structure to select it on the canvas too.",
    "개 모두 따로 저장": " — save each separately",

    // ---- 물성 ----
    "분자량 (MW)": "Molecular weight (MW)",
    "정확 질량 (M)": "Exact mass (M)",
    "약물 유사성": "Drug-likeness",
    "통과": "pass",
    "미달": "fail",
    "위반 ": "violations: ",
    "개": "",
    "회전 결합": "Rotatable bonds",
    "회전 결합 ": "Rotatable bonds ",
    "분자량": "Molecular weight",
    "Lipinski 수소결합 주개": "Lipinski H-bond donors",
    "Lipinski 수소결합 받개": "Lipinski H-bond acceptors",
    "몰 굴절률 (MR)": "Molar refractivity (MR)",
    "무거운 원자 수": "Heavy atoms",
    "고리 / 방향족 고리": "Rings / aromatic rings",
    "헤테로 원자": "Heteroatoms",
    "전하": "Formal charge",
    "참고": "info",
    "계산할 수 없어요": "could not be calculated",

    // ---- CNS MPO ----
    "실험값이 있으면 직접 넣기": "Enter measured values instead",
    "pKa (가장 센 염기)": "pKa (most basic)",
    "cLogD (pH 7.4)": "cLogD (pH 7.4)",
    "pKa (염기)": "pKa (base)",
    "예상값 사용": "using estimate",
    "예상": "estimated",
    "직접 입력": "entered",
    "넣어주신 값을 씁니다": "using the value you entered",
    "이온화되는 염기가 없어요": "no ionizable base",
    " · 염기 pKa ": " · base pKa ",
    " · 산 pKa ": " · acid pKa ",
    " 로 계산": "",
    "찾은 이온화 작용기: ": "Ionizable groups found: ",
    "≥ 4 좋음": "≥ 4 good",
    "보통": "moderate",
    "낮음": "low",
    "≤ 3 만점, ≥ 5 0점": "1.0 at ≤ 3, 0 at ≥ 5",
    "≤ 2 만점, ≥ 4 0점": "1.0 at ≤ 2, 0 at ≥ 4",
    "≤ 360 만점, ≥ 500 0점": "1.0 at ≤ 360, 0 at ≥ 500",
    "40–90 만점, < 20 또는 > 120 0점": "1.0 at 40–90, 0 below 20 or above 120",
    "0 만점, ≥ 3.5 0점": "1.0 at 0, 0 at ≥ 3.5",
    "≤ 8 만점, ≥ 10 0점": "1.0 at ≤ 8, 0 at ≥ 10",

    // ---- 입체화학 ----
    "입체화학": "Stereochemistry",
    "입체 요소가 없어요": "No stereocentres",
    "비칭 탄소": "stereocentre",
    "이중결합": "double bond",
    "지정됨": "defined",
    "지정 안 됨": "undefined",
    "가능한 입체이성질체": "possible stereoisomers",
    "모든 이성질체 보기": "Show all isomers", "그리던 구조 선택": "Select the structure being drawn", "지금 구조": "current",
    "이 구조의 입체이성질체": "Stereoisomers of this structure",
    "편집기로 보내기": "Send to editor",
    "입체 지정이 없어 이성질체를 셀 수 없어요": "No stereo information to enumerate",
    "개 중 ": " of ",
    "개만 보여요": " shown",

    // ---- 3D ----
    "3D 보기": "View in 3D",
    "3D 좌표를 만드는 중…": "Building 3D coordinates…",
    "무거운 원자 {n}개 · 평면에서 {f}Å 벗어남 · {ms}ms":
      "{n} heavy atoms · {f} Å out of plane · {ms} ms",
    "수소 표시": "Show hydrogens",
    "공-막대": "Ball & stick",
    "다시 계산": "Recalculate",
    "구조 {n}": "Structure {n}",
    "먼저 구조를 그려주세요": "Draw a structure first",

    // ---- 복사·저장 ----
    "복사": "Copy",
    "HRMS 문구": "HRMS text",
    "계산 중…": "calculating…",
    "복사했어요": "Copied",
    "Ctrl+C로 복사해주세요": "Press Ctrl+C to copy",
    "MOL 저장": "Save MOL",
    "ChemDraw(.cdxml) 저장": "Save ChemDraw (.cdxml)",
    "화합물 목록에 저장": "Save to compound list",
    "화합물 코드": "Compound code",
    "타깃": "Target",
    "활성값": "Activity",
    "메모": "Note",
    "목록에 저장": "Save to list",
    "모두 따로 저장": "Save each separately",
    "예: NVR-001": "e.g. NVR-001",
    "예: EGFR": "e.g. EGFR",
    "예: 35": "e.g. 35",
    "예: IC50 (EGFR), 3단계 수율 42%": "e.g. IC50 (EGFR), 42% over 3 steps",
    "이름": "Name",
    "취소": "Cancel",
    "닫기": "Close",
    "저장했어요": " saved",
    "지웠어요": "Deleted",
    "모두 지웠어요": "All deleted",
    "정말 삭제": "Really delete",
    "삭제": "Delete",
    "편집기로": "Open",
    "한 번 더 누르면 모두 삭제": "Press again to delete everything",

    // ---- 불러오기 ----
    "SMILES · MOL · CDXML 불러오기": "Load SMILES · MOL · CDXML",
    "편집기에 넣기": "Load into editor",
    "파일 열기": "Open file",
    "예: CC(=O)Oc1ccccc1C(=O)O": "e.g. CC(=O)Oc1ccccc1C(=O)O",

    // ---- 목록 탭 ----
    "유사도 기준": "Similarity reference",
    "지금 그린 구조로": "Use current structure",
    "기준 지우기": "Clear",
    "부분 구조 검색 (SMILES, 비우면 전체)": "Substructure search (SMILES, blank = all)",
    "그린 구조": "Current structure",
    "지금 그린 구조를 검색어로": "Search with the structure you drew",
    "예: c1ccc2[nH]ccc2c1 (인돌)": "e.g. c1ccc2[nH]ccc2c1 (indole)",
    "파일": "File",
    "CSV 내보내기": "Export CSV",
    "CSV·SMILES 불러오기": "Import CSV / SMILES",
    "모두 지우기": "Clear all",
    "구조": "Structure",
    "코드": "Code",
    "활성": "Activity",
    "유사도": "Similarity",
    "경고": "Alerts",
    "열 이름 바꾸기": "Rename column",
    "새 이름": "New name",
    "바꾸기": "Rename",
    "계산할 수 없는 구조예요": "This structure could not be calculated",
    "열 제목을 누르면 정렬돼요": "click a column header to sort",
    "아직 저장한 화합물이 없어요. 구조 탭에서 구조를 그리고 ":
      "No compounds yet. Draw a structure on the Draw tab and press ",
    "을 누르거나, 위의 ": ", or use ",
    "로 한 번에 넣을 수 있어요.": " above to load several at once.",
    "전체 ": "",
    "개 찾음 (전체 ": " found (of ",
    "개)": ")",

    // ---- 위험 신호 ----
    "구조 그리기 탭에서 구조를 그리면 여기에 위험 신호가 표시돼요.":
      "Draw a structure on the Draw tab and alerts will appear here.",
    "주의": "high",
    "개 · 보통 ": " · moderate ",
    "개 (누르면 위험 신호 탭에서 자세히)": " (click for details on the Alerts tab)",

    // ---- 준비 중 탭 ----
    "NMR 처리 (다음 단계)": "NMR processing (next)",
    "Mnova에서 매일 쓰는 1D ¹H NMR 기능을 붙일 예정이에요.":
      "The everyday 1D ¹H NMR workflow from Mnova.",
    "Bruker fid, JCAMP-DX 불러오기": "Read Bruker fid and JCAMP-DX",
    "푸리에 변환, 위상·기준선 보정, 용매 피크로 ppm 맞추기":
      "Fourier transform, phase and baseline correction, referencing on the solvent peak",
    "적분, 다중선 분석(J 값), 보고서 문구 자동 생성":
      "Integration, multiplet analysis (J values), report text",
    "구조 탭에서 그린 구조와 함께 보기": "Side by side with the structure you drew",
    "3D 구조 · 도킹 (다음 단계)": "3D structure · docking (next)",
    "그린 구조의 3D 모양 보기, 단백질(PDB) 불러오기": "3D view of your structure, load a protein (PDB)",
    "AutoDock Vina(브라우저용)로 소수 화합물 도킹, 결합 자세 보기":
      "Dock a few compounds with AutoDock Vina, inspect binding poses",
    "MD는 입력 파일 생성과 결과 궤적 보기": "For MD: build input files and view the trajectory",

    // ---- 단축키 ----
    "내 단축키 모두 지우기": "Clear all my shortcuts",
    "기능 찾기 (예: 이중, 벤젠, 지우개)": "Find a function (e.g. double, benzene, eraser)",
    "키를 누르세요… (Esc 취소)": "Press a key… (Esc to cancel)",
    "단축키 정하기": "Set shortcut",
    "지우기": "Clear",
    "내 단축키를 모두 지웠어요": "All your shortcuts were cleared",
    "한 번 더 누르면 모두 지워져요": "Press again to clear them all",
    "기능 옆의 칸을 누르고 원하는 키를 누르세요. Ctrl, Alt, Shift를 함께 눌러도 돼요. 한글 입력 상태에서도 동작해요. 내가 정한 단축키가 편집기 기본 단축키보다 우선해요.":
      "Click the box next to a function and press the key you want. Ctrl, Alt and Shift work too, and it works with a Korean keyboard layout. Your shortcut overrides the editor's default.",
    "ChemDraw 기본값으로": "Use ChemDraw defaults",
    "ChemDraw 단축키로 맞췄어요": "Switched to ChemDraw shortcuts",
    "도구는 지금 화면에 없어요": " is not on screen right now",

    // ---- 염 ----
    "물": "water", "메실산": "mesylate", "토실산": "tosylate", "푸마르산": "fumaric acid",
    "숙신산": "succinic acid", "타르타르산": "tartaric acid", "시트르산": "citric acid",
    "옥살산": "oxalic acid", "아세트산": "acetic acid", "메탄올": "methanol", "에탄올": "ethanol",


    // ---- pKa 근거 문구 (pka.js) ----
    "1차 아민": "primary amine", "2차 아민": "secondary amine", "3차 아민": "tertiary amine",
    "아닐린": "aniline", "아미딘": "amidine", "구아니딘": "guanidine",
    "피리딘": "pyridine", "다이아진": "diazine", "이미다졸": "imidazole", "피라졸": "pyrazole",
    "티아졸": "thiazole", "옥사졸": "oxazole", "트라이아졸": "triazole",
    "카복실산": "carboxylic acid", "페놀": "phenol", "설폰아마이드": "sulfonamide",
    "아실설폰아마이드": "acylsulfonamide", "테트라졸": "tetrazole", "싸이올": "thiol",
    "산소": "oxygen", "질소": "nitrogen", "나이트릴": "nitrile",
    "카보닐·설포닐": "carbonyl/sulfonyl", "방향족 고리": "aromatic ring",
    "전자 끄는 기 ": "electron-withdrawing groups: ",
    // ---- 남은 조각 ----
    "전체 ": "Total ",
    " · 열 제목을 누르면 정렬돼요": " · click a column header to sort",
    "유사도 기준: ": "similarity vs ",
    "찾는 기능이 없어요.": "No matching function.",
    "기본 ": "default ",
    "이 설정이 우선해요.": "Your setting takes priority.",
    "편집기 기본 단축키 ": "The editor's default shortcut ",
    "와 겹쳐요. ": " conflicts. ",

    // ---- 단축키 창의 기능 이름 ----
    "파일·편집": "File · Edit",
    "새로 그리기 (모두 지우기)": "New (clear all)", "다른 형식으로 저장": "Save as…",
    "붙여넣기": "Paste", "잘라내기": "Cut", "실행 취소": "Undo", "다시 실행": "Redo",
    "방향족으로 표시": "Aromatize", "케쿨레로 표시": "Dearomatize (Kekulé)",
    "자동 배치 (Layout)": "Layout", "구조 다듬기 (Clean Up)": "Clean up",
    "R/S·E/Z 표시": "Calculate CIP (R/S, E/Z)", "구조 검사": "Check structure",
    "계산값 창": "Calculated values", "수소 보이기/숨기기": "Show/hide hydrogens",
    "편집기 설정": "Editor settings",
    "선택·지우기": "Select · Erase",
    "손 (화면 이동)": "Hand (pan)", "사각형 선택": "Rectangle select",
    "올가미 선택": "Lasso select", "조각 선택": "Fragment select", "지우개": "Eraser",
    "결합": "Bonds", "단일 결합": "Single bond", "이중 결합": "Double bond",
    "삼중 결합": "Triple bond", "방향족 결합": "Aromatic bond",
    "쐐기 결합 (위)": "Wedge (up)", "점선 쐐기 (아래)": "Hash (down)",
    "물결 결합 (위/아래)": "Wavy (up or down)",
    "이중 결합 (cis/trans 모름)": "Double bond (cis/trans unknown)",
    "아무 결합": "Any bond", "배위 결합": "Dative bond", "수소 결합": "Hydrogen bond", "사슬": "Chain",
    "원자": "Atoms", "탄소 C": "Carbon C", "질소 N": "Nitrogen N", "산소 O": "Oxygen O",
    "황 S": "Sulfur S", "플루오린 F": "Fluorine F", "염소 Cl": "Chlorine Cl",
    "브로민 Br": "Bromine Br", "아이오딘 I": "Iodine I", "인 P": "Phosphorus P", "수소 H": "Hydrogen H",
    "주기율표": "Periodic table", "아무 원자": "Any atom",
    "고리": "Rings", "벤젠": "Benzene", "사이클로헥산": "Cyclohexane", "사이클로펜탄": "Cyclopentane",
    "사이클로펜타다이엔": "Cyclopentadiene", "사이클로프로판": "Cyclopropane", "사이클로부탄": "Cyclobutane",
    "사이클로헵탄": "Cycloheptane", "사이클로옥탄": "Cyclooctane", "구조 라이브러리": "Template library",
    "기타 도구": "Other tools", "전하 +": "Charge +", "전하 −": "Charge −",
    "입체 화학": "Stereochemistry", "S-그룹": "S-group", "R-그룹 표시": "R-group label",
    "R-그룹 조각": "R-group fragment", "연결점": "Attachment point",
    "반응 +": "Reaction plus", "반응 화살표": "Reaction arrow",
    "역합성 화살표": "Retrosynthetic arrow", "평형 화살표": "Equilibrium arrow",
    "자동 원자 매핑": "Automap", "글상자": "Text", "사각형": "Rectangle", "타원": "Ellipse",
    "의약화학 오피스": "MedChem Office",
    "SMILES 복사": "Copy SMILES", "HRMS 문구 복사": "Copy HRMS text",
    "pKa 입력칸으로 이동": "Jump to pKa field", "화합물 목록 탭 열기": "Open Compounds tab",
    "위험 신호 탭 열기": "Open Alerts tab", "구조 그리기 탭 열기": "Open Draw tab",
    "저장된 화합물 ": "Saved compounds: ", "도구 1": "Tool 1", "개 · ": " · ", "편집기": "editor", "기본값으로": "defaults",
    // ---- 메시지 ----
    "먼저 구조를 그려주세요": "Draw a structure first",
    "이 구조를 읽지 못했어요. 다시 그려주세요.": "Could not read this structure. Please redraw it.",
    "구조 탭에서 먼저 구조를 그려주세요": "Draw a structure on the Draw tab first",
    "저장 공간이 부족해요. CSV로 내보내 두세요.": "Storage is full. Export to CSV to keep your data.",
    "내보낼 화합물이 없어요": "No compounds to export"
  };

  var LANGS = { ko: "한국어", en: "English" };
  var KEY = "medchem-office-lang";
  var lang = "ko";
  try { var v = localStorage.getItem(KEY); if (v && LANGS[v]) lang = v; } catch (e) {}

  function T(s) {
    if (lang === "ko" || s == null) return s;
    var k = String(s);
    if (EN[k] != null) return EN[k];
    var t = k.trim();
    if (t !== k && EN[t] != null) return k.replace(t, EN[t]);   // 앞뒤 공백만 다른 경우
    return s;
  }

  // 원문 보관 — 한국어로 되돌릴 때 씁니다
  var orig = new WeakMap();
  function textOf(node) { if (!orig.has(node)) orig.set(node, node.nodeValue); return orig.get(node); }
  var ATTRS = ["placeholder", "title", "aria-label", "value"];

  function apply(where) {
    where = where || document.body;
    if (!where) return;
    var w = document.createTreeWalker(where, NodeFilter.SHOW_TEXT, null);
    var n, jobs = [];
    while ((n = w.nextNode())) jobs.push(n);
    jobs.forEach(function (node) {
      var p = node.parentNode;
      if (p && (p.tagName === "SCRIPT" || p.tagName === "STYLE")) return;
      var src = textOf(node);
      if (!src || !src.trim()) return;
      var out = lang === "ko" ? src : translateRun(src);
      if (node.nodeValue !== out) node.nodeValue = out;
    });
    var els = where.querySelectorAll("[placeholder],[title],[aria-label]");
    Array.prototype.forEach.call(els, function (el) {
      ATTRS.forEach(function (a) {
        if (!el.hasAttribute(a)) return;
        var key = "__i18n_" + a;
        if (el[key] == null) el[key] = el.getAttribute(a);
        var src = el[key];
        var out = lang === "ko" ? src : translateRun(src);
        if (el.getAttribute(a) !== out) el.setAttribute(a, out);
      });
    });
  }
  // 문구 전체가 사전에 있으면 그걸로 바꾸고, 없으면 "Lipinski 통과" 처럼
  // 조립된 글 안에서 사전 항목을 찾아 바꿔요. 다만 "물" 이 "물성" 안에서 바뀌면 안 되니
  // 앞뒤가 한글이면 건너뜁니다.
  var SUBS = null;
  function subsList() {
    if (SUBS) return SUBS;
    SUBS = Object.keys(EN).filter(function (k) { return k.trim().length >= 2 && /[가-힣]/.test(k); })
      .sort(function (a, b) { return b.length - a.length; });
    return SUBS;
  }
  function isHangul(ch) { return !!ch && /[가-힣]/.test(ch); }
  function translateRun(s) {
    var t = s.trim();
    if (!t) return s;
    if (EN[t] != null) return s.replace(t, EN[t]);
    if (!/[가-힣]/.test(s)) return s;
    var out = s;
    subsList().forEach(function (k) {
      var from = 0, i;
      while ((i = out.indexOf(k, from)) >= 0) {
        var before = out.charAt(i - 1), after = out.charAt(i + k.length);
        if (isHangul(before) || isHangul(after)) { from = i + 1; continue; }
        out = out.slice(0, i) + EN[k] + out.slice(i + k.length);
        from = i + EN[k].length;
      }
    });
    return out;
  }

  function set(l) {
    if (!LANGS[l]) return;
    lang = l;
    try { localStorage.setItem(KEY, l); } catch (e) {}
    document.documentElement.setAttribute("lang", l === "en" ? "en" : "ko");
    if (root.onLangChange) root.onLangChange(l);
    apply();
  }

  root.I18N = { get lang() { return lang; }, set: set, T: T, apply: apply, langs: LANGS, dict: EN };
  root.T = T;
})(window);
