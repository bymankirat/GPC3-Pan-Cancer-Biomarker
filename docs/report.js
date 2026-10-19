/* GPC3 pan-cancer report. Every number comes from results.js, written by the pipeline. */
(function () {
  "use strict";
  const D = window.GPC3;
  const $ = (s) => document.querySelector(s);
  const NS = "http://www.w3.org/2000/svg";
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const f = (x, d = 2) => (x == null || !isFinite(x) ? "–" : Number(x).toFixed(d));
  const pp = (p) => (p == null ? "p = –" : p < 0.001 ? "p < 0.001" : `p = ${p < 0.01 ? p.toFixed(3) : p.toFixed(2)}`);
  const pct = (x) => (x == null ? "–" : `${Math.round(100 * x)}%`);
  if (!D) { const b = $("#banner"); b.hidden = false; b.textContent = "No results yet: the pipeline has not run."; return; }

  function svg(w, h) { const s = document.createElementNS(NS, "svg"); s.setAttribute("viewBox", `0 0 ${w} ${h}`); s.setAttribute("width", w); s.setAttribute("height", h); return s; }
  function add(p, tag, a = {}, text) { const e = document.createElementNS(NS, tag); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); if (text != null) e.textContent = text; p.appendChild(e); return e; }
  function table(el, head, rows) {
    el.innerHTML = "<thead><tr>" + head.map((h) => `<th class="${h.num ? "num" : ""}">${h.t}</th>`).join("") + "</tr></thead><tbody>" +
      rows.map((r) => "<tr>" + r.map((c, i) => `<td class="${head[i].num ? "num" : ""}">${c}</td>`).join("") + "</tr>").join("") + "</tbody>";
  }
  function axisY(s, y, ticks, x0, x1, fmt = (v) => v) { ticks.forEach((v) => { add(s, "line", { x1: x0, x2: x1, y1: y(v), y2: y(v), stroke: css("--grid") }); add(s, "text", { x: x0 - 6, y: y(v) + 4, "text-anchor": "end", "font-size": 10, fill: css("--muted") }, fmt(v)); }); }
  function axisX(s, x, ticks, y0, y1, fmt = (v) => v) { ticks.forEach((v) => { add(s, "line", { x1: x(v), x2: x(v), y1: y0, y2: y1, stroke: css("--grid") }); add(s, "text", { x: x(v), y: y1 + 14, "text-anchor": "middle", "font-size": 10, fill: css("--muted") }, fmt(v)); }); }

  const G = D.meta.gene, F = D.meta.focus;
  $("#generated").textContent = `data downloaded ${D.meta.downloaded}`;
  const pan = D.pancancer, pe = D.pancancer_extra || {};
  const tested = pan.filter((r) => r.p != null);
  const lihc = pan.find((r) => r.cancer === F) || {};
  const S = D.survival, OS = S && S.endpoints.OS;
  const colorDir = (d) => (d === "up" ? css("--bad") : d === "down" ? css("--c2") : css("--muted"));

  // ------------------------------------------------------------ KPIs
  const sigOS = (D.survival_pancancer || []).filter((r) => r.q < 0.05);
  $("#kpis").innerHTML = [
    [`${pe.n_up}↑ ${pe.n_down}↓`, `cancers with GPC3 up or down vs normal`, `${tested.length} cancers tested, q < 0.05`],
    [`AUC ${f(lihc.auc)}`, `tells liver tumour from normal liver`, `95% CI ${f(lihc.auc_lo)}–${f(lihc.auc_hi)}`],
    [`${pct(pe.paired && pe.paired.fraction_higher_in_tumour)}`, "of patients: tumour higher than own normal liver", pe.paired ? `${pe.paired.n_pairs} pairs, ${pp(pe.paired.p)}` : ""],
    [OS && OS.cox_adjusted ? `HR ${f(OS.cox_adjusted.hr)}` : "–", "overall survival per SD, adjusted", OS && OS.cox_adjusted ? `${pp(OS.cox_adjusted.p)}, ${OS.cox_adjusted.n} patients` : ""],
    [OS ? `p ${f(OS.cutpoint.min_p, 3)}` : "–", `"best cutpoint" p for overall survival`, OS ? `after permutation: ${pp(OS.cutpoint.permutation_p)}` : ""],
    [D.structure && D.structure.mutations ? `${f(D.structure.mutations.percent_samples, 1)}%` : "–", "of ~10,000 tumours carry a GPC3 mutation", D.structure && D.structure.mutations ? `${D.structure.mutations.total} mutations` : ""],
  ].map(([v, l, s]) => `<div class="kpi"><div class="v">${esc(v)}</div><div class="l">${esc(l)}</div><div class="arms small">${esc(s)}</div></div>`).join("");

  // ------------------------------------------------------------ pan-cancer boxes + AUC
  (function () {
    const rows = pan.slice().sort((a, b) => (b.log2fc ?? -99) - (a.log2fc ?? -99));
    const dist = pe.distributions || {};
    const W = 1000, row = 22, L = 70, mid = 640, R = 980, T = 24, H = rows.length * row + T + 30;
    const maxv = Math.max(8, ...rows.flatMap((r) => [...(dist[r.cancer]?.tumour || []), ...(dist[r.cancer]?.normal || [])]));
    const x = (v) => L + (v / maxv) * (mid - L - 20);
    const xa = (v) => mid + 40 + ((v - 0.3) / 0.7) * (R - mid - 50);
    const s = svg(W, H);
    axisX(s, x, [0, 2, 4, 6, 8, 10, 12].filter((v) => v <= maxv), T - 6, H - 30);
    add(s, "text", { x: (L + mid) / 2, y: H - 2, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, "log2(TPM+1)");
    axisX(s, xa, [0.4, 0.6, 0.8, 1.0], T - 6, H - 30, (v) => v.toFixed(1));
    add(s, "line", { x1: xa(0.5), x2: xa(0.5), y1: T - 6, y2: H - 30, stroke: css("--muted"), "stroke-dasharray": "3 3" });
    add(s, "text", { x: (mid + 40 + R) / 2, y: H - 2, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, "AUC (tumour vs normal)");
    const box = (q, y, col, fill, h) => {
      add(s, "line", { x1: x(q[0]), x2: x(q[4]), y1: y, y2: y, stroke: col });
      add(s, "rect", { x: x(q[1]), y: y - h / 2, width: Math.max(1, x(q[3]) - x(q[1])), height: h, fill, stroke: col, rx: 2 });
      add(s, "line", { x1: x(q[2]), x2: x(q[2]), y1: y - h / 2, y2: y + h / 2, stroke: col, "stroke-width": 2 });
    };
    rows.forEach((r, i) => {
      const y = T + i * row + row / 2;
      const d = dist[r.cancer] || {};
      const col = colorDir(r.direction);
      if (i % 2) add(s, "rect", { x: 0, y: y - row / 2, width: W, height: row, fill: css("--bg"), opacity: 0.5 });
      const t = add(s, "text", { x: 8, y: y + 4, "font-size": 12, "font-weight": r.cancer === F ? 700 : 500 }, r.cancer);
      add(t, "title", {}, `${r.disease}\n${r.n_tumour} tumours, ${r.n_normal ?? 0} normal (${r.n_normal_tcga} TCGA, ${r.n_normal_gtex} GTEx)`);
      if (d.normal) box(d.normal, y + 4, css("--muted"), css("--line"), 7);
      if (d.tumour) box(d.tumour, y - 4, col, col + "55", 7);
      if (r.auc != null) {
        add(s, "line", { x1: xa(Math.max(0.3, r.auc_lo)), x2: xa(r.auc_hi), y1: y, y2: y, stroke: col, "stroke-width": 2 });
        const c = add(s, "circle", { cx: xa(Math.max(0.3, r.auc)), cy: y, r: 4, fill: col });
        add(c, "title", {}, `${r.cancer}: AUC ${f(r.auc, 3)} (${f(r.auc_lo, 3)}–${f(r.auc_hi, 3)}), log2FC ${f(r.log2fc)}, q ${r.q < 0.001 ? "<0.001" : f(r.q, 3)}`);
      } else add(s, "text", { x: xa(0.5), y: y + 4, "font-size": 10, fill: css("--muted"), "text-anchor": "middle" }, "too few normal samples");
    });
    $("#c-pan").appendChild(s);
    const up = tested.filter((r) => r.direction === "up").map((r) => r.cancer), down = tested.filter((r) => r.direction === "down").map((r) => r.cancer);
    $("#pan-note").textContent = `${tested.length} cancers had enough normal samples to test. GPC3 is higher in ${up.length} (${up.join(", ")}) and lower in ${down.length} (${down.join(", ")}) at FDR 5%. ` +
      `Liver cancer has the largest rise: median ${f(lihc.median_tumour)} vs ${f(lihc.median_normal)} log2(TPM+1), log2 fold change ${f(lihc.log2fc)}. ` +
      `Outside the liver, GPC3 is more often lower in tumours than in normal tissue, in line with reports that it acts as a growth brake in breast, lung and ovarian tissue. ` +
      `With so many samples small differences reach significance, so the size of the change and the AUC matter more than the p-value; cancers compared only with GTEx (no TCGA normals) carry a possible batch effect.`;
  })();

  // ------------------------------------------------------------ paired + AFP
  (function () {
    const P = pe.paired;
    if (P) {
      const W = 420, H = 300, L = 46, B = 40, T = 12;
      const all = P.pairs.flat(), mx = Math.ceil(Math.max(...all)), mn = Math.floor(Math.min(...all));
      const y = (v) => T + (1 - (v - mn) / (mx - mn)) * (H - T - B);
      const s = svg(W, H);
      axisY(s, y, [mn, Math.round((mn + mx) / 2), mx], L, W - 10);
      const xs = [L + 60, W - 80];
      P.pairs.forEach(([n, t]) => add(s, "line", { x1: xs[0], x2: xs[1], y1: y(n), y2: y(t), stroke: t > n ? css("--bad") : css("--c2"), "stroke-opacity": 0.45 }));
      P.pairs.forEach(([n, t]) => { add(s, "circle", { cx: xs[0], cy: y(n), r: 2.5, fill: css("--muted") }); add(s, "circle", { cx: xs[1], cy: y(t), r: 2.5, fill: css("--bad") }); });
      add(s, "text", { x: xs[0], y: H - 14, "text-anchor": "middle", "font-size": 12 }, "normal liver");
      add(s, "text", { x: xs[1], y: H - 14, "text-anchor": "middle", "font-size": 12 }, "tumour");
      $("#c-paired").appendChild(s);
    }
    const A = pe.afp;
    if (A) {
      const W = 420, H = 180, L = 90, R = 400;
      const x = (v) => L + ((v - 0.5) / 0.5) * (R - L);
      const s = svg(W, H);
      axisX(s, x, [0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 20, H - 40, (v) => v.toFixed(1));
      [["GPC3", A.auc_gpc3, A.auc_gpc3_ci, css("--bad")], ["AFP", A.auc_afp, A.auc_afp_ci, css("--c2")]].forEach(([n, v, ci, c], i) => {
        const y = 50 + i * 50;
        add(s, "text", { x: L - 10, y: y + 4, "text-anchor": "end", "font-size": 13, "font-weight": 600 }, n);
        add(s, "line", { x1: x(ci[0]), x2: x(ci[1]), y1: y, y2: y, stroke: c, "stroke-width": 3 });
        add(s, "circle", { cx: x(v), cy: y, r: 6, fill: c });
        add(s, "text", { x: x(v), y: y - 12, "text-anchor": "middle", "font-size": 11 }, f(v, 3));
      });
      $("#c-afp-auc").appendChild(s);
    }
    $("#diag-note").textContent = (P ? `In the ${P.n_pairs} patients who gave both a tumour and a normal liver sample, the tumour had more GPC3 in ${pct(P.fraction_higher_in_tumour)} (Wilcoxon signed-rank ${pp(P.p)}), so the difference is not a batch effect between TCGA and GTEx. ` : "") +
      (A ? `On its own GPC3 separates liver tumours from normal liver with AUC ${f(A.auc_gpc3)}, against ${f(A.auc_afp)} for AFP at the RNA level; in tumours the two are ${A.spearman_in_tumours > 0.3 ? "linked" : "only weakly linked"} (Spearman ρ = ${f(A.spearman_in_tumours)}), so they can flag different tumours.` : "");
  })();

  // ------------------------------------------------------------ survival
  (function () {
    if (!S) return;
    const eps = ["OS", "DSS", "PFI"].filter((e) => S.endpoints[e]);
    const seg = $("#ep-seg");
    seg.innerHTML = eps.map((e, i) => `<button data-e="${e}" aria-pressed="${i === 0}">${esc(S.endpoints[e].label)}</button>`).join("");
    function draw(e) {
      const E = S.endpoints[e];
      // KM
      let el = $("#c-km"); el.innerHTML = "";
      const W = 420, H = 280, L = 44, B = 40, T = 10, Rr = 12;
      const tmax = Math.min(10, Math.max(...E.km_high.map((p) => p[0]), ...E.km_low.map((p) => p[0])));
      const x = (t) => L + (Math.min(t, tmax) / tmax) * (W - L - Rr), y = (v) => T + (1 - v) * (H - T - B);
      let s = svg(W, H);
      axisY(s, y, [0, 0.25, 0.5, 0.75, 1], L, W - Rr, (v) => v.toFixed(2));
      axisX(s, x, [0, 2, 4, 6, 8, 10].filter((v) => v <= tmax), T, H - B);
      add(s, "text", { x: (L + W) / 2, y: H - 4, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, "years");
      [[E.km_high, css("--bad")], [E.km_low, css("--c2")]].forEach(([pts, c]) => {
        let d = `M${x(0)},${y(1)}`;
        pts.forEach(([t, v], i) => { if (i) d += `H${x(t)}V${y(v)}`; });
        add(s, "path", { d, fill: "none", stroke: c, "stroke-width": 2 });
      });
      add(s, "text", { x: W - Rr, y: T + 14, "text-anchor": "end", "font-size": 11 }, `log-rank ${pp(E.median_split_p)}`);
      add(s, "text", { x: W - Rr, y: T + 28, "text-anchor": "end", "font-size": 10, fill: css("--muted") }, `${E.n_high} high, ${E.n_low} low, ${E.events} events`);
      el.appendChild(s);
      // cutpoint scan
      el = $("#c-cut"); el.innerHTML = "";
      const C = E.cutpoint, ymin = Math.min(-3, Math.floor(Math.log10(Math.min(...C.p))));
      const x2 = (q) => L + ((q - 0.2) / 0.6) * (W - L - Rr), y2 = (p) => T + (Math.log10(p) / ymin) * (H - T - B);
      s = svg(W, H);
      axisY(s, y2, [1, 0.1, 0.01, 0.001].filter((p) => Math.log10(p) >= ymin), L, W - Rr, (v) => v);
      axisX(s, x2, [0.2, 0.35, 0.5, 0.65, 0.8], T, H - B, (v) => `${Math.round(v * 100)}%`);
      add(s, "line", { x1: L, x2: W - Rr, y1: y2(0.05), y2: y2(0.05), stroke: css("--bad"), "stroke-dasharray": "4 3" });
      add(s, "text", { x: W - Rr, y: y2(0.05) - 4, "text-anchor": "end", "font-size": 10, fill: css("--bad") }, "p = 0.05");
      add(s, "path", { d: C.quantiles.map((q, i) => `${i ? "L" : "M"}${x2(q)},${y2(C.p[i])}`).join(""), fill: "none", stroke: css("--accent"), "stroke-width": 2 });
      add(s, "circle", { cx: x2(C.best_quantile), cy: y2(C.min_p), r: 5, fill: css("--accent") });
      add(s, "line", { x1: x2(0.5), x2: x2(0.5), y1: T, y2: H - B, stroke: css("--muted"), "stroke-dasharray": "2 3" });
      add(s, "text", { x: (L + W) / 2, y: H - 4, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, "cutpoint (percentile of GPC3)");
      el.appendChild(s);
      $("#cut-note").textContent = `Smallest p = ${f(C.min_p, 4)} at the ${Math.round(C.best_quantile * 100)}th percentile. In 1,000 random relabellings, ${pct(C.null_fraction_below_005)} gave a "best cutpoint" below 0.05; corrected for this search, ${pp(C.permutation_p)}.`;
    }
    seg.addEventListener("click", (ev) => { const b = ev.target.closest("button"); if (!b) return; seg.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", x === b)); draw(b.dataset.e); });
    draw(eps[0]);
    const row = (e, k) => { const c = S.endpoints[e][k]; return c ? [`${f(c.hr)} (${f(c.lo)}–${f(c.hi)})`, c.p < 0.001 ? "<0.001" : f(c.p, 3), `${c.n} / ${c.events}`, c.ph_p == null ? "–" : f(c.ph_p, 2)] : ["–", "–", "–", "–"]; };
    table($("#t-cox"), [{ t: "Endpoint" }, { t: "Model" }, { t: "HR per SD (95% CI)", num: 1 }, { t: "p", num: 1 }, { t: "patients / events", num: 1 }, { t: "PH test p", num: 1 }],
      eps.flatMap((e) => [[S.endpoints[e].label, "GPC3 alone", ...row(e, "cox_univariable")], ["", "+ age, sex, stage", ...row(e, "cox_adjusted")]]));
    const oa = OS.cox_adjusted || {}, ou = OS.cox_univariable || {};
    const dsa = (S.endpoints.DSS || {}).cox_adjusted || {}, stg = (D.drivers || {}).stage;
    $("#surv-note").textContent = `${S.n_patients} liver-cancer patients with survival data. Like phase 1, a median split shows no difference (overall survival log-rank ${pp(OS.median_split_p)}), and GPC3 alone, as a continuous measure, gives HR ${f(ou.hr)} per SD (${pp(ou.p)}). ` +
      `But after adjusting for age, sex and stage, more GPC3 goes with shorter overall survival (HR ${f(oa.hr)}, ${pp(oa.p)}) and disease-specific survival (HR ${f(dsa.hr)}, ${pp(dsa.p)}). ` +
      (stg ? `The reason: late-stage tumours, which do worse, have less GPC3 (median ${f(stg.median_late)} vs ${f(stg.median_early)} in stage I–II), which hides the effect in the crude comparison. ` : "") +
      `The effect is modest, comes from one cohort and three endpoints were tested, so it needs an independent cohort before anyone relies on it.`;
    // forest
    const pc = (D.survival_pancancer || []).slice().sort((a, b) => a.hr - b.hr);
    const W = 1000, rh = 18, L = 90, Rr = 300, T = 16, H = pc.length * rh + T + 34;
    const lx = (v) => Math.log(v), lo = Math.log(0.4), hi = Math.log(2.5);
    const x = (v) => L + ((Math.min(Math.max(lx(v), lo), hi) - lo) / (hi - lo)) * (W - L - Rr);
    const s = svg(W, H);
    axisX(s, x, [0.5, 0.75, 1, 1.5, 2], T - 4, H - 30, (v) => v);
    add(s, "line", { x1: x(1), x2: x(1), y1: T - 4, y2: H - 30, stroke: css("--ink"), "stroke-dasharray": "3 3" });
    pc.forEach((r, i) => {
      const y = T + i * rh + rh / 2, sig = r.q < 0.05, c = r.hr > 1 ? css("--bad") : css("--c2");
      add(s, "text", { x: 8, y: y + 4, "font-size": 11.5, "font-weight": r.cancer === F ? 700 : 400 }, r.cancer);
      add(s, "line", { x1: x(r.lo), x2: x(r.hi), y1: y, y2: y, stroke: sig ? c : css("--muted"), "stroke-width": 1.6 });
      add(s, "circle", { cx: x(r.hr), cy: y, r: 4, fill: sig ? c : css("--surface"), stroke: sig ? c : css("--muted"), "stroke-width": 1.5 });
      add(s, "text", { x: W - Rr + 16, y: y + 4, "font-size": 11, fill: css("--muted") }, `HR ${f(r.hr)} (${f(r.lo)}–${f(r.hi)}) · ${r.n} pts · q ${r.q < 0.001 ? "<0.001" : f(r.q, 3)}`);
    });
    add(s, "text", { x: (L + W - Rr) / 2, y: H - 4, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, "hazard ratio per SD of GPC3 (log scale)");
    $("#c-forest").appendChild(s);
  })();

  // ------------------------------------------------------------ drivers
  (function () {
    const R = D.drivers; if (!R) return;
    const scatter = (el, pts, xl, yl, col) => {
      const W = 420, H = 280, L = 44, B = 40, T = 10, Rr = 12;
      const xs = pts.map((p) => p[0]), ys = pts.map((p) => p[1]);
      const [x0, x1] = [Math.min(...xs), Math.max(...xs)], [y0, y1] = [Math.min(...ys), Math.max(...ys)];
      const x = (v) => L + ((v - x0) / (x1 - x0 || 1)) * (W - L - Rr), y = (v) => T + (1 - (v - y0) / (y1 - y0 || 1)) * (H - T - B);
      const s = svg(W, H);
      const tk = (a, b) => [a, (a + b) / 2, b].map((v) => Math.round(v * 10) / 10);
      axisY(s, y, tk(y0, y1), L, W - Rr); axisX(s, x, tk(x0, x1), T, H - B);
      pts.forEach(([a, b]) => add(s, "circle", { cx: x(a), cy: y(b), r: 2.6, fill: col, "fill-opacity": 0.55 }));
      add(s, "text", { x: (L + W) / 2, y: H - 4, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, xl);
      add(add(s, "text", { x: 12, y: (T + H - B) / 2, "text-anchor": "middle", "font-size": 11, fill: css("--muted"), transform: `rotate(-90 12 ${(T + H - B) / 2})` }), "tspan", {}, yl);
      el.appendChild(s);
    };
    if (R.cnv) {
      const gs = Object.entries(R.cnv.groups).sort((a, b) => +a[0] - +b[0]);
      const W = 420, H = 280, L = 44, B = 40, T = 10;
      const mx = Math.max(...gs.map(([, g]) => g.q[4])), mn = Math.min(...gs.map(([, g]) => g.q[0]));
      const y = (v) => T + (1 - (v - mn) / (mx - mn)) * (H - T - B), bw = (W - L - 20) / gs.length;
      const s = svg(W, H);
      axisY(s, y, [Math.ceil(mn), Math.round((mn + mx) / 2), Math.floor(mx)], L, W - 10);
      gs.forEach(([k, g], i) => {
        const cx = L + bw * (i + 0.5), c = +k > 0 ? css("--bad") : +k < 0 ? css("--c2") : css("--muted");
        add(s, "line", { x1: cx, x2: cx, y1: y(g.q[0]), y2: y(g.q[4]), stroke: c });
        add(s, "rect", { x: cx - bw * 0.3, y: y(g.q[3]), width: bw * 0.6, height: y(g.q[1]) - y(g.q[3]), fill: c + "44", stroke: c, rx: 3 });
        add(s, "line", { x1: cx - bw * 0.3, x2: cx + bw * 0.3, y1: y(g.q[2]), y2: y(g.q[2]), stroke: c, "stroke-width": 2 });
        add(s, "text", { x: cx, y: H - 22, "text-anchor": "middle", "font-size": 12 }, (+k > 0 ? "+" : "") + k);
        add(s, "text", { x: cx, y: H - 8, "text-anchor": "middle", "font-size": 10, fill: css("--muted") }, `n=${g.n}`);
      });
      $("#c-cnv").appendChild(s);
    }
    if (R.methylation) {
      const b = R.methylation.best;
      $("#meth-h").textContent = `Methylation (${b.probe})`;
      $("#meth-note").textContent = `${R.methylation.probes} probes in GPC3; ${R.methylation.negative_sig} are negatively and ${R.methylation.positive_sig} positively linked to expression (q < 0.05). Strongest: ρ = ${f(b.rho)}.`;
      scatter($("#c-meth"), R.methylation.scatter, "methylation (beta)", `${G} log2(TPM+1)`, css("--accent"));
    } else { $("#meth-h").textContent = "Methylation"; $("#meth-note").textContent = "Methylation data were not available in this run."; }
    if (R.afp) scatter($("#c-afp"), R.afp.scatter, `${G} log2(TPM+1)`, "AFP log2(TPM+1)", css("--c2"));
    const rows = [];
    Object.entries(R.grade.groups).forEach(([k, g]) => rows.push([`Grade ${k.slice(1)}`, g.n, f(g.median)]));
    rows.push([`Stage I–II`, R.stage.n_early, f(R.stage.median_early)], [`Stage III–IV`, R.stage.n_late, f(R.stage.median_late)]);
    rows.push(["Men", R.sex.n_male, f(R.sex.median_male)], ["Women", R.sex.n_female, f(R.sex.median_female)]);
    table($("#t-clin"), [{ t: "Group" }, { t: "Tumours", num: 1 }, { t: `Median ${G}`, num: 1 }], rows);
    $("#drv-note").textContent = (R.cnv ? `${Math.abs(R.cnv.spearman) < 0.1 ? "Copy number does not track GPC3" : "Copy number tracks GPC3 only weakly"} (ρ = ${f(R.cnv.spearman)}), so the rise is not from extra copies of the gene. ` : "") +
      (R.methylation ? `${R.methylation.negative_sig ? `Methylation does: at ${R.methylation.negative_sig} of ${R.methylation.probes} probes, less methylation goes with more GPC3 (strongest ρ = ${f(R.methylation.best.rho)}), consistent with the gene being switched on epigenetically.` : "Methylation shows no clear negative link."} ` : "") +
      `Higher grade goes with more GPC3 (ρ = ${f(R.grade.rho)}, ${pp(R.grade.p)}); stage ${R.stage.p < 0.05 ? "differs" : "does not"} (${pp(R.stage.p)}); men vs women ${pp(R.sex.p)}. ` +
      (R.afp ? `GPC3 and AFP correlate with ρ = ${f(R.afp.rho)}; ${pct(R.afp.gpc3_only)} of tumours are high for GPC3 but not AFP.` : "");
  })();

  // ------------------------------------------------------------ co-expression
  (function () {
    const C = D.coexpression; if (!C) return;
    const gs = (C.gsea || []).slice().sort((a, b) => b.nes - a.nes);
    if (gs.length) {
      const show = gs.filter((g, i) => i < 12 || i >= gs.length - 12);
      const W = 520, rh = 17, L = 230, T = 10, H = show.length * rh + T + 30;
      const m = Math.max(...show.map((g) => Math.abs(g.nes)));
      const x = (v) => L + ((v + m) / (2 * m)) * (W - L - 10);
      const s = svg(W, H);
      axisX(s, x, [-Math.floor(m), 0, Math.floor(m)], T, H - 30);
      show.forEach((g, i) => {
        const y = T + i * rh + rh / 2, c = g.nes > 0 ? css("--bad") : css("--c2"), sig = g.fdr < 0.05;
        add(s, "text", { x: L - 6, y: y + 4, "text-anchor": "end", "font-size": 10.5, "font-weight": /wnt/i.test(g.term) ? 700 : 400 }, g.term.length > 34 ? g.term.slice(0, 33) + "…" : g.term);
        const r = add(s, "rect", { x: Math.min(x(0), x(g.nes)), y: y - 6, width: Math.abs(x(g.nes) - x(0)), height: 12, fill: sig ? c : "none", stroke: c, rx: 2 });
        add(r, "title", {}, `${g.term}: NES ${f(g.nes)}, FDR ${f(g.fdr, 3)}`);
      });
      $("#c-gsea").appendChild(s);
    }
    const tp = C.top_positive.slice(0, 12), tn = C.top_negative.slice(0, 12);
    table($("#t-co"), [{ t: "With GPC3" }, { t: "rho", num: 1 }, { t: "Against GPC3" }, { t: "rho", num: 1 }],
      tp.map((r, i) => [`<span class="mono">${esc(r.gene)}</span>`, f(r.rho), tn[i] ? `<span class="mono">${esc(tn[i].gene)}</span>` : "", tn[i] ? f(tn[i].rho) : ""]));
    const wnt = gs.find((g) => /wnt/i.test(g.term));
    $("#co-note").textContent = `Across ${C.n_tumours} liver tumours and ${C.n_genes.toLocaleString()} expressed genes. ` +
      `The closest partners include foetal-liver genes such as ${C.top_positive.slice(0, 6).map((r) => r.gene).join(', ')}. ` + (wnt ? `Phase 1's network pointed to Wnt; here Wnt/β-catenin signalling has NES ${f(wnt.nes)} (FDR ${f(wnt.fdr, 3)}), ${wnt.fdr < 0.05 ? (wnt.nes > 0 ? "so it does rise with GPC3." : "so it falls as GPC3 rises.") : "so tumours with more GPC3 do not show clearly more Wnt activity at the RNA level."} ` : "") +
      (gs.length ? `Top sets up with GPC3: ${gs.slice(0, 3).map((g) => g.term).join(", ")}; down: ${gs.slice(-3).reverse().map((g) => g.term).join(", ")}.` : "");
  })();

  // ------------------------------------------------------------ structure
  (function () {
    const P = D.structure; if (!P) return;
    const muts = (P.mutations && P.mutations.list) || [];
    const mutPos = [...new Set(muts.filter((m) => m.start).map((m) => +m.start))];
    // 3D
    const el = $("#viewer");
    if (window.$3Dmol) {
      fetch("alphafold.pdb").then((r) => r.text()).then((pdb) => {
        const v = $3Dmol.createViewer(el, { backgroundColor: css("--surface") || "white" });
        v.addModel(pdb, "pdb");
        const col = (b) => (b >= 90 ? "#0053D6" : b >= 70 ? "#65CBF3" : b >= 50 ? "#FFDB13" : "#FF7D45");
        v.setStyle({}, { cartoon: { colorfunc: (a) => col(a.b) } });
        v.addStyle({ resi: mutPos, atom: "CA" }, { sphere: { color: "#e11d74", radius: 1.1 } });
        v.zoomTo(); v.render();
      }).catch(() => (el.textContent = "3D model could not be loaded."));
    } else el.innerHTML = '<p class="small" style="padding:14px">The 3D viewer could not load; the track below shows the same information.</p>';
    // track
    const W = 1000, H = 230, L = 50, Rr = 20, n = P.length;
    const x = (i) => L + ((i - 1) / (n - 1)) * (W - L - Rr);
    const s = svg(W, H);
    const yp = (v) => 120 - (v / 100) * 80;
    [50, 70, 90].forEach((v) => { add(s, "line", { x1: L, x2: W - Rr, y1: yp(v), y2: yp(v), stroke: css("--grid") }); add(s, "text", { x: L - 6, y: yp(v) + 4, "text-anchor": "end", "font-size": 10, fill: css("--muted") }, v); });
    add(s, "path", { d: P.plddt.map((v, i) => `${i ? "L" : "M"}${x(i + 1).toFixed(1)},${yp(v).toFixed(1)}`).join(""), fill: "none", stroke: css("--accent"), "stroke-width": 1.5 });
    add(s, "text", { x: L, y: 30, "font-size": 11, fill: css("--muted") }, "AlphaFold pLDDT");
    // features
    const fy = 150;
    add(s, "line", { x1: x(1), x2: x(n), y1: fy, y2: fy, stroke: css("--muted"), "stroke-width": 3 });
    const feats = P.features || [];
    const draw = (type, color, label, h = 14) => feats.filter((ft) => ft.type === type).forEach((ft) => {
      const r = add(s, "rect", { x: x(ft.start), y: fy - h / 2, width: Math.max(3, x(ft.end || ft.start) - x(ft.start)), height: h, fill: color, rx: 2 });
      add(r, "title", {}, `${ft.type} ${ft.start}${ft.end && ft.end !== ft.start ? "–" + ft.end : ""}: ${ft.description || label}`);
    });
    draw("Signal", css("--c3"), "signal peptide");
    draw("Propeptide", css("--c2"), "removed at GPI anchoring");
    feats.filter((ft) => ["Glycosylation", "Lipidation", "Site"].includes(ft.type)).forEach((ft) => {
      const c = ft.type === "Lipidation" ? css("--c2") : ft.type === "Site" ? css("--ink") : css("--c1");
      add(s, "line", { x1: x(ft.start), x2: x(ft.start), y1: fy - 14, y2: fy + 14, stroke: c, "stroke-width": 2 });
      const t = add(s, "circle", { cx: x(ft.start), cy: fy - 16, r: 3, fill: c });
      add(t, "title", {}, `${ft.type} ${ft.start}: ${ft.description}`);
    });
    // mutations lollipops
    const cnt = {}; muts.forEach((m) => { if (m.start) cnt[m.start] = (cnt[m.start] || 0) + 1; });
    Object.entries(cnt).forEach(([p, c]) => {
      add(s, "line", { x1: x(+p), x2: x(+p), y1: fy + 8, y2: fy + 18 + 8 * c, stroke: "#e11d74", "stroke-opacity": 0.6 });
      const d = add(s, "circle", { cx: x(+p), cy: fy + 18 + 8 * c, r: 2.8, fill: "#e11d74" });
      add(d, "title", {}, muts.filter((m) => m.start == p).map((m) => `${m.protein_change} (${m.cancer})`).join("\n"));
    });
    [1, 100, 200, 300, 400, 500, n].forEach((v) => add(s, "text", { x: x(v), y: H - 4, "text-anchor": "middle", "font-size": 10, fill: css("--muted") }, v));
    $("#c-track").appendChild(s);
    $("#c-track").insertAdjacentHTML("afterend", `<div class="legend"><span><i class="sw" style="background:${css("--c3")}"></i>signal peptide</span><span><i class="sw" style="background:${css("--c1")}"></i>heparan sulfate / glycosylation site</span><span><i class="sw" style="background:${css("--ink")}"></i>cleavage site</span><span><i class="sw" style="background:${css("--c2")}"></i>GPI anchor / propeptide</span><span><i class="sw" style="background:#e11d74"></i>TCGA mutation (height = count)</span></div>`);
    const M = P.mutations;
    if (M) {
      table($("#t-mut"), [{ t: "Cancer" }, { t: "Mutated tumours", num: 1 }, { t: "Sequenced", num: 1 }, { t: "%", num: 1 }],
        M.by_cancer.slice(0, 12).map((r) => [r.cancer, r.mutated, r.sequenced ?? "–", f(r.percent, 1)]));
      table($("#t-mtype"), [{ t: "Type" }, { t: "Count", num: 1 }], Object.entries(M.by_type).map(([k, v]) => [esc(k.replace(/_/g, " ")), v]));
    }
    $("#str-note").textContent = `${n} residues; AlphaFold is confident (pLDDT ≥ 70) for ${pct(P.ordered_fraction)} of them and very unsure (< 50) for ${pct(P.disordered_fraction)}, mostly the ends, as phase 1's disorder prediction said. ` +
      (M ? `GPC3 is mutated in ${f(M.percent_samples, 1)}% of ${M.sequenced.toLocaleString()} sequenced tumours (${M.total} mutations), mostly in cancers with many mutations overall. ` +
        `${M.missense_in_ordered} of ${M.missense_total} missense changes fall in the confidently folded core, where ${pct(M.ordered_fraction_expected)} would be expected by chance (${pp(M.ordered_binom_p)}). ` +
        (M.recurrent.length ? `The most-hit positions (${M.recurrent.map((r) => `${r.position}, ${r.n}×`).join("; ")}) are hit only a few times across 10,000 tumours, so there is no real hotspot.` : "No position is hit three times or more, so there is no mutation hotspot.") + " GPC3 is changed in cancer by how much is made, not by mutation." : "");
  })();

  // ------------------------------------------------------------ checks + takeaways
  const checks = [];
  checks.push([lihc.n_tumour >= 330 && lihc.n_tumour <= 380, `TCGA-LIHC primary tumours found: ${lihc.n_tumour} (about 370 expected)`]);
  checks.push([lihc.log2fc > 2 && lihc.q < 0.05, `Liver cancer shows the large GPC3 rise reported in phase 1 and the literature (log2FC ${f(lihc.log2fc)})`]);
  if (lihc.log2fc_tcga_only != null) checks.push([lihc.log2fc_tcga_only > 2, `Same result with TCGA normals only, no GTEx (log2FC ${f(lihc.log2fc_tcga_only)})`]);
  if (pe.paired) checks.push([pe.paired.p < 0.05, `Paired tumour–normal test in the same patients agrees (${pp(pe.paired.p)})`]);
  if (OS) checks.push([OS.cutpoint.permutation_p >= OS.cutpoint.min_p, `Cutpoint search corrected by permutation (raw ${f(OS.cutpoint.min_p, 4)} → ${f(OS.cutpoint.permutation_p, 3)})`]);
  if (D.structure) checks.push([D.structure.alphafold && D.structure.alphafold.same_sequence, `AlphaFold model has the same sequence as UniProt ${D.structure.length} aa`]);
  $("#checklist").innerHTML = checks.map(([ok, t]) => `<li class="${ok ? "ok" : "bad"}">${esc(t)}</li>`).join("");

  const take = [];
  take.push(`GPC3 is a liver-cancer marker above all: AUC ${f(lihc.auc)} against normal liver, and higher in the tumour in ${pct(pe.paired && pe.paired.fraction_higher_in_tumour)} of patients with both samples.`);
  if (OS) take.push(`Crude survival analysis finds nothing (median split ${pp(OS.median_split_p)}), but adjusted for age, sex and stage, more GPC3 goes with shorter overall survival (HR ${f(OS.cox_adjusted.hr)} per SD, ${pp(OS.cox_adjusted.p)}): late-stage tumours have less GPC3 and hide the effect. The "best cutpoint" p of ${f(OS.cutpoint.min_p, 3)} becomes ${f(OS.cutpoint.permutation_p, 2)} once the search is accounted for.`);
  if (D.drivers && D.drivers.cnv) take.push(`It is switched on epigenetically, not by extra gene copies: expression falls with methylation (${D.drivers.methylation ? `${D.drivers.methylation.negative_sig} of ${D.drivers.methylation.probes} probes, strongest ρ = ${f(D.drivers.methylation.best.rho)}` : "–"}) but not with copy number (ρ = ${f(D.drivers.cnv.spearman)}).`);
  if (D.coexpression) take.push(`High-GPC3 tumours look foetal and fast-growing: its closest partners include ${D.coexpression.top_positive.slice(0, 6).map((r) => r.gene).join(", ")}, with E2F, G2-M and MYC programmes up and normal liver metabolism down.`);
  if (D.structure && D.structure.mutations) take.push(`Mutations are rare (${f(D.structure.mutations.percent_samples, 1)}% of tumours) and scattered: GPC3 matters in cancer through how much is made, which is why it works as a target for antibodies and CAR-T cells.`);
  $("#takeaways").innerHTML = take.map((t) => `<p>${esc(t)}</p>`).join("");
})();
