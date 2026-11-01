"use strict";

const $ = (id) => document.getElementById(id);
const COLORS = ["var(--s1)", "var(--s2)", "var(--s3)", "var(--s4)", "var(--s5)", "var(--s6)", "var(--s7)", "var(--s8)"];
const plot = new LinePlot($("plot"));
let ST = null, VIEW = "mass";

function phaseNames(steps) {
  const seen = [];
  steps.forEach((s) => s.solids.forEach((p) => { if (!seen.includes(p.name)) seen.push(p.name); }));
  return seen;
}
function label(steps, name) {
  for (const s of steps) for (const p of s.solids) if (p.name === name) return p.label;
  return name;
}
function grams(step, name) {
  const p = step.solids.find((q) => q.name === name);
  return p ? p.grams : 0;
}
function color(cssVar) { return getComputedStyle(document.body).getPropertyValue(cssVar.slice(4, -1)).trim() || cssVar; }

function series(steps) {
  const x = steps.map((s) => s.t);
  if (VIEW === "mass") {
    const out = [{ label: "Melt", color: color("var(--s1)"), y: steps.map((s) => s.liquid_grams) }];
    phaseNames(steps).forEach((n, i) => out.push({ label: label(steps, n), color: color(COLORS[(i + 1) % 8]), y: steps.map((s) => grams(s, n)) }));
    return { x, s: out, opts: { yLabel: "grams per 100 g of starting melt" } };
  }
  if (VIEW === "liquid") {
    const ox = Object.keys(steps[0].liquid).filter((o) => o !== "H2O");
    return { x, s: ox.map((o, i) => ({ label: o, color: color(COLORS[i % 8]), dash: i >= 8 ? "5 4" : null, y: steps.map((s) => s.liquid[o] ?? null) })),
      opts: { yLabel: "melt, wt%" } };
  }
  const out = [];
  const num = (name, key) => steps.map((s) => { const p = s.solids.find((q) => q.name === name); return p && p.numbers && p.numbers[key] !== undefined && p.grams > 1e-9 ? p.numbers[key] : null; });
  const add = (name, key, lab, col, dash) => { const y = num(name, key); if (y.some((v) => v !== null)) out.push({ label: lab, color: color(col), dash, y }); };
  add("olivine", "Fo", "Olivine Fo", "var(--s3)");
  add("feldspar", "An", "Plagioclase An", "var(--s2)");
  add("clinopyroxene", "Mg#", "Clinopyroxene Mg#", "var(--s4)");
  add("clinopyroxene", "Wo", "Clinopyroxene Wo", "var(--s5)", "5 4");
  out.push({ label: "Melt Mg# ", color: color("var(--s6)"), dash: "5 4", y: steps.map((s) => s.mg_number) });
  return { x, s: out, opts: { yLabel: "mol%" } };
}

function draw() {
  if (!ST || !ST.steps.length) return;
  const { x, s, opts } = series(ST.steps);
  plot.set(x, s, Object.assign({ xLabel: "temperature, °C", xName: "", xUnit: " °C", hoverFixed: 2 }, opts));
  $("legend").replaceChildren(...s.map((q) => h("span", {}, h("span", { class: "sw", style: `border-color:${q.color}${q.dash ? ";border-top-style:dashed" : ""}` }), q.label)));
  $("axisnote").textContent = ST.inputs.mode === "fractional" && VIEW === "mass" ? "Crystal masses are the totals removed so far." : "";
}

function story(st) {
  const steps = st.steps, out = [];
  const first = steps.find((s) => s.solids.some((p) => p.grams > 1e-6));
  const t0 = steps[0].t, last = steps[steps.length - 1];
  if (!first) out.push(`No mineral crystallizes between ${t0.toFixed(0)} and ${last.t.toFixed(0)} °C; the melt stays fully liquid.`);
  else {
    out.push(`The first crystals appear at ${first.t.toFixed(0)} °C: ${first.solids.map((p) => `${p.label.toLowerCase()} (${p.formula})`).join(" and ")}.`);
    const order = [];
    steps.forEach((s) => s.solids.forEach((p) => { if (p.grams > 1e-6 && !order.find((o) => o.n === p.name)) order.push({ n: p.name, t: s.t, l: p.label }); }));
    if (order.length > 1) out.push(`Order of crystallization: ${order.map((o) => `${o.l.toLowerCase()} at ${o.t.toFixed(0)} °C`).join(", then ")}.`);
    out.push(`At ${last.t.toFixed(0)} °C the system is ${last.melt_pct.toFixed(1)}% melt by mass` +
      (st.inputs.mode === "fractional" ? ` (${last.liquid_grams.toFixed(1)} g of melt left from 100 g).` : ".") +
      ` The melt has ${last.liquid.SiO2.toFixed(1)} wt% SiO2 and Mg# ${last.mg_number.toFixed(0)}.`);
  }
  if (first && !(st.inputs.phases || []).includes("clinopyroxene")) out.push("Clinopyroxene was not allowed in this run.");
  return out;
}

function table(steps) {
  const names = phaseNames(steps);
  const ox = Object.keys(steps[0].liquid);
  const head = h("tr", {}, h("th", {}, "T, °C"), h("th", {}, "Melt, %"), ...names.map((n) => h("th", { class: "l" }, `${label(steps, n)}, g`)), ...ox.map((o) => h("th", {}, o)), h("th", {}, "Mg#"));
  const rows = steps.map((s) => h("tr", {}, h("td", {}, s.t.toFixed(1)), h("td", {}, s.melt_pct.toFixed(2)),
    ...names.map((n) => { const p = s.solids.find((q) => q.name === n); return h("td", { class: "l" }, p ? `${p.grams.toFixed(3)} ${p.formula}` : ""); }),
    ...ox.map((o) => h("td", {}, (s.liquid[o] ?? 0).toFixed(2))), h("td", {}, s.mg_number.toFixed(1))));
  $("table").replaceChildren(h("thead", {}, head), h("tbody", {}, rows));
}

function csv() {
  const steps = ST.steps, names = phaseNames(steps), ox = Object.keys(steps[0].liquid);
  const lines = [["T_C", "P_bar", "melt_wt_pct", "melt_g", ...names.flatMap((n) => [`${n}_g`, `${n}_composition`]), ...ox.map((o) => `melt_${o}`), "melt_Mg#"].join(",")];
  for (const s of steps) {
    lines.push([s.t, s.p, s.melt_pct.toFixed(4), s.liquid_grams.toFixed(4),
      ...names.flatMap((n) => { const p = s.solids.find((q) => q.name === n); return p ? [p.grams.toFixed(5), p.formula] : ["0", ""]; }),
      ...ox.map((o) => (s.liquid[o] ?? 0).toFixed(4)), s.mg_number.toFixed(2)].join(","));
  }
  return lines.join("\r\n");
}
function download(name, text, type) {
  const a = h("a", { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.append(a); a.click(); a.remove();
}

function show(st) {
  ST = st;
  const has = st.steps.length > 0;
  $("empty").style.display = has ? "none" : "";
  $("content").style.display = has ? "flex" : "none";
  const pill = $("pill-run");
  pill.style.display = st.status === "idle" ? "none" : "";
  pill.textContent = { running: `running ${st.progress[0]}/${st.progress[1]}`, done: "done", stopped: "stopped", error: "error" }[st.status] || "";
  pill.className = `pill ${st.status === "error" ? "bad" : st.status === "running" ? "warn" : "ok"}`;
  if (!has) return;
  const i = st.inputs;
  $("title").textContent = `${i.name || "Melt"}: ${i.mode} crystallization`;
  $("meta").textContent = `${i.t_start} to ${i.t_end} °C in ${i.step} °C steps, ${i.p_bar} bar, oxygen ${i.fo2 === "closed" ? "closed system" : i.fo2}`;
  $("story").replaceChildren(...story(st).map((t) => h("p", {}, t)), ...(st.status === "error" ? [h("div", { class: "errbox" }, st.message)] : []));
  $("notes").replaceChildren(...st.notes.map((n) => h("li", {}, n)));
  table(st.steps);
  draw();
}

document.querySelectorAll("#view button").forEach((b) => b.addEventListener("click", () => {
  VIEW = b.dataset.view;
  document.querySelectorAll("#view button").forEach((x) => x.classList.toggle("on", x === b));
  draw();
}));
$("dl-csv").addEventListener("click", () => download("magma_run.csv", csv(), "text/csv"));
$("dl-svg").addEventListener("click", () => download("magma_plot.svg", plot.svgText(), "image/svg+xml"));
connect(show, (ok) => { $("pill-link").textContent = ok ? "connected" : "not connected"; $("pill-link").className = `pill ${ok ? "ok" : "bad"}`; });
