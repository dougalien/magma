"use strict";

const $ = (id) => document.getElementById(id);
const COLORS = ["var(--s1)", "var(--s2)", "var(--s3)", "var(--s4)", "var(--s5)", "var(--s6)", "var(--s7)", "var(--s8)"];
const plot = new LinePlot($("plot"));
let ST = null, VIEW = "mass", INFO = null;
fetch("/api/presets").then((r) => r.json()).then((j) => { INFO = j; draw(); });

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
    if (steps.some((s) => s.pooled_grams > 0)) out.push({ label: "Melt removed (total)", color: color("var(--s1)"), dash: "5 4", y: steps.map((s) => s.pooled_grams || 0) });
    phaseNames(steps).forEach((n, i) => out.push({ label: label(steps, n), color: color(COLORS[(i + 1) % 8]), y: steps.map((s) => grams(s, n)) }));
    return { x, s: out, opts: { yLabel: "grams per 100 g of starting material" } };
  }
  if (VIEW === "tas") return tasSeries(steps);
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
  add("alkali-feldspar", "Or", "K-feldspar Or", "var(--s5)");
  add("hornblende", "Mg#", "Hornblende Mg#", "var(--s3)", "5 4");
  add("biotite", "Mg#", "Biotite Mg#", "var(--s4)", "5 4");
  add("clinopyroxene", "Mg#", "Augite Mg#", "var(--s4)");
  add("clinopyroxene", "Wo", "Augite Wo", "var(--s5)", "5 4");
  add("pigeonite", "Mg#", "Pigeonite Mg#", "var(--s8)");
  add("orthopyroxene", "Mg#", "Orthopyroxene Mg#", "var(--s7)");
  add("spinel", "Cr#", "Spinel Cr#", "var(--s8)");
  add("spinel", "Usp", "Spinel ulvöspinel", "var(--s1)", "2 3");
  add("rhombohedral", "Ilm", "Ilmenite (FeTiO3)", "var(--s7)", "2 3");
  out.push({ label: "Melt Mg# ", color: color("var(--s6)"), dash: "5 4", y: steps.map((s) => s.mg_number) });
  return { x, s: out, opts: { yLabel: "mol%" } };
}

/* The melt's path on the total-alkali-silica diagram (Le Bas et al. 1986), with the field boundaries. */
function tasSeries(steps) {
  const x = steps.map((s) => s.tas_xy[0]);
  const out = [{ label: "Melt", color: color("var(--s2)"), y: steps.map((s) => s.tas_xy[1]) }];
  const faint = getComputedStyle(document.body).getPropertyValue("--faint").trim() || "#888";
  if (INFO) INFO.tas.lines.forEach((ln) => out.push({ label: "", noHover: true, noLegend: true, color: faint, dash: "none",
    x: ln.map((p) => p[0]), y: ln.map((p) => p[1]) }));
  const labels = INFO ? Object.entries(INFO.tas.fields).map(([name, pts]) => {
    const n = pts.length;
    return { text: name.replace(" or ", "/"), x: pts.reduce((a, p) => a + p[0], 0) / n, y: pts.reduce((a, p) => a + p[1], 0) / n };
  }) : [];
  const marks = [];
  steps.forEach((s, i) => { if (i === 0 || s.magma_type !== steps[i - 1].magma_type) marks.push({ x: s.tas_xy[0], y: s.tas_xy[1], text: `${s.t.toFixed(0)} °C`, color: color("var(--s2)") }); });
  if (ST.source && ST.source.primary) {
    const pr = ST.source.tas_xy;
    marks.push({ x: pr[0], y: pr[1], text: "primary melt", color: color("var(--s1)") });
  }
  return { x, s: out, opts: { xLabel: "SiO2, wt% (volatile-free)", yLabel: "Na2O + K2O, wt%", xName: "SiO2", xUnit: " wt%",
    xRange: [40, 80], yRange: [0, 16], labels, marks, sub: (i) => `${steps[i].t.toFixed(0)} °C, ${steps[i].liquid_grams.toFixed(1)} g melt: ${steps[i].tas}` } };
}

/* Bowen's reaction series, with what this run grew and the temperature each mineral appeared at. */
function bowenView(steps) {
  const first = {};
  steps.forEach((s) => s.solids.forEach((p) => {
    const g = p.step_grams !== undefined ? p.step_grams : p.grams;
    if (g > 1e-6 && !(p.name in first)) first[p.name] = { t: s.t, label: p.label };
  }));
  const plag = steps.map((s) => ({ t: s.t, p: s.solids.find((q) => q.name === "feldspar" && (q.step_grams !== undefined ? q.step_grams : q.grams) > 1e-6) }))
    .filter((o) => o.p && o.p.numbers && o.p.numbers.An !== undefined);
  const box = (title, names) => {
    const seen = names.filter((n) => n in first).sort((a, b) => first[b].t - first[a].t);
    if (!names.length) return h("div", { class: "box off" }, h("b", {}, title), "not in Magma's model yet");
    if (!seen.length) return h("div", { class: "box off" }, h("b", {}, title), "did not form in this run");
    return h("div", { class: "box on" }, h("b", {}, title), seen.map((n) => `${first[n].label.toLowerCase()} from ${first[n].t.toFixed(0)} °C`).join("; "));
  };
  const B = INFO ? INFO.bowen : { discontinuous: [], continuous: [], late: [] };
  const disc = [];
  B.discontinuous.forEach(([title, names], i) => { if (i) disc.push(h("div", { class: "arrow" }, "↓")); disc.push(box(title, names)); });
  let cont;
  if (!plag.length) cont = [box("Ca-rich plagioclase", ["feldspar"]), h("div", { class: "arrow" }, "↓"), h("div", { class: "box off" }, h("b", {}, "Na-rich plagioclase"), "did not form in this run")];
  else {
    const a = plag[0], b = plag[plag.length - 1];
    cont = [h("div", { class: "box on" }, h("b", {}, "Ca-rich plagioclase"), `An${a.p.numbers.An.toFixed(0)} from ${a.t.toFixed(0)} °C`),
      h("div", { class: "arrow" }, "↓ continuous: the same mineral, less calcic"),
      b.p.numbers.An < 50 ? h("div", { class: "box on" }, h("b", {}, "Na-rich plagioclase"), `An${b.p.numbers.An.toFixed(0)} at ${b.t.toFixed(0)} °C`)
        : h("div", { class: "box off" }, h("b", {}, "Na-rich plagioclase"), `still An${b.p.numbers.An.toFixed(0)} at ${b.t.toFixed(0)} °C`)];
  }
  const late = B.late.map(([title, names]) => box(title, names));
  const oxides = ["spinel", "titanomagnetite", "rhombohedral", "garnet", "water"].filter((n) => n in first);
  const t0 = steps[0].t, t1 = steps[steps.length - 1].t;
  return [h("div", { class: "bowen" },
    h("div", { class: "col" }, h("h4", {}, "Discontinuous series"), ...disc),
    h("div", { class: "col" }, h("h4", {}, "Continuous series"), ...cont),
    h("div", { class: "late" }, ...late),
    h("div", { class: "temp" }, h("span", {}, `↑ hot: ${t0.toFixed(0)} °C`), h("span", {}, `cool: ${t1.toFixed(0)} °C ↓`))),
  h("p", { class: "hint", style: "text-align:center" }, oxides.length
    ? `Also formed, outside Bowen's scheme: ${oxides.map((n) => `${first[n].label.toLowerCase()} from ${first[n].t.toFixed(0)} °C`).join(", ")}.`
    : "Dashed boxes: minerals that did not form, or that Magma cannot make yet.")];
}

function draw() {
  if (!ST || !ST.steps.length) return;
  const bowen = VIEW === "bowen";
  $("plot").style.display = bowen ? "none" : "";
  $("bowen").style.display = bowen ? "" : "none";
  if (bowen) {
    $("bowen").replaceChildren(...bowenView(ST.steps));
    $("legend").replaceChildren();
    $("axisnote").textContent = "";
    return;
  }
  const { x, s, opts } = series(ST.steps);
  plot.set(x, s, Object.assign({ xLabel: "temperature, °C", xName: "", xUnit: " °C", hoverFixed: 2 }, opts));
  $("legend").replaceChildren(...s.filter((q) => !q.noLegend).map((q) => h("span", {}, h("span", { class: "sw", style: `border-color:${q.color}${q.dash ? ";border-top-style:dashed" : ""}` }), q.label)));
  const melting = ST.inputs.process === "melt";
  $("axisnote").textContent = VIEW === "tas" ? "Total alkali-silica diagram (Le Bas et al. 1986); dots mark a new magma type."
    : ST.inputs.mode === "fractional" && VIEW === "mass"
    ? (melting ? "Mineral masses are what is left in the residue; the melt left in it is 1%."
      : ST.inputs.process === "example" ? "Grams per 100 g of primary melt; crystal masses are the totals removed so far." : "Crystal masses are the totals removed so far.") : "";
}

function meltStory(st) {
  const steps = st.steps, out = [], first = steps[0], last = steps[steps.length - 1];
  const frac = st.inputs.mode === "fractional";
  out.push(`At ${first.t.toFixed(0)} °C the rock is ${first.melt_pct.toFixed(1)}% melt; its residue is ${first.solids.map((p) => p.label.toLowerCase()).join(", ")}.`);
  const gone = [];
  for (let i = 1; i < steps.length; i++) steps[i - 1].solids.forEach((p) => {
    if (p.grams > 1e-3 && !steps[i].solids.some((q) => q.name === p.name && q.grams > 1e-6) && !gone.find((g) => g.n === p.name)) gone.push({ n: p.name, l: p.label, t: steps[i].t });
  });
  if (gone.length) out.push(`Used up: ${gone.map((g) => `${g.l.toLowerCase()} by ${g.t.toFixed(0)} °C`).join(", ")}.`);
  if (frac) {
    const pm = last.pooled;
    out.push(`By ${last.t.toFixed(0)} °C, ${last.pooled_grams.toFixed(1)} g of melt has been removed from 100 g of rock` +
      (pm ? `; the pooled melt has ${pm.SiO2.toFixed(1)} wt% SiO2, ${pm.MgO.toFixed(1)} wt% MgO and Mg# ${last.pooled_mg.toFixed(0)}.` : "."));
  } else {
    out.push(`At ${last.t.toFixed(0)} °C the rock is ${last.melt_pct.toFixed(1)}% melt, with ${last.liquid.SiO2.toFixed(1)} wt% SiO2, ${(last.liquid.MgO ?? 0).toFixed(1)} wt% MgO and Mg# ${last.mg_number.toFixed(0)}.`);
  }
  return out;
}

/* The magma types the melt passes through: the first step of each. */
function typeSteps(steps) {
  return steps.filter((s, i) => i === 0 || s.magma_type !== steps[i - 1].magma_type);
}

/* The step after which the melt's silica falls while its alkalis climb past 10 wt%: the sodium-rich
   runaway of the MELTS liquid (in ThermoEngine's MELTS too), not a real trend.  -1 if none. */
function runawayIndex(steps) {
  let best = -Infinity;
  for (let i = 0; i < steps.length; i++) {
    const [si, alk] = steps[i].tas_xy;
    if (si < best - 1.5 && alk > 10) return i;
    best = Math.max(best, si);
  }
  return -1;
}
function runawayNote(steps) {
  const k = runawayIndex(steps);
  if (k < 0) return null;
  const s = steps[k];
  return `Below about ${s.t.toFixed(0)} °C the melt's silica falls while its Na2O + K2O climbs (${s.tas_xy[1].toFixed(1)} wt% here). ` +
    "That is the MELTS liquid model holding too much sodium, which ThermoEngine's own MELTS does as well; it is not a real trend, so read the path no further than this.";
}

function exampleStory(st) {
  const src = st.source, steps = st.steps, i = st.inputs, ex = i.example || {}, out = [];
  if (!src) { out.push(st.status === "running" ? (st.message || "Stage 1: melting the mantle rock...") : "Stage 1 did not finish."); return out; }
  const pr = src.primary;
  out.push(`Stage 1, melting the mantle: at ${src.t.toFixed(0)} °C and ${(src.p / 1000).toFixed(0)} kbar the ${i.name || "rock"} is ` +
    `${src.melt_pct.toFixed(1)}% melt, and the solid left behind is ${src.solids.map((p) => p.label.toLowerCase()).join(", ")}. ` +
    `That melt is the primary magma: ${pr.SiO2.toFixed(1)} wt% SiO2, ${(pr.MgO ?? 0).toFixed(1)} wt% MgO, Mg# ${src.primary_mg.toFixed(0)}, a ${src.magma_type} (TAS: ${src.tas}).`);
  if (!steps.length) return out;
  const last = steps[steps.length - 1];
  const water = ex.h2o > 0 ? `with ${ex.h2o} wt% water added` : "dry";
  out.push(`Stage 2, cooling in the crust: the melt rises to ${(i.p_bar / 1000).toFixed(1)} kbar, ${water}, at the ${i.fo2} oxygen buffer, ` +
    "and the crystals are taken out as they form (fractional crystallization).");
  const order = [];
  steps.forEach((s) => s.solids.forEach((p) => { const g = p.step_grams !== undefined ? p.step_grams : p.grams; if (g > 1e-6 && !order.find((o) => o.n === p.name)) order.push({ n: p.name, t: s.t, l: p.label }); }));
  if (order.length) out.push(`Minerals in the order they appear: ${order.map((o) => `${o.l.toLowerCase()} (${o.t.toFixed(0)} °C)`).join(", ")}. ` +
    "See Bowen's series for how this compares with Bowen's order.");
  const ts = typeSteps(steps);
  if (ts.length > 1) out.push(`The melt changes from ${ts.map((s, k) => k === 0 ? s.magma_type : `${s.magma_type} at ${s.t.toFixed(0)} °C (${s.liquid_grams.toFixed(0)} g of melt left)`).join(", to ")}.`);
  out.push(`At ${last.t.toFixed(0)} °C, ${last.liquid_grams.toFixed(1)} g of melt is left from 100 g of primary magma: ${last.tas_xy[0].toFixed(1)} wt% SiO2 and ` +
    `${last.tas_xy[1].toFixed(1)} wt% Na2O + K2O, a ${last.magma_type} by silica (TAS: ${last.tas}).`);
  const rn = runawayNote(steps);
  if (rn) out.push(rn);
  if (/trachy|phono|tephri|foid/.test(last.tas)) out.push("The TAS name is alkaline because the model's liquid keeps too much sodium (a known failing of MELTS as well); a real melt with this silica would be the subalkaline type named by silica.");
  const grew = new Set(order.map((o) => o.n));
  const late = ["hornblende", "biotite", "alkali-feldspar", "muscovite", "quartz"].filter((n) => !grew.has(n));
  if (late.length) out.push(`Late minerals of Bowen's series that did not form in this run: ${late.map((n) => ({ "alkali-feldspar": "K-feldspar" }[n] || n)).join(", ")}. ` +
    "Magma has them (as MELTS does), but MELTS's hornblende and biotite are simple models that seldom saturate, and muscovite needs higher pressure (above about 3 kbar).");
  if (last.magma_type !== "rhyolite") out.push("To follow Bowen's series to the end, run the \"Arc andesite\" or \"Dacite\" card: a crustal andesite, richer in potassium, goes on to rhyolite, " +
    "then grows quartz, biotite and K-feldspar together at about 720 °C, the granite minimum.");
  return out;
}

function story(st) {
  if (st.inputs.process === "melt") return meltStory(st);
  if (st.inputs.process === "example") return exampleStory(st);
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
  const rn = runawayNote(steps);
  if (rn) out.push(rn);
  return out;
}

function table(steps) {
  const names = phaseNames(steps);
  const ox = Object.keys(steps[0].liquid);
  const head = h("tr", {}, h("th", {}, "T, °C"), h("th", {}, "Melt, %"), h("th", { class: "l" }, "Magma type"), ...names.map((n) => h("th", { class: "l" }, `${label(steps, n)}, g`)), ...ox.map((o) => h("th", {}, o)), h("th", {}, "Mg#"));
  const rows = steps.map((s) => h("tr", {}, h("td", {}, s.t.toFixed(1)), h("td", {}, s.melt_pct.toFixed(2)), h("td", { class: "l", title: `TAS: ${s.tas}` }, s.magma_type),
    ...names.map((n) => { const p = s.solids.find((q) => q.name === n); return h("td", { class: "l" }, p ? `${p.grams.toFixed(3)} ${p.formula}` : ""); }),
    ...ox.map((o) => h("td", {}, (s.liquid[o] ?? 0).toFixed(2))), h("td", {}, s.mg_number.toFixed(1))));
  $("table").replaceChildren(h("thead", {}, head), h("tbody", {}, rows));
}

function csv() {
  const steps = ST.steps, names = phaseNames(steps), ox = Object.keys(steps[0].liquid);
  const lines = [["T_C", "P_bar", "melt_wt_pct", "melt_g", "magma_type", "TAS_name", ...names.flatMap((n) => [`${n}_g`, `${n}_composition`]), ...ox.map((o) => `melt_${o}`), "melt_Mg#"].join(",")];
  for (const s of steps) {
    lines.push([s.t, s.p, s.melt_pct.toFixed(4), s.liquid_grams.toFixed(4), s.magma_type, `"${s.tas}"`,
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
  const has = st.steps.length > 0 || (st.inputs.process === "example" && st.status !== "idle");
  $("empty").style.display = has ? "none" : "";
  $("content").style.display = has ? "flex" : "none";
  const pill = $("pill-run");
  pill.style.display = st.status === "idle" ? "none" : "";
  pill.textContent = { running: `running ${st.progress[0]}/${st.progress[1]}`, done: "done", stopped: "stopped", error: "error" }[st.status] || "";
  pill.className = `pill ${st.status === "error" ? "bad" : st.status === "running" ? "warn" : "ok"}`;
  if (!has) return;
  const i = st.inputs;
  const melting = i.process === "melt";
  const ex = i.process === "example";
  $("title").textContent = ex ? `${i.name || "Mantle"}: from mantle to magma` : melting ? `${i.name || "Rock"}: ${i.mode === "fractional" ? "fractional" : "batch"} melting`
    : `${i.name || "Melt"}: ${i.mode} crystallization`;
  $("meta").textContent = ex
    ? `Melted at ${i.example.t_source} °C and ${i.example.p_source} bar; the melt cooled at ${i.p_bar} bar from ${i.t_start} to ${i.t_end} °C in ${i.step} °C steps, ${i.example.h2o} wt% water, oxygen ${i.fo2}`
    : `${i.t_start} to ${i.t_end} °C in ${i.step} °C steps, ${i.p_bar} bar, oxygen ${i.fo2 === "closed" ? "closed system" : i.fo2}`;
  $("types").replaceChildren(...(melting ? [] : typeSteps(st.steps).map((s) => h("span", { title: `TAS: ${s.tas}` }, `${s.magma_type} · ${s.t.toFixed(0)} °C`))));
  $("story").replaceChildren(...story(st).map((t) => h("p", {}, t)), ...(st.status === "error" ? [h("div", { class: "errbox" }, st.message)] : st.message && st.status === "done" ? [h("p", {}, st.message)] : []));
  $("notes").replaceChildren(...st.notes.map((n) => h("li", {}, n)));
  if (!st.steps.length) { $("table").replaceChildren(); return; }
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
