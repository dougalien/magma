"use strict";

const $ = (id) => document.getElementById(id);
let INFO = null, MODE = "equilibrium", PROCESS = "crystallize", ENGINE = "hp";
const HINTS = {
  crystallize: {
    equilibrium: "The crystals stay with the melt and keep reacting with it (a closed system).",
    fractional: "The crystals are taken out at each step; the next step starts from the remaining melt.",
  },
  example: {
    fractional: "Stage 1 is batch melting of the mantle rock; stage 2 is fractional crystallization of its melt.",
  },
  melt: {
    equilibrium: "Batch melting: the melt stays with the rock and keeps reacting with it. Starts at the solidus.",
    fractional: "The melt above 1% of the rock is taken out at each step and pooled; the residue keeps melting.",
  },
};
function hint() { $("mode-hint").textContent = HINTS[PROCESS][PROCESS === "example" ? "fractional" : MODE]; }
function setProcess(p) {
  PROCESS = p;
  document.querySelectorAll("#process button").forEach((x) => x.classList.toggle("on", x.dataset.process === p));
  $("mode-eq").textContent = p === "melt" ? "Batch" : "Equilibrium";
  $("start-legend").textContent = p === "melt" ? "Starting rock" : p === "example" ? "Mantle rock" : "Starting melt";
  const ex = p === "example";
  $("example-box").style.display = ex ? "" : "none";
  $("conditions-box").style.display = ex ? "none" : "";
  $("mode").style.display = ex ? "none" : "";
  if (ex) {
    if (INFO.presets["KLB-1 peridotite"]) { $("preset").value = "KLB-1 peridotite"; fillPreset("KLB-1 peridotite"); }
    hint();
    return;
  }
  const a = +$("t_start").value, b = +$("t_end").value;
  if ((p === "melt") === (a > b)) { $("t_start").value = b; $("t_end").value = a; }
  if (p === "melt" && INFO.presets["KLB-1 peridotite"] && $("preset").value !== "custom") {
    $("preset").value = "KLB-1 peridotite"; fillPreset("KLB-1 peridotite");
    $("t_start").value = 1250; $("t_end").value = 1500; $("step").value = 10; $("p_bar").value = 10000;
  }
  hint();
}

function openResults() { window.open("/results", "magma-results"); }

function fillPreset(name) {
  const p = INFO.presets[name];
  for (const ox of INFO.oxides) $(`ox-${ox}`).value = p && p.wt[ox] ? p.wt[ox] : "";
  $("preset-source").replaceChildren(...chem(p ? p.source : "Type your own composition, in oxide wt%."));
  total();
}

function total() {
  let s = 0;
  for (const ox of INFO.oxides) s += parseFloat($(`ox-${ox}`).value) || 0;
  $("total").textContent = s.toFixed(2);
}

function inputs() {
  const bulk = {};
  for (const ox of INFO.oxides) {
    const v = $(`ox-${ox}`).value.trim();
    if (v) bulk[ox] = parseFloat(v);
  }
  const b = $("buffer").value, off = parseFloat($("offset").value) || 0;
  const fo2 = b === "closed" ? "closed" : off ? `${b}${off > 0 ? "+" : ""}${off}` : b;
  const phases = [...document.querySelectorAll("#phases .chip.on")].map((c) => c.dataset.name);
  const name = $("preset").value === "custom" ? (PROCESS === "crystallize" ? "Custom melt" : "Custom rock") : $("preset").value;
  if (PROCESS === "example") {
    const v = (id) => +$(`ex-${id}`).value;
    return { name, bulk, process: "example", mode: "fractional", phases, t_start: v("t_source"), t_end: v("t_end"),
      step: v("step"), p_bar: v("p_crust"), fo2: $("ex-fo2").value,
      engine: ENGINE, example: { t_source: v("t_source"), p_source: v("p_source"), h2o: v("h2o"), fo2_source: "QFM" } };
  }
  return { name, bulk, t_start: +$("t_start").value, t_end: +$("t_end").value, step: +$("step").value,
    p_bar: +$("p_bar").value, fo2, mode: MODE, process: PROCESS, engine: ENGINE, phases };
}

async function run() {
  $("err").style.display = "none";
  const r = await api("/api/run", inputs());
  if (r.error) { $("err").textContent = r.error; $("err").style.display = ""; return; }
  openResults();
}

function show(st) {
  const running = st.status === "running";
  $("run").disabled = running;
  $("stop").style.display = running ? "" : "none";
  const [k, n] = st.progress;
  $("prog").style.display = running ? "" : "none";
  $("prog").max = n || 1; $("prog").value = k;
  const last = st.steps.length ? st.steps[st.steps.length - 1] : null;
  const msg = {
    idle: "", running: st.message || `Running: step ${k} of ${n}${last ? `, ${last.t.toFixed(0)} C` : ""}`,
    done: `Done: ${st.steps.length} steps. The results window has the plots.`,
    stopped: `Stopped after ${st.steps.length} steps.`, error: "",
  }[st.status];
  $("status").textContent = msg;
  if (st.status === "error") { $("err").textContent = st.message; $("err").style.display = ""; }
}

function applyCard(c) {
  setProcess(c.kind === "rock" ? "melt" : c.kind === "example" ? "example" : "crystallize");
  if (c.kind === "example") return;
  $("preset").value = c.preset; fillPreset(c.preset);
  $("t_start").value = c.t_start; $("t_end").value = c.t_end; $("step").value = c.step; $("p_bar").value = c.p_bar;
  const m = /^([A-Z]+)([+-][\d.]+)?$/.exec(c.fo2 || "QFM");
  $("buffer").value = c.fo2 === "closed" ? "closed" : (m ? m[1] : "QFM");
  $("offset").value = m && m[2] ? parseFloat(m[2]) : 0;
  $("offset-wrap").style.display = c.fo2 === "closed" ? "none" : "";
  MODE = c.mode;
  document.querySelectorAll("#mode button").forEach((x) => x.classList.toggle("on", x.dataset.mode === MODE));
  hint();
  document.querySelectorAll(".cards .scard").forEach((x) => x.classList.toggle("on", x.dataset.title === c.title));
}

function cardEl(c) {
  const p = INFO.presets[c.preset];
  const si = p ? (100 * p.wt.SiO2 / Object.values(p.wt).reduce((a, b) => a + b, 0)).toFixed(0) : null;
  const meta = c.kind === "example" ? "KLB-1, then cool its melt" :
    `${si}% SiO2 · ${c.kind === "rock" ? `${c.p_bar / 1000} kbar, ${c.t_start} to ${c.t_end} °C` : `${c.t_start} to ${c.t_end} °C, ${c.p_bar / 1000} kbar`}`;
  return h("button", { class: "scard", "data-title": c.title, onclick: () => applyCard(c) },
    h("b", {}, c.title), h("span", { class: "blurb" }, c.blurb), h("span", { class: "meta" }, ...chem(meta)));
}

async function init() {
  INFO = await (await fetch("/api/presets")).json();
  const sel = $("preset");
  for (const name of Object.keys(INFO.presets)) sel.append(h("option", { value: name }, INFO.presets[name].name || name));
  sel.append(h("option", { value: "custom" }, "Custom"));
  sel.addEventListener("change", () => fillPreset(sel.value));
  for (const ox of INFO.oxides) {
    const inp = h("input", { id: `ox-${ox}`, type: "number", min: "0", step: "any", oninput: () => { total(); sel.value = "custom"; $("preset-source").textContent = "Your own composition."; } });
    $("oxides").append(h("label", { class: "f" }, ...chem(ox), inp));
  }
  for (const p of INFO.phases) {
    const c = h("span", { class: "chip on", "data-name": p.name, role: "checkbox", "aria-checked": "true", tabindex: "0" }, p.label);
    c.addEventListener("click", () => { c.classList.toggle("on"); c.setAttribute("aria-checked", c.classList.contains("on")); });
    $("phases").append(c);
  }
  document.querySelectorAll("#mode button").forEach((b) => b.addEventListener("click", () => {
    MODE = b.dataset.mode;
    document.querySelectorAll("#mode button").forEach((x) => x.classList.toggle("on", x === b));
    hint();
  }));
  document.querySelectorAll("#process button").forEach((b) => b.addEventListener("click", () => setProcess(b.dataset.process)));
  document.querySelectorAll("#engine button").forEach((b) => b.addEventListener("click", () => {
    ENGINE = b.dataset.engine;
    document.querySelectorAll("#engine button").forEach((x) => x.classList.toggle("on", x === b));
  }));
  hint();
  $("buffer").addEventListener("change", () => { $("offset-wrap").style.display = $("buffer").value === "closed" ? "none" : ""; });
  $("run").addEventListener("click", run);
  $("stop").addEventListener("click", () => api("/api/stop"));
  $("open-results").addEventListener("click", openResults);
  fillPreset(sel.value);
  for (const c of INFO.cards || []) $(c.kind === "rock" ? "cards-rock" : "cards-melt").append(cardEl(c));
  $("cards-example").append(cardEl({ kind: "example", title: "Mantle to magma: Bowen's series",
    blurb: "Partly melt the mantle, then cool that melt with fractional crystallization through Bowen's series." }));
  connect(show, (ok) => { $("pill-link").textContent = ok ? "running" : "not connected"; $("pill-link").className = `pill ${ok ? "ok" : "bad"}`; });
}
init();
