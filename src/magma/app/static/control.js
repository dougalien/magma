"use strict";

const $ = (id) => document.getElementById(id);
let INFO = null, MODE = "equilibrium";
const HINTS = {
  equilibrium: "The crystals stay with the melt and keep reacting with it (a closed system).",
  fractional: "The crystals are taken out at each step; the next step starts from the remaining melt.",
};

function openResults() { window.open("/results", "magma-results"); }

function fillPreset(name) {
  const p = INFO.presets[name];
  for (const ox of INFO.oxides) $(`ox-${ox}`).value = p && p.wt[ox] ? p.wt[ox] : "";
  $("preset-source").textContent = p ? p.source : "Type your own composition, in oxide wt%.";
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
  const name = $("preset").value === "custom" ? "Custom melt" : $("preset").value;
  return { name, bulk, t_start: +$("t_start").value, t_end: +$("t_end").value, step: +$("step").value,
    p_bar: +$("p_bar").value, fo2, mode: MODE, phases };
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
    idle: "", running: `Running: step ${k} of ${n}${last ? `, ${last.t.toFixed(0)} C` : ""}`,
    done: `Done: ${st.steps.length} steps. The results window has the plots.`,
    stopped: `Stopped after ${st.steps.length} steps.`, error: "",
  }[st.status];
  $("status").textContent = msg;
  if (st.status === "error") { $("err").textContent = st.message; $("err").style.display = ""; }
}

async function init() {
  INFO = await (await fetch("/api/presets")).json();
  const sel = $("preset");
  for (const name of Object.keys(INFO.presets)) sel.append(h("option", { value: name }, INFO.presets[name].name || name));
  sel.append(h("option", { value: "custom" }, "Custom"));
  sel.addEventListener("change", () => fillPreset(sel.value));
  for (const ox of INFO.oxides) {
    const inp = h("input", { id: `ox-${ox}`, type: "number", min: "0", step: "any", oninput: () => { total(); sel.value = "custom"; $("preset-source").textContent = "Your own composition."; } });
    $("oxides").append(h("label", { class: "f" }, ox, inp));
  }
  for (const p of INFO.phases) {
    const c = h("span", { class: "chip on", "data-name": p.name, role: "checkbox", "aria-checked": "true", tabindex: "0" }, p.label);
    c.addEventListener("click", () => { c.classList.toggle("on"); c.setAttribute("aria-checked", c.classList.contains("on")); });
    $("phases").append(c);
  }
  document.querySelectorAll("#mode button").forEach((b) => b.addEventListener("click", () => {
    MODE = b.dataset.mode;
    document.querySelectorAll("#mode button").forEach((x) => x.classList.toggle("on", x === b));
    $("mode-hint").textContent = HINTS[MODE];
  }));
  $("mode-hint").textContent = HINTS[MODE];
  $("buffer").addEventListener("change", () => { $("offset-wrap").style.display = $("buffer").value === "closed" ? "none" : ""; });
  $("run").addEventListener("click", run);
  $("stop").addEventListener("click", () => api("/api/stop"));
  $("open-results").addEventListener("click", openResults);
  fillPreset(sel.value);
  connect(show, (ok) => { $("pill-link").textContent = ok ? "running" : "not connected"; $("pill-link").className = `pill ${ok ? "ok" : "bad"}`; });
}
init();
