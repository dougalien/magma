"use strict";

async function api(path, body) {
  const r = await fetch(path, {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body || {}),
  });
  return r.json();
}

function connect(onState, onStatus) {
  let es = null, retry = null;
  function open() {
    es = new EventSource("/events");
    es.onmessage = (e) => { onStatus && onStatus(true); onState(JSON.parse(e.data)); };
    es.onerror = () => {
      onStatus && onStatus(false);
      es.close();
      clearTimeout(retry);
      retry = setTimeout(open, 2000);
    };
  }
  open();
}

function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "class") e.className = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (v !== false && v !== null && v !== undefined) e.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat()) {
    if (kid === null || kid === undefined || kid === false) continue;
    e.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return e;
}

// Chemical formulas with real subscripts and superscripts: SiO2, Fe2O3, Ca0.85Mg0.90Si1.96O6, Fe2+0.50,
// O10(OH)2, and end-member shorthand such as Fo85, An60 Ab38 Or2 or Ilm69 Hem12.  Only words that look like formulas
// change; other words (QFM, Mg#, NNO+1, °C) are left alone, so the same function serves prose.
const FORMULA = /^(?:[A-Z][a-z]?(?:[1-9][+-])?(?:\d+(?:\.\d+)?)?|\(OH\)(?:\d+)?)+$/;
const PART = /([A-Z][a-z]?|\(OH\))([1-9][+-])?(\d+(?:\.\d+)?)?/g;
function chemParts(text) {
  const out = [];
  const plain = (t) => { if (t) out.push({ t, k: "" }); };
  for (const word of String(text).split(/(\s+)/)) {
    const m = /^([(\[“"']*)(.*?)([)\].,;:!?”"']*)$/.exec(word);
    const core = m[2];
    const em = /^([A-Z][a-z]{2,3})(\d+(?:\.\d+)?)$/.exec(core);     // end member in percent: Ilm69, Hem12, Grs5
    if (em) { plain(m[1]); plain(em[1]); out.push({ t: em[2], k: "sub" }); plain(m[3]); continue; }
    if (!core || !FORMULA.test(core) || !/\d|[+-]$/.test(core) || /^[A-Z][a-z]?[+-]?$/.test(core)) { plain(word); continue; }
    plain(m[1]);
    for (const [, el, charge, n] of core.matchAll(PART)) {
      plain(el);
      if (charge) out.push({ t: charge, k: "sup" });
      if (n) out.push({ t: n, k: "sub" });
    }
    plain(m[3]);
  }
  return out;
}
// ... as DOM nodes, for h(tag, attrs, ...chem(text))
function chem(text) {
  return chemParts(text).map((p) => (p.k ? h(p.k, {}, p.t) : document.createTextNode(p.t)));
}
const ESC = { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" };
function esc(t) { return String(t).replace(/[&<>"]/g, (c) => ESC[c]); }
// ... as HTML text (for innerHTML), and as SVG text (tspans; dy moves a subscript down and back up)
function chemHTML(text) {
  return chemParts(text).map((p) => (p.k ? `<${p.k}>${esc(p.t)}</${p.k}>` : esc(p.t))).join("");
}
function chemSVG(text) {
  let shift = 0, out = "";
  for (const p of chemParts(text)) {
    const want = p.k === "sub" ? 3 : p.k === "sup" ? -5 : 0;
    const dy = want - shift;
    shift = want;
    out += p.k || dy ? `<tspan${dy ? ` dy="${dy}"` : ""}${p.k ? ' font-size="75%"' : ""}>${esc(p.t)}</tspan>` : esc(p.t);
  }
  return out;
}

function setIfIdle(input, value) {
  if (document.activeElement !== input) input.value = value;
}

const SUP = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵",
  "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻" };
function sup(n) { return String(n).split("").map((c) => SUP[c] || c).join(""); }

function fmtVal(v) {
  if (v === null || v === undefined || !isFinite(v)) return "–";
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a < 1e-3 || a >= 1e5) {
    const [m, e] = v.toExponential(2).split("e");
    return `${m}×10${sup(parseInt(e, 10))}`;
  }
  return Number(v.toPrecision(3)).toString();
}
