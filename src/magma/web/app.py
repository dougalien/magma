"""Magma on Streamlit: an Input window and a Results window, as in the desktop app.

The run itself is the desktop app's ``State`` (``magma.app.server``), unchanged: it runs in a
background thread one temperature step at a time.  Each user has a private workspace id in the page
address (``?w=...``); the Input window and the Results window opened from it share that id, so they
show the same run, live, in two browser windows.  Nothing here computes thermodynamics.
"""
from __future__ import annotations

import base64
import csv
import io
import secrets
from pathlib import Path
import threading
import time

import plotly.graph_objects as go
import streamlit as st

from .. import __version__
from ..app.server import PHASE_LABEL, State
from ..classify import TAS_FIELDS, TAS_LINES
from ..composition import CARDS, OXIDES, PRESETS
from ..differentiation import DEFAULTS as EXAMPLE_DEFAULTS
from ..thermo.phases import PHASES
from .story import bowen, story, type_steps

FEEDBACK_URL = "https://github.com/dougalien/magma/issues"
# Where the two windows live; a host app (Protolith) that mounts these pages elsewhere changes them,
# and can add panels to the windows through the hooks (each a function of the State).
PATHS = {"input": "./", "results": "./results"}
RESULTS_EXTRAS: list = []
APP_NAME = "Magma"
WORKSPACE_TTL = 12 * 3600          # a workspace untouched for 12 hours is dropped
# Muted series colours that sit with the Dougalien sand, slate and rust
PALETTE = ["#3E5C6B", "#B4532A", "#C49A55", "#5E8C6A", "#8A6A9A", "#4F86A0", "#A8644A", "#6B6F45"]
STATIC = Path(__file__).resolve().parents[1] / "app" / "static"
ICON = str(STATIC / "magma_512.png")
STRIPE = ("linear-gradient(90deg, #D6B27C 0 22%, #3E5C6B 22% 34%, #D6B27C 34% 62%, #B4532A 62% 70%, "
          "#D6B27C 70% 88%, #3E5C6B 88% 100%)")
OXIDE_INPUTS = [o for o in OXIDES if o != "CO2"]
BUFFERS = ["QFM", "NNO", "IW", "MH", "closed"]
MODE_HINTS = {
    ("crystallize", "equilibrium"): "The crystals stay with the melt and keep reacting with it (a closed system).",
    ("crystallize", "fractional"): "The crystals are taken out at each step; the next step starts from the remaining melt.",
    ("melt", "equilibrium"): "Batch melting: the melt stays with the rock and keeps reacting with it. Starts at the solidus.",
    ("melt", "fractional"): "The melt above 1% of the rock is taken out at each step and pooled; the residue keeps melting.",
    ("example", "fractional"): "Stage 1 is batch melting of the mantle rock; stage 2 is fractional crystallization of its melt.",
}


# ---- workspaces: one State per user, shared by that user's two windows ------------------------------
@st.cache_resource
def _registry() -> dict:
    return {"lock": threading.Lock(), "ws": {}}


def workspace() -> tuple[str, State]:
    """This user's workspace id (from the address, or a new one) and its run State."""
    reg = _registry()
    wid = st.query_params.get("w")
    if not wid or len(wid) > 64:
        wid = secrets.token_urlsafe(16)
        st.query_params["w"] = wid
    now = time.time()
    with reg["lock"]:
        for k in [k for k, (_, t) in reg["ws"].items() if now - t > WORKSPACE_TTL and k != wid]:
            if reg["ws"][k][0].status != "running":
                del reg["ws"][k]
        state = reg["ws"][wid][0] if wid in reg["ws"] else State()
        reg["ws"][wid] = (state, now)
    return wid, state


def follow(state: State, key: str) -> None:
    """Re-run the page when the run changes (a step lands, it finishes, or the other window starts a run)."""
    seen_key = f"seen-{key}"
    st.session_state.setdefault(seen_key, state.version)

    @st.fragment(run_every=1.0 if state.status == "running" else 2.0)
    def _watch():
        if state.version != st.session_state[seen_key]:
            st.session_state[seen_key] = state.version
            st.rerun()
    _watch()
    st.session_state[seen_key] = state.version


def _svg(name: str) -> str:
    return "data:image/svg+xml;base64," + base64.b64encode((STATIC / name).read_bytes()).decode()


def header(subtitle: str) -> None:
    """The Dougalien Systems bar: logo, app name, maker, and the logo's strata as a stripe."""
    st.markdown(
        f"""<div style="display:flex;align-items:center;gap:12px;padding:4px 0 10px;position:relative">
<img src="{_svg('dougalien_mark.svg')}" width="34" height="34" alt="Dougalien Systems">
<img src="{_svg('magma_mark.svg')}" width="30" height="30" alt="">
<span style="font-size:1.55rem;font-weight:600;color:#2B3134">Magma</span>
<a href="https://dougalien.com" target="_blank" style="color:#5b6468;text-decoration:none">by Dougalien Systems</a>
<span style="color:#6f777a;font-size:0.9rem">· {subtitle}</span></div>
<div style="height:4px;border-radius:2px;background:{STRIPE};margin:0 0 14px"></div>""",
        unsafe_allow_html=True)


def footer() -> None:
    st.divider()
    st.caption(f"Magma {__version__} · alpha · [Report a problem]({FEEDBACK_URL}) · Every number comes from "
               "Magma's own solver; every datum is cited in the notes.")


# ---- the Input window -------------------------------------------------------------------------------
def _defaults() -> None:
    ss = st.session_state
    if ss.get("_init"):
        return
    ss["_init"] = True
    _apply_card(CARDS[0])


def _set_composition(preset: str) -> None:
    ss = st.session_state
    comp = PRESETS[preset]
    for ox in OXIDE_INPUTS:
        ss[f"ox-{ox}"] = float(comp.wt.get(ox, 0.0))
    ss["preset"] = preset


def _apply_card(c: dict) -> None:
    ss = st.session_state
    _set_composition(c["preset"])
    ss["process"] = {"rock": "melt", "example": "example"}.get(c["kind"], "crystallize")
    if c["kind"] == "example":
        return
    ss["t_start"], ss["t_end"], ss["step"], ss["p_bar"] = (float(c["t_start"]), float(c["t_end"]),
                                                           float(c["step"]), float(c["p_bar"]))
    fo2 = c.get("fo2", "QFM")
    if fo2 == "closed":
        ss["buffer"], ss["offset"] = "closed", 0.0
    else:
        base = fo2.rstrip("0123456789.+-") or "QFM"
        ss["buffer"], ss["offset"] = base, float(fo2[len(base):] or 0)
    ss["mode"] = c["mode"]


def _example_card() -> None:
    ss = st.session_state
    _set_composition("KLB-1 peridotite")
    ss["process"] = "example"
    for k in ("t_source", "p_source", "p_crust", "h2o", "t_end", "step"):
        ss[f"ex-{k}"] = float(EXAMPLE_DEFAULTS[k])
    ss["ex-fo2"] = "NNO"


def _cards() -> None:
    st.markdown("**Start from a card.** Each card fills in a composition and conditions that show it well. "
                "Change anything before you press Run.")
    for title, kind in (("Melts to cool", "melt"), ("Rocks to melt", "rock")):
        st.caption(title)
        group = [c for c in CARDS if c["kind"] == kind]
        cols = st.columns(len(group))
        for col, c in zip(cols, group):
            col.button(c["title"], key=f"card-{c['title']}", help=c["blurb"], on_click=_apply_card, args=(c,),
                       width="stretch")
    st.caption("The whole story")
    st.button("Mantle to magma: Bowen's series", help="Melt the mantle, then cool the melt: Bowen's series and "
              "the magma types.", on_click=_example_card)


def _inputs() -> dict:
    """The run request, in the same form the desktop control window sends."""
    ss = st.session_state
    proc = ss["process"]
    bulk = {ox: ss[f"ox-{ox}"] for ox in OXIDE_INPUTS if ss.get(f"ox-{ox}")}
    phases = ss["phases"]
    name = PRESETS[ss["preset"]].name if ss.get("preset") in PRESETS else ""
    if proc == "example":
        g = lambda k: ss[f"ex-{k}"]
        return {"name": name, "bulk": bulk, "process": "example", "mode": "fractional", "phases": phases,
                "t_start": g("t_source"), "t_end": g("t_end"), "step": g("step"), "p_bar": g("p_crust"),
                "fo2": ss["ex-fo2"], "engine": ss["engine"],
                "example": {"t_source": g("t_source"), "p_source": g("p_source"), "h2o": g("h2o"), "fo2_source": "QFM"}}
    off = ss["offset"]
    fo2 = "closed" if ss["buffer"] == "closed" else (f"{ss['buffer']}{'+' if off > 0 else ''}{off:g}" if off else ss["buffer"])
    return {"name": name, "bulk": bulk, "t_start": ss["t_start"], "t_end": ss["t_end"], "step": ss["step"],
            "p_bar": ss["p_bar"], "fo2": fo2, "mode": ss["mode"], "process": proc, "engine": ss["engine"],
            "phases": phases}


def input_page() -> None:
    st.set_page_config(page_title="Magma | Input", page_icon=ICON, layout="wide")
    wid, state = workspace()
    _defaults()
    ss = st.session_state
    ss.setdefault("engine", "hp")
    ss.setdefault("mode", "equilibrium")
    ss.setdefault("phases", list(PHASES))
    for k in ("t_source", "p_source", "p_crust", "h2o", "t_end", "step"):
        ss.setdefault(f"ex-{k}", float(EXAMPLE_DEFAULTS[k]))
    ss.setdefault("ex-fo2", "NNO")

    top = st.columns([4, 1.4])
    with top[0]:
        header("cooling a melt, melting a rock")
    with top[1]:
        st.link_button("Open results window", f"{PATHS['results']}?w={wid}", width="stretch",
                       help="Opens the Results window in a new browser window; put it on a second screen. "
                            "It follows this window's runs.")

    with st.container(border=True):
        _cards()

    proc = ss["process"]
    left, right = st.columns([1, 1], gap="large")
    with left:
        with st.container(border=True):
            st.markdown("**" + {"melt": "Starting rock", "example": "Mantle rock"}.get(proc, "Starting melt") + "**")
            st.selectbox("Composition", list(PRESETS), key="preset", on_change=lambda: _set_composition(ss["preset"]))
            st.caption(PRESETS[ss["preset"]].source)
            cols = st.columns(4)
            for i, ox in enumerate(OXIDE_INPUTS):
                cols[i % 4].number_input(f"{ox}, wt%", key=f"ox-{ox}", min_value=0.0, step=0.1, format="%.3f")
            total = sum(ss[f"ox-{ox}"] or 0 for ox in OXIDE_INPUTS)
            st.caption(f"Total **{total:.2f}** wt%. Values are normalized to 100. Give iron as FeO and Fe2O3, "
                       "or all as FeO.")
    with right:
        with st.container(border=True):
            st.markdown("**Process**")
            st.radio("Engine", ["hp", "melts"], key="engine", horizontal=True,
                     format_func=lambda e: {"hp": "Holland & Powell", "melts": "MELTS"}[e],
                     help="Holland and Powell (as in MAGEMin) finds all phases at once and follows paths below "
                          "the solidus; MELTS uses the Ghiorso and Sack (1995) models.")
            st.radio("Process", ["crystallize", "melt", "example"], key="process", horizontal=True,
                     format_func=lambda p: {"crystallize": "Cool a melt", "melt": "Heat a rock",
                                            "example": "Mantle to magma"}[p])
            if proc != "example":
                st.radio("Mode", ["equilibrium", "fractional"], key="mode", horizontal=True,
                         format_func=str.capitalize)
            st.caption(MODE_HINTS[(proc, "fractional" if proc == "example" else ss["mode"])])
            st.multiselect("Minerals allowed", list(PHASES), key="phases",
                           format_func=lambda n: PHASE_LABEL.get(n, n))
            st.caption("Water (fluid) appears when the melt is saturated; at 1 bar any water leaves the melt at once.")
        with st.container(border=True):
            if proc == "example":
                st.markdown("**Mantle to magma**")
                st.caption("Two stages. First the mantle rock partly melts at the source (batch melting): that melt "
                           "is the primary basalt. Then the melt rises to the crust and cools, and the crystals are "
                           "taken out as they form (fractional crystallization).")
                st.markdown("1. Melting the mantle")
                a, b = st.columns(2)
                a.number_input("Melting temperature, °C", key="ex-t_source", step=10.0)
                b.number_input("Pressure, bar", key="ex-p_source", min_value=1.0, step=500.0)
                st.markdown("2. Cooling the melt in the crust")
                a, b = st.columns(2)
                a.number_input("Pressure, bar", key="ex-p_crust", min_value=1.0, step=500.0)
                b.number_input("Water added, wt%", key="ex-h2o", min_value=0.0, step=0.5)
                a.number_input("Cool to, °C", key="ex-t_end", step=10.0)
                b.number_input("Step, °C", key="ex-step", min_value=1.0, step=1.0)
                st.selectbox("Oxygen", ["NNO", "QFM", "NNO+1"], key="ex-fo2")
            else:
                st.markdown("**Conditions**")
                a, b, c = st.columns(3)
                a.number_input("Start, °C", key="t_start", step=10.0)
                b.number_input("End, °C", key="t_end", step=10.0)
                c.number_input("Step, °C", key="step", min_value=0.5, step=1.0)
                a, b, c = st.columns(3)
                a.number_input("Pressure, bar", key="p_bar", min_value=1.0, step=100.0)
                b.selectbox("Oxygen", BUFFERS, key="buffer",
                            format_func=lambda x: "Closed to oxygen" if x == "closed" else f"{x} buffer")
                if ss["buffer"] != "closed":
                    c.number_input("Offset, log units", key="offset", step=0.5)

    run_bar(state, input_page=True)
    footer()


def run_bar(state: State, input_page: bool) -> None:
    running = state.status == "running"
    cols = st.columns([1, 1, 6])
    if input_page:
        if cols[0].button("Run", type="primary", disabled=running, width="stretch"):
            try:
                state.start(_inputs())
                st.session_state[f"seen-input"] = state.version
                st.rerun()
            except (ValueError, TypeError) as exc:
                st.error(str(exc))
    if running and cols[1].button("Stop", width="stretch"):
        state.stop()
    if running:
        k, n = state.progress
        cols[2].progress(min(k / n, 1.0) if n else 0.0, text=state.message or f"Running: step {k} of {n}")
    elif state.status in ("done", "stopped"):
        cols[2].success(f"Run {state.run_id}: {state.status}. " + (state.message or ""), icon="✅")
    elif state.status == "error":
        cols[2].error(state.message)
    follow(state, "input" if input_page else "results")


# ---- the Results window -----------------------------------------------------------------------------
def _phase_names(steps):
    seen = []
    for s in steps:
        for p in s["solids"]:
            if p["name"] not in seen:
                seen.append(p["name"])
    return seen


def _label(steps, name):
    for s in steps:
        for p in s["solids"]:
            if p["name"] == name:
                return p["label"]
    return PHASE_LABEL.get(name, name)


def _fig(xlabel, ylabel, **kw) -> go.Figure:
    fig = go.Figure()
    layout = dict(height=480, margin=dict(l=10, r=10, t=30, b=10), xaxis_title=xlabel, yaxis_title=ylabel,
                  hovermode="x unified", legend=dict(orientation="h", y=-0.18), colorway=PALETTE,
                  paper_bgcolor="#F3F1EB", plot_bgcolor="#F3F1EB", font=dict(color="#2B3134"))
    fig.update_layout(**{**layout, **kw})          # a caller's own settings (hovermode for TAS) win
    fig.update_xaxes(gridcolor="#DDD8CC", zerolinecolor="#CFC8B8")
    fig.update_yaxes(gridcolor="#DDD8CC", zerolinecolor="#CFC8B8")
    return fig


def _mass_plot(steps) -> go.Figure:
    x = [s["t"] for s in steps]
    fig = _fig("temperature, °C", "grams per 100 g of starting material")
    fig.add_scatter(x=x, y=[s["liquid_grams"] for s in steps], name="Melt", line=dict(width=3))
    if any(s.get("pooled_grams", 0) > 0 for s in steps):
        fig.add_scatter(x=x, y=[s.get("pooled_grams", 0) for s in steps], name="Melt removed (total)",
                        line=dict(dash="dash"))
    for n in _phase_names(steps):
        fig.add_scatter(x=x, y=[next((p["grams"] for p in s["solids"] if p["name"] == n), 0) for s in steps],
                        name=_label(steps, n))
    return fig


def _liquid_plot(steps) -> go.Figure:
    x = [s["t"] for s in steps]
    fig = _fig("temperature, °C", "melt, wt%")
    for i, ox in enumerate(o for o in steps[0]["liquid"] if o != "H2O"):
        fig.add_scatter(x=x, y=[s["liquid"].get(ox) for s in steps], name=ox,
                        line=dict(dash="dash" if i >= 8 else None))
    return fig


MINERAL_NUMBERS = [("olivine", "Fo", "Olivine Fo"), ("feldspar", "An", "Plagioclase An"),
                   ("alkali-feldspar", "Or", "K-feldspar Or"), ("hornblende", "Mg#", "Hornblende Mg#"),
                   ("biotite", "Mg#", "Biotite Mg#"), ("clinopyroxene", "Mg#", "Augite Mg#"),
                   ("clinopyroxene", "Wo", "Augite Wo"), ("pigeonite", "Mg#", "Pigeonite Mg#"),
                   ("orthopyroxene", "Mg#", "Orthopyroxene Mg#"), ("spinel", "Cr#", "Spinel Cr#"),
                   ("spinel", "Usp", "Spinel ulvöspinel"), ("rhombohedral", "Ilm", "Ilmenite (FeTiO3)")]


def _mineral_plot(steps) -> go.Figure:
    x = [s["t"] for s in steps]
    fig = _fig("temperature, °C", "mol%")
    for name, key, lab in MINERAL_NUMBERS:
        y = []
        for s in steps:
            p = next((q for q in s["solids"] if q["name"] == name), None)
            y.append(p["numbers"][key] if p and key in (p.get("numbers") or {}) and p["grams"] > 1e-9 else None)
        if any(v is not None for v in y):
            fig.add_scatter(x=x, y=y, name=lab)
    fig.add_scatter(x=x, y=[s["mg_number"] for s in steps], name="Melt Mg#", line=dict(dash="dash"))
    return fig


def _tas_plot(st_) -> go.Figure:
    steps = st_["steps"]
    fig = _fig("SiO2, wt% (volatile-free)", "Na2O + K2O, wt%", hovermode="closest")
    for ln in TAS_LINES:
        fig.add_scatter(x=[p[0] for p in ln], y=[p[1] for p in ln], mode="lines", showlegend=False,
                        hoverinfo="skip", line=dict(color="#aaaaaa", width=1))
    for name, pts in TAS_FIELDS.items():
        fig.add_annotation(x=sum(p[0] for p in pts) / len(pts), y=sum(p[1] for p in pts) / len(pts),
                           text=name.replace(" or ", "/"), showarrow=False, font=dict(size=10, color="#888888"))
    if steps:
        fig.add_scatter(x=[s["tas_xy"][0] for s in steps], y=[s["tas_xy"][1] for s in steps], name="Melt",
                        mode="lines+markers", marker=dict(size=4), line=dict(color=PALETTE[1], width=3),
                        text=[f"{s['t']:.0f} °C, {s['liquid_grams']:.1f} g melt: {s['tas']}" for s in steps],
                        hovertemplate="%{text}<extra></extra>")
        ts = type_steps(steps)
        fig.add_scatter(x=[s["tas_xy"][0] for s in ts], y=[s["tas_xy"][1] for s in ts], mode="markers+text",
                        text=[f"{s['t']:.0f} °C" for s in ts], textposition="top right", showlegend=False,
                        marker=dict(size=9, color=PALETTE[1]), hoverinfo="skip")
    src = st_.get("source")
    if src and src.get("primary"):
        fig.add_scatter(x=[src["tas_xy"][0]], y=[src["tas_xy"][1]], mode="markers+text", text=["primary melt"],
                        textposition="bottom right", name="Primary melt", marker=dict(size=10, color=PALETTE[0]))
    fig.update_xaxes(range=[40, 80])
    fig.update_yaxes(range=[0, 16])
    return fig


def _bowen_view(steps) -> None:
    b = bowen(steps)

    def boxes(items):
        for i, (title, text, formed) in enumerate(items):
            if i:
                st.markdown("<div style='text-align:center'>↓</div>", unsafe_allow_html=True)
            with st.container(border=True):
                st.markdown(f"**{title}**" if formed else f"<span style='color:gray'><b>{title}</b></span>",
                            unsafe_allow_html=True)
                st.caption(text)
    t0, t1 = steps[0]["t"], steps[-1]["t"]
    st.caption(f"↑ hot: {t0:.0f} °C · cool: {t1:.0f} °C ↓")
    a, c = st.columns(2)
    with a:
        st.markdown("**Discontinuous series**")
        boxes(b["discontinuous"])
    with c:
        st.markdown("**Continuous series**")
        boxes(b["continuous"])
    st.markdown("**Late minerals**")
    for col, item in zip(st.columns(len(b["late"])), b["late"]):
        with col:
            boxes([item])
    st.caption(f"Also formed, outside Bowen's scheme: {b['also']}." if b["also"]
               else "Grey boxes: minerals that did not form, or that Magma cannot make yet.")


def _table_rows(steps) -> tuple[list[str], list[list]]:
    names, ox = _phase_names(steps), list(steps[0]["liquid"])
    head = (["T_C", "P_bar", "melt_wt_pct", "melt_g", "magma_type", "TAS_name"]
            + [f for n in names for f in (f"{n}_g", f"{n}_composition")] + [f"melt_{o}" for o in ox] + ["melt_Mg#"])
    rows = []
    for s in steps:
        r = [s["t"], s["p"], round(s["melt_pct"], 4), round(s["liquid_grams"], 4), s["magma_type"], s["tas"]]
        for n in names:
            p = next((q for q in s["solids"] if q["name"] == n), None)
            r += [round(p["grams"], 5), p["formula"]] if p else [0.0, ""]
        r += [round(s["liquid"].get(o, 0.0), 4) for o in ox] + [round(s["mg_number"], 2)]
        rows.append(r)
    return head, rows


def results_page() -> None:
    st.set_page_config(page_title="Magma | Results", page_icon=ICON, layout="wide")
    wid, state = workspace()
    top = st.columns([4, 1.4])
    with top[0]:
        header("results")
    with top[1]:
        st.link_button("Back to the input window", f"{PATHS['input']}?w={wid}", width="stretch")
    snap = state.snapshot()
    has = snap["steps"] or (snap["inputs"].get("process") == "example" and snap["status"] != "idle")
    if not has:
        st.info("No run yet. Set up a run in the Input window and press Run; this window follows it live.")
        run_bar(state, input_page=False)
        footer()
        return
    i = snap["inputs"]
    proc = i.get("process")
    if proc == "example":
        ex = i.get("example", {})
        title = f"{i.get('name') or 'Mantle'}: from mantle to magma"
        meta = (f"Melted at {ex.get('t_source')} °C and {ex.get('p_source')} bar; the melt cooled at {i['p_bar']} bar "
                f"from {i['t_start']} to {i['t_end']} °C in {i['step']} °C steps, {ex.get('h2o')} wt% water, oxygen {i['fo2']}")
    else:
        title = (f"{i.get('name') or 'Rock'}: {'fractional' if i['mode'] == 'fractional' else 'batch'} melting"
                 if proc == "melt" else f"{i.get('name') or 'Melt'}: {i['mode']} crystallization")
        meta = (f"{i['t_start']} to {i['t_end']} °C in {i['step']} °C steps, {i['p_bar']} bar, oxygen "
                + ("closed system" if i["fo2"] == "closed" else str(i["fo2"])))
    st.subheader(title)
    st.caption(meta + f" · engine: {'Holland & Powell' if i.get('engine') == 'hp' else 'MELTS'} · run {snap['run_id']}")
    run_bar(state, input_page=False)
    steps = snap["steps"]
    if steps and proc != "melt":
        st.markdown(" → ".join(f"`{s['magma_type']} · {s['t']:.0f} °C`" for s in type_steps(steps)))
    with st.container(border=True):
        st.markdown("**What we found**")
        for t in story(snap):
            st.markdown(t)
        if snap["status"] == "done" and snap["message"]:
            st.markdown(snap["message"])
    if not steps:
        footer()
        return
    views = {"Masses": "mass", "Melt composition": "liquid", "Mineral compositions": "minerals",
             "TAS diagram": "tas", "Bowen's series": "bowen"}
    view = views[st.segmented_control("View", list(views), default="Masses", key="view") or "Masses"]
    if view == "bowen":
        _bowen_view(steps)
    else:
        fig = {"mass": _mass_plot, "liquid": _liquid_plot, "minerals": _mineral_plot}.get(view)
        st.plotly_chart(fig(steps) if fig else _tas_plot(snap), width="stretch")
        if view == "tas":
            st.caption("Total alkali-silica diagram (Le Bas et al. 1986); dots mark a new magma type.")
        elif view == "mass" and i.get("mode") == "fractional":
            st.caption("Mineral masses are what is left in the residue; the melt left in it is 1%." if proc == "melt"
                       else "Crystal masses are the totals removed so far.")
    head, rows = _table_rows(steps)
    st.markdown("**Every step**")
    st.dataframe([dict(zip(head, r)) for r in rows], width="stretch", hide_index=True)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(head)
    w.writerows(rows)
    st.download_button("Download the table (CSV)", buf.getvalue(), file_name="magma_run.csv", mime="text/csv")
    for extra in RESULTS_EXTRAS:
        extra(state, snap)
    with st.expander("Notes on the models"):
        for n in snap["notes"]:
            st.markdown(f"- {n}")
    footer()
