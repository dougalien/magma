"""A small local web server for the two windows (control, results).

Standard library only, like Reaction's.  Both windows follow one shared state through a
server-sent-event stream.  The server listens on the loopback address only and rejects requests
whose Host or Origin is not the loopback address, so a web page you happen to have open cannot
drive it.  A run goes in a background thread, one temperature step at a time, so the windows
show progress and Stop takes effect between steps.
"""
from __future__ import annotations

import json
import mimetypes
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

from .. import __version__
import numpy as np

from ..classify import BOWEN_CONTINUOUS, BOWEN_DISCONTINUOUS, BOWEN_LATE, TAS_FIELDS, TAS_LINES, magma_type, tas_name, tas_point
from ..composition import CARDS, OXIDES, PRESETS, Composition
from ..differentiation import DEFAULTS as EXAMPLE_DEFAULTS, iter_differentiate, melt_source, primary_with_water
from ..melting import iter_melt
from ..path import iter_crystallize
from ..hp import engine as hp_engine
from ..thermo.phases import PHASES

STATIC = Path(__file__).resolve().parent / "static"
PAGES = {"/": "control.html", "/control": "control.html", "/results": "results.html"}
DEFAULT_PORT = 8766            # Reaction uses its own port; both apps can run at once
MAX_BODY = 200_000
PHASE_LABEL = {"olivine": "Olivine", "feldspar": "Plagioclase", "clinopyroxene": "Augite", "pigeonite": "Pigeonite", "orthopyroxene": "Orthopyroxene",
               "spinel": "Spinel", "titanomagnetite": "Titanomagnetite", "rhombohedral": "Ilmenite", "garnet": "Garnet",
               "hornblende": "Hornblende", "biotite": "Biotite", "alkali-feldspar": "K-feldspar",
               "muscovite": "Muscovite", "quartz": "Quartz", "tridymite": "Tridymite",
               "water": "Water (fluid)"}
PHASE_LABEL.update(hp_engine.LABEL)
ENGINES = ("melts", "hp")

HP_NOTES = [
    "Engine: Holland and Powell. End-member data from the Holland and Powell (2011) data set (ds636); "
    "solution models of the igneous set of Holland, Green and Powell (2018) with the amphibole of Green "
    "et al. (2016) and the micas, cordierite, epidote and chlorite of White et al. (2014); liquid and "
    "aqueous fluid of Holland, Green and Powell (2018). Parameter values as published and as used in "
    "MAGEMin (Riel et al. 2022); the equations are Magma's own code, checked against MAGEMin.",
    "Each step is a Gibbs energy minimization over all the phases at once, so the path can go below the "
    "solidus and phases MELTS lacks can appear: cordierite, epidote, chlorite, rutile, sphene, the "
    "Al2SiO5 minerals, corundum and nepheline. A feldspar with Or above 30% is shown as K-feldspar, a "
    "clinopyroxene with Wo below 25% as pigeonite, a spinel with more than half magnetite and "
    "ulvöspinel as titanomagnetite.",
    "Not in this engine: MnO, P2O5, NiO, CoO and CO2 (left out), apatite and zircon. The liquid model is "
    "calibrated for 1 to 30 kbar and is least certain at very low pressure.",
    "Iron: the system is open to oxygen at the chosen buffer (or closed, keeping its starting Fe2O3).",
    "Masses are grams per 100 g of starting material. In fractional crystallization the solids of each "
    "step are removed, and the solid masses shown are the totals removed so far. Melting shows the steps "
    "from the first melt (the solidus) upward.",
]

NOTES = [
    "Liquid: Ghiorso and Sack (1995) regular-solution model, with the parameter values used in MELTS.",
    "Solids: olivine, plagioclase (Elkins and Grove 1990), augite, pigeonite and orthopyroxene (the "
    "MELTS pyroxene model of Sack and Ghiorso 1994; augite and pigeonite are two copies of one model, as "
    "in MELTS, and pigeonite is left out above 15 kbar), spinel and titanomagnetite (two copies of the "
    "spinel model of Sack and Ghiorso 1991, the second without chromite, as MELTS allows two spinels), "
    "ilmenite (Ghiorso and Evans 2008), garnet (Berman 1990), and the late minerals as in MELTS v1.0.2 "
    "(rhyolite-MELTS): hornblende, biotite, K-feldspar (a second copy of the feldspar model, Or above 30%), "
    "muscovite, quartz and tridymite. MELTS's hornblende and biotite are simple models and seldom "
    "saturate; treat them as rough. Equilibrium paths end near the solidus.",
    "Not in the model: apatite, zircon, sphene, cordierite, garnet in granites (it is the mantle garnet), "
    "and the Na-rich amphiboles.",
    "Magma types: by SiO2 (basalt 45-52, basaltic andesite 52-57, andesite 57-63, dacite 63-69, rhyolite "
    "above 69 wt%), and the TAS name of Le Bas et al. (1986) from SiO2 and Na2O + K2O, volatile-free. "
    "As in MELTS, the liquid holds sodium too well, so evolved melts plot high in alkalis (the trachy- "
    "fields of TAS) where real ones are andesite or dacite.",
    "Iron: Fe3+/Fe2+ of the liquid from Kress and Carmichael (1991) at the chosen oxygen buffer "
    "(Frost 1991), or fixed by the starting Fe2O3 when the system is closed to oxygen.",
    "MnO, NiO, CoO and CO2 are not part of the liquid model and are left out of the calculation.",
    "Water: dissolved H2O uses the MELTS v1.0.2 standard state; a pure-water fluid (IAPWS-95) forms when "
    "the melt is saturated. At 1 bar any water leaves the melt at once; MORB at 1200 C holds about 3 wt% "
    "H2O at 1 kbar.",
    "Masses are grams per 100 g of starting material. In fractional crystallization the solids of each "
    "step are removed, and the solid masses shown are the totals removed so far.",
    "Melting: batch melting is equilibrium, computed by cooling the rock's composition from above its "
    "liquidus and shown from low to high temperature; it starts at the solidus, which Magma cannot cross "
    "yet. Fractional melting removes the melt above 1% of the residue at each step and pools it. "
    "Garnet is not in the model, so peridotite melting is right only in the spinel field (below about 2 GPa). "
    "The lowest-degree melts are too rich in Na2O, as in MELTS.",
]


def _fo2_arg(v):
    if v in (None, "", "closed"):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return str(v)


class State:
    """What both windows show.  Every change bumps `version` and wakes the event streams."""

    def __init__(self):
        self.cond = threading.Condition()
        self.version = 0
        self.status = "idle"            # idle | running | done | stopped | error
        self.message = ""
        self.progress = [0, 0]
        self.inputs: dict = {}
        self.steps: list[dict] = []
        self.run_id = 0
        self.source: dict | None = None
        self._stop = threading.Event()

    def snapshot(self) -> dict:
        return {"version": self.version, "status": self.status, "message": self.message,
                "progress": self.progress, "inputs": self.inputs, "steps": self.steps,
                "run_id": self.run_id,
                "notes": HP_NOTES if self.inputs.get("engine") == "hp" else NOTES, "app_version": __version__,
                "source": self.source}

    def changed(self):
        with self.cond:
            self.version += 1
            self.cond.notify_all()

    # ---- runs ----------------------------------------------------------------------------------
    def start(self, inputs: dict) -> None:
        if self.status == "running":
            raise ValueError("a run is already going; stop it first")
        bulk, kw = parse_inputs(inputs)
        self.source = None
        self._stop.clear()
        self.run_id += 1
        self.inputs = inputs
        self.steps = []
        self.status, self.message = "running", ""
        n = int(round(abs(kw["t_start"] - kw["t_end"]) / kw["step"])) + 1
        self.progress = [0, n]
        self.changed()
        threading.Thread(target=self._run, args=(bulk, kw), daemon=True).start()

    def stop(self) -> None:
        self._stop.set()

    def _run(self, bulk: Composition, kw: dict) -> None:
        removed: dict[str, float] = {}
        try:
            if kw["process"] == "melt":
                return self._run_melt(bulk, kw)
            if kw["process"] == "example":
                return self._run_example(bulk, kw)
            kw = {k: v for k, v in kw.items() if k not in ("process", "engine")}
            cool = hp_engine.iter_crystallize if self.inputs.get("engine") == "hp" else iter_crystallize
            for r in cool(bulk, **kw):
                step = step_json(r)
                if kw["mode"] == "fractional":
                    for s in step["solids"]:
                        removed[s["name"]] = removed.get(s["name"], 0.0) + s["grams"]
                        s["step_grams"] = s["grams"]
                        s["grams"] = removed[s["name"]]
                    for name, g in removed.items():
                        if not any(s["name"] == name for s in step["solids"]):
                            step["solids"].append({"name": name, "label": PHASE_LABEL.get(name, name),
                                                   "grams": g, "step_grams": 0.0, "x": {}, "numbers": {},
                                                   "formula": ""})
                self.steps.append(step)
                self.progress = [len(self.steps), self.progress[1]]
                self.changed()
                if self._stop.is_set():
                    self.status, self.message = "stopped", "Stopped."
                    break
            else:
                self.status = "done"
                if self.steps and self.steps[-1]["t"] > kw["t_end"] + 1e-6:
                    if self.inputs.get("engine") == "hp":
                        self.message = f"Ended at {self.steps[-1]['t']:.0f} °C: the melt is used up (the solidus)."
                    else:
                        self.message = (f"Ended at {self.steps[-1]['t']:.0f} °C: the melt is nearly gone (near the "
                                        "solidus), which Magma cannot follow further yet.")
                elif not self.steps:
                    self.message = "No melt at these conditions: the composition is below its solidus."
        except Exception as exc:          # shown in both windows
            self.status, self.message = "error", f"{type(exc).__name__}: {exc}"
        self.changed()

    def _run_melt(self, bulk: Composition, kw: dict) -> None:
        def progress(k, n):
            self.progress = [k, n]
            self.changed()
            if self._stop.is_set():
                raise _Stopped()
        args = dict(kw, mode="batch" if kw["mode"] == "equilibrium" else "fractional")
        args.pop("process")
        hp = args.pop("engine") == "hp"
        try:
            melt = hp_engine.iter_melt if hp else iter_melt
            for m in melt(bulk, on_progress=progress if (hp or args["mode"] == "batch") else None, **args):
                step = step_json(m.result)
                step["extracted_grams"] = m.extracted_grams
                step["pooled_grams"] = m.pooled_grams
                if m.pooled is not None:
                    step["pooled"] = {ox: m.pooled.get(ox) for ox in OXIDES if m.pooled.get(ox) > 0}
                    step["pooled_mg"] = m.pooled.mg_number
                self.steps.append(step)
                if args["mode"] == "fractional" and not hp:
                    self.progress = [len(self.steps), self.progress[1]]
                self.changed()
                if self._stop.is_set():
                    raise _Stopped()
            self.status = "done"
            if not self.steps:
                self.message = "The rock does not melt at these temperatures (below its solidus)."
            elif self.steps[0]["t"] > kw["t_start"] + 1e-6:
                self.message = (f"Melting starts near {self.steps[0]['t']:.0f} °C (the solidus)." if hp else
                                f"Melting starts near {self.steps[0]['t']:.0f} °C, the lowest temperature "
                                "Magma could follow (the solidus).")
        except _Stopped:
            self.status, self.message = "stopped", "Stopped."
        except Exception as exc:
            self.status, self.message = "error", f"{type(exc).__name__}: {exc}"
        self.changed()


    def _run_example(self, rock: Composition, kw: dict) -> None:
        """Mantle to magma: batch-melt the rock at the source, then cool the melt fractionally."""
        try:
            self.message = (f"Melting the mantle at {kw['t_source']:.0f} °C and {kw['p_source'] / 1000:g} kbar "
                            "(this takes a minute or two)...")
            self.changed()
            hp = kw["engine"] == "hp"
            if hp:
                src_result, primary = hp_engine.melt_source(rock, kw["t_source"], kw["p_source"],
                                                            fo2=kw["fo2_source"], phases=kw["phases"])
            else:
                src = melt_source(rock, kw["t_source"], kw["p_source"], fo2=kw["fo2_source"], phases=kw["phases"])
                src_result, primary = src.result, src.primary
            if self._stop.is_set():
                raise _Stopped()
            sj = step_json(src_result)
            sj["primary"] = {ox: primary.get(ox) for ox in OXIDES if primary.get(ox) > 0}
            sj["primary_mg"] = primary.mg_number
            self.source = sj
            self.message = ""
            self.changed()
            removed: dict[str, float] = {}
            if hp:
                path = hp_engine.iter_crystallize(primary_with_water(primary, kw["h2o"]), kw["t_source"],
                                                  kw["t_end"], kw["step"], kw["p_bar"], kw["fo2"], "fractional",
                                                  kw["phases"])
            else:
                path = iter_differentiate(primary, kw["t_source"], kw["t_end"], kw["step"], kw["p_bar"],
                                          kw["h2o"], kw["fo2"], kw["phases"])
            for r in path:
                step = step_json(r)
                for s in step["solids"]:
                    removed[s["name"]] = removed.get(s["name"], 0.0) + s["grams"]
                    s["step_grams"] = s["grams"]
                    s["grams"] = removed[s["name"]]
                for name, g in removed.items():
                    if not any(s["name"] == name for s in step["solids"]):
                        step["solids"].append({"name": name, "label": PHASE_LABEL.get(name, name),
                                               "grams": g, "step_grams": 0.0, "x": {}, "numbers": {},
                                               "formula": ""})
                self.steps.append(step)
                self.progress = [len(self.steps), self.progress[1]]
                self.changed()
                if self._stop.is_set():
                    raise _Stopped()
            self.status = "done"
            if self.steps and self.steps[-1]["t"] > kw["t_end"] + 1e-6:
                self.message = (f"Ended at {self.steps[-1]['t']:.0f} °C: the melt is used up (the solidus)." if hp else
                                f"Ended at {self.steps[-1]['t']:.0f} °C: the melt is nearly gone, which Magma "
                                "cannot follow further yet.")
        except _Stopped:
            self.status, self.message = "stopped", "Stopped."
        except Exception as exc:
            self.status, self.message = "error", f"{type(exc).__name__}: {exc}"
        self.changed()


class _Stopped(Exception):
    pass


def parse_inputs(d: dict) -> tuple[Composition, dict]:
    wt = {}
    for ox, v in (d.get("bulk") or {}).items():
        if ox not in OXIDES:
            raise ValueError(f"unknown oxide {ox}")
        if v in (None, ""):
            continue
        v = float(v)
        if v < 0:
            raise ValueError(f"{ox} is negative")
        if v > 0:
            wt[ox] = v
    if sum(wt.values()) <= 0:
        raise ValueError("enter a composition")
    if wt.get("SiO2", 0) <= 0:
        raise ValueError("the composition needs SiO2")
    phases = [p for p in d.get("phases", list(PHASES)) if p in PHASES]
    if not phases:
        raise ValueError("pick at least one mineral")
    kw = dict(t_start=float(d.get("t_start", 1250)), t_end=float(d.get("t_end", 1150)),
              step=float(d.get("step", 5)), p_bar=float(d.get("p_bar", 1)),
              fo2=_fo2_arg(d.get("fo2", "QFM")), mode=d.get("mode", "equilibrium"), phases=tuple(phases),
              process=d.get("process", "crystallize"), engine=d.get("engine", "hp"))
    if kw["engine"] not in ENGINES:
        raise ValueError("engine must be melts or hp")
    if kw["mode"] not in ("equilibrium", "fractional"):
        raise ValueError("mode must be equilibrium or fractional")
    if kw["process"] not in ("crystallize", "melt", "example"):
        raise ValueError("process must be crystallize, melt or example")
    if kw["process"] == "example":
        ex = d.get("example") or {}
        num = lambda k: float(ex.get(k, EXAMPLE_DEFAULTS[k]))
        kw.update(t_source=num("t_source"), p_source=num("p_source"), h2o=num("h2o"),
                  fo2_source=_fo2_arg(ex.get("fo2_source", "QFM")), mode="fractional")
        kw["t_start"] = kw["t_source"]
        if not 1100 <= kw["t_source"] <= 1700:
            raise ValueError("the melting temperature must be between 1100 and 1700 C")
        if not 1 <= kw["p_source"] <= 30000:
            raise ValueError("the melting pressure must be between 1 and 30000 bar")
        if not 0 <= kw["h2o"] <= 10:
            raise ValueError("the water added must be between 0 and 10 wt%")
        if not 600 <= kw["t_end"] < kw["t_source"]:
            raise ValueError("the melt must cool to between 600 C and the melting temperature")
    if kw["process"] == "crystallize" and not 600 <= kw["t_end"] <= kw["t_start"] <= 2000:
        raise ValueError("cooling: temperatures must run downward, between 600 and 2000 C")
    if kw["process"] == "melt" and not 600 <= kw["t_start"] <= kw["t_end"] <= 2000:
        raise ValueError("heating: temperatures must run upward, between 600 and 2000 C")
    if kw["step"] <= 0 or abs(kw["t_start"] - kw["t_end"]) / kw["step"] > 400:
        raise ValueError("the step must be positive, with at most 400 steps")
    if not 1 <= kw["p_bar"] <= 30000:
        raise ValueError("pressure must be between 1 and 30000 bar")
    return Composition(wt, name=d.get("name", "")), kw


def _numbers(name: str, x: dict) -> dict:
    """Composition numbers for the plots: olivine Fo, plagioclase An, pyroxene Wo and Mg#, spinel
    ulvospinel, Mg# and Cr#."""
    if name == "olivine":
        return {"Fo": 100 * x.get("forsterite", 0)}
    if name == "feldspar":
        return {"An": 100 * x.get("anorthite", 0)}
    if name == "alkali-feldspar":
        return {"Or": 100 * x.get("sanidine", 0), "An": 100 * x.get("anorthite", 0)}
    if name == "biotite":
        return {"Mg#": 100 * x.get("phlogopite", 0)}
    if name == "hornblende":
        return {"Mg#": 100 * (1 - x.get("ferropargasite", 0))}
    if name in ("clinopyroxene", "pigeonite", "orthopyroxene"):
        from ..thermo.pyroxene import r_from_x
        r = r_from_x([x[e] for e in PHASES[name].end_members])
        ca, fe = 1 - r[4] - r[5], r[0]
        mg = 1 - r[0] - r[3] - 0.5 * (r[1] + r[2] + r[4]) + r[5]
        return {"Wo": 100 * ca / (ca + mg + fe), "Mg#": 100 * mg / (mg + fe)}
    if name == "rhombohedral":
        return {"Ilm": 100 * x.get("ilmenite", 0), "Hem": 100 * x.get("hematite", 0)}
    if name in ("spinel", "titanomagnetite"):
        cr, hc, mt, sp, usp = (x.get(e, 0) for e in PHASES[name].end_members)
        mg, fe2, al, crr = sp, hc + cr + mt + 2 * usp, 2 * (hc + sp), 2 * cr
        out = {"Usp": 100 * usp, "Mg#": 100 * mg / (mg + fe2) if mg + fe2 > 0 else 0.0}
        if crr > 1e-6:
            out["Cr#"] = 100 * crr / (crr + al)
        return out
    return {}


def _solid_json(s) -> dict:
    if isinstance(s, hp_engine.HPPhase):
        return {"name": s.name, "label": PHASE_LABEL.get(s.name, s.name), "grams": s.grams, "x": s.x,
                "numbers": s.numbers, "formula": s.formula}
    return {"name": s.name, "label": PHASE_LABEL.get(s.name, s.name), "grams": s.grams, "x": s.x,
            "numbers": _numbers(s.name, s.x), "formula": PHASES[s.name].formula(np.array(list(s.x.values())))}


def step_json(r) -> dict:
    liq = r.liquid
    return {
        "t": r.t_c, "p": r.p_bar, "melt_pct": 100 * r.melt_fraction, "liquid_grams": r.liquid_grams,
        "liquid": {ox: liq.get(ox) for ox in OXIDES if liq.get(ox) > 0},
        "mg_number": liq.mg_number,
        "magma_type": magma_type(liq), "tas": tas_name(liq), "tas_xy": tas_point(liq),
        "solids": [_solid_json(s) for s in r.solids],
        # JSON has no infinity: a phase that cannot form at all (its affinity is -inf) is sent as null
        "affinities": {k: (float(v) if np.isfinite(v) else None) for k, v in r.affinities.items()},
    }


class Handler(BaseHTTPRequestHandler):
    server_version = "Magma"
    protocol_version = "HTTP/1.1"

    @property
    def state(self) -> State:
        return self.server.state            # type: ignore[attr-defined]

    def log_message(self, fmt, *args):
        pass

    def _allowed(self) -> bool:
        port = self.server.server_address[1]
        ok_hosts = {f"127.0.0.1:{port}", f"localhost:{port}", f"[::1]:{port}"}
        if self.headers.get("Host", "") not in ok_hosts:
            return False
        origin = self.headers.get("Origin")
        return origin is None or origin in {f"http://{h}" for h in ok_hosts}

    def _send(self, status: int, body: bytes, ctype: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, status: int = 200) -> None:
        self._send(status, json.dumps(obj).encode(), "application/json; charset=utf-8")

    def _error(self, status: int, message: str) -> None:
        self._json({"error": message}, status)

    def _file(self, p: Path) -> None:
        p = p.resolve()
        if STATIC not in p.parents or not p.is_file():
            return self._error(404, "not found")
        ctype = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype in ("application/javascript", "image/svg+xml"):
            ctype += "; charset=utf-8"
        self._send(200, p.read_bytes(), ctype)

    def do_GET(self) -> None:
        if not self._allowed():
            return self._error(HTTPStatus.FORBIDDEN, "forbidden")
        path = urlparse(self.path).path
        if path in PAGES:
            return self._file(STATIC / PAGES[path])
        if path.startswith("/static/"):
            return self._file(STATIC / path[len("/static/"):])
        if path == "/favicon.ico":                      # the Magma M, for browsers and pinned shortcuts that ask for it
            return self._file(STATIC / "magma.ico")
        if path == "/api/state":
            return self._json(self.state.snapshot())
        if path == "/api/presets":
            return self._json({"oxides": [o for o in OXIDES if o != "CO2"],
                               "presets": {k: {"wt": c.wt, "source": c.source, "name": c.name}
                                           for k, c in PRESETS.items()},
                               "phases": [{"name": n, "label": PHASE_LABEL.get(n, n)} for n in PHASES],
                               "tas": {"fields": TAS_FIELDS, "lines": TAS_LINES},
                               "bowen": {"discontinuous": BOWEN_DISCONTINUOUS, "continuous": BOWEN_CONTINUOUS,
                                         "late": BOWEN_LATE},
                               "example": EXAMPLE_DEFAULTS, "cards": CARDS})
        if path == "/events":
            return self._events()
        return self._error(404, "not found")

    def _events(self) -> None:
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "close")
        self.end_headers()
        self.close_connection = True
        seen = -1
        st = self.state
        try:
            while True:
                with st.cond:
                    if st.version == seen:
                        st.cond.wait(timeout=15)
                if st.version == seen:
                    self.wfile.write(b": keep-alive\n\n")
                else:
                    seen = st.version
                    self.wfile.write(b"data: " + json.dumps(st.snapshot()).encode() + b"\n\n")
                self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError, OSError):
            return

    def do_POST(self) -> None:
        if not self._allowed():
            return self._error(HTTPStatus.FORBIDDEN, "forbidden")
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_BODY:
            return self._error(413, "request too large")
        try:
            body = json.loads(self.rfile.read(n) or b"{}")
        except json.JSONDecodeError:
            return self._error(400, "bad JSON")
        path = urlparse(self.path).path
        try:
            if path == "/api/run":
                self.state.start(body)
                return self._json({"ok": True, "run_id": self.state.run_id})
            if path == "/api/stop":
                self.state.stop()
                return self._json({"ok": True})
        except (ValueError, TypeError) as exc:
            return self._error(400, str(exc))
        return self._error(404, "not found")


def make_server(port: int = DEFAULT_PORT) -> ThreadingHTTPServer:
    try:
        httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    except OSError:                      # port taken: let the system pick one
        httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    httpd.daemon_threads = True
    httpd.state = State()                # type: ignore[attr-defined]
    return httpd
