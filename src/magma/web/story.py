"""The written account of a run, as the desktop results window gives it (ported from results.js).

Every function takes the run snapshot (``magma.app.server.State.snapshot()``) and returns plain
sentences; nothing here computes thermodynamics.
"""
from __future__ import annotations

import re

from ..classify import BOWEN_DISCONTINUOUS, BOWEN_LATE


def _g(p: dict) -> float:
    return p["step_grams"] if "step_grams" in p else p["grams"]


def _order(steps, use_step_grams=True) -> list[dict]:
    order: list[dict] = []
    for s in steps:
        for p in s["solids"]:
            g = _g(p) if use_step_grams else p["grams"]
            if g > 1e-6 and not any(o["n"] == p["name"] for o in order):
                order.append({"n": p["name"], "t": s["t"], "l": p["label"]})
    return order


def type_steps(steps):
    """The first step of each magma type the melt passes through."""
    return [s for i, s in enumerate(steps) if i == 0 or s["magma_type"] != steps[i - 1]["magma_type"]]


def runaway_note(steps) -> str | None:
    """The sodium-rich runaway of the MELTS liquid (silica falls while alkalis pass 10 wt%)."""
    best = float("-inf")
    for s in steps:
        si, alk = s["tas_xy"]
        if si < best - 1.5 and alk > 10:
            return (f"Below about {s['t']:.0f} °C the melt's silica falls while its Na2O + K2O climbs "
                    f"({alk:.1f} wt% here). That is the MELTS liquid model holding too much sodium, which "
                    "ThermoEngine's own MELTS does as well; it is not a real trend, so read the path no "
                    "further than this.")
        best = max(best, si)
    return None


def _melt_story(st) -> list[str]:
    steps, out = st["steps"], []
    first, last = steps[0], steps[-1]
    out.append(f"At {first['t']:.0f} °C the rock is {first['melt_pct']:.1f}% melt; its residue is "
               f"{', '.join(p['label'].lower() for p in first['solids'])}.")
    gone: list[dict] = []
    for i in range(1, len(steps)):
        for p in steps[i - 1]["solids"]:
            if (p["grams"] > 1e-3 and not any(q["name"] == p["name"] and q["grams"] > 1e-6 for q in steps[i]["solids"])
                    and not any(g["n"] == p["name"] for g in gone)):
                gone.append({"n": p["name"], "l": p["label"], "t": steps[i]["t"]})
    if gone:
        out.append("Used up: " + ", ".join(f"{g['l'].lower()} by {g['t']:.0f} °C" for g in gone) + ".")
    if st["inputs"].get("mode") == "fractional":
        pm = last.get("pooled")
        s = f"By {last['t']:.0f} °C, {last.get('pooled_grams', 0):.1f} g of melt has been removed from 100 g of rock"
        if pm:
            s += (f"; the pooled melt has {pm.get('SiO2', 0):.1f} wt% SiO2, {pm.get('MgO', 0):.1f} wt% MgO "
                  f"and Mg# {last['pooled_mg']:.0f}.")
        out.append(s if pm else s + ".")
    else:
        out.append(f"At {last['t']:.0f} °C the rock is {last['melt_pct']:.1f}% melt, with "
                   f"{last['liquid'].get('SiO2', 0):.1f} wt% SiO2, {last['liquid'].get('MgO', 0):.1f} wt% MgO "
                   f"and Mg# {last['mg_number']:.0f}.")
    return out


def _example_story(st) -> list[str]:
    src, steps, i = st["source"], st["steps"], st["inputs"]
    ex, out = i.get("example") or {}, []
    if not src:
        out.append((st["message"] or "Stage 1: melting the mantle rock...") if st["status"] == "running"
                   else "Stage 1 did not finish.")
        return out
    pr = src["primary"]
    out.append(f"Stage 1, melting the mantle: at {src['t']:.0f} °C and {src['p'] / 1000:.0f} kbar the "
               f"{i.get('name') or 'rock'} is {src['melt_pct']:.1f}% melt, and the solid left behind is "
               f"{', '.join(p['label'].lower() for p in src['solids'])}. That melt is the primary magma: "
               f"{pr.get('SiO2', 0):.1f} wt% SiO2, {pr.get('MgO', 0):.1f} wt% MgO, Mg# {src['primary_mg']:.0f}, "
               f"a {src['magma_type']} (TAS: {src['tas']}).")
    if not steps:
        return out
    last = steps[-1]
    water = f"with {ex.get('h2o')} wt% water added" if float(ex.get("h2o") or 0) > 0 else "dry"
    out.append(f"Stage 2, cooling in the crust: the melt rises to {i['p_bar'] / 1000:.1f} kbar, {water}, at the "
               f"{i['fo2']} oxygen buffer, and the crystals are taken out as they form (fractional crystallization).")
    order = _order(steps)
    if order:
        out.append("Minerals in the order they appear: "
                   + ", ".join(f"{o['l'].lower()} ({o['t']:.0f} °C)" for o in order)
                   + ". See Bowen's series for how this compares with Bowen's order.")
    ts = type_steps(steps)
    if len(ts) > 1:
        parts = [ts[0]["magma_type"]] + [f"{s['magma_type']} at {s['t']:.0f} °C ({s['liquid_grams']:.0f} g of melt left)"
                                         for s in ts[1:]]
        out.append("The melt changes from " + ", to ".join(parts) + ".")
    out.append(f"At {last['t']:.0f} °C, {last['liquid_grams']:.1f} g of melt is left from 100 g of primary magma: "
               f"{last['tas_xy'][0]:.1f} wt% SiO2 and {last['tas_xy'][1]:.1f} wt% Na2O + K2O, a "
               f"{last['magma_type']} by silica (TAS: {last['tas']}).")
    rn = runaway_note(steps)
    if rn:
        out.append(rn)
    if re.search(r"trachy|phono|tephri|foid", last["tas"]):
        out.append("The TAS name is alkaline because the model's liquid keeps too much sodium (a known failing of "
                   "MELTS as well); a real melt with this silica would be the subalkaline type named by silica.")
    grew = {o["n"] for o in order}
    late = [n for n in ("hornblende", "biotite", "alkali-feldspar", "muscovite", "quartz") if n not in grew]
    if late:
        out.append("Late minerals of Bowen's series that did not form in this run: "
                   + ", ".join("K-feldspar" if n == "alkali-feldspar" else n for n in late)
                   + ". Magma has them (as MELTS does), but MELTS's hornblende and biotite are simple models that "
                   "seldom saturate, and muscovite needs higher pressure (above about 3 kbar).")
    if last["magma_type"] != "rhyolite":
        out.append("To follow Bowen's series to the end, run the \"Arc andesite\" or \"Dacite\" card: a crustal "
                   "andesite, richer in potassium, goes on to rhyolite, then grows quartz, biotite and K-feldspar "
                   "together at about 720 °C, the granite minimum.")
    return out


def story(st) -> list[str]:
    """What happened in the run, in sentences."""
    proc = st["inputs"].get("process")
    if proc == "example":
        return _example_story(st)
    if not st["steps"]:
        return []
    if proc == "melt":
        return _melt_story(st)
    steps, out = st["steps"], []
    first = next((s for s in steps if any(p["grams"] > 1e-6 for p in s["solids"])), None)
    t0, last = steps[0]["t"], steps[-1]
    if not first:
        out.append(f"No mineral crystallizes between {t0:.0f} and {last['t']:.0f} °C; the melt stays fully liquid.")
    else:
        out.append(f"The first crystals appear at {first['t']:.0f} °C: "
                   + " and ".join(f"{p['label'].lower()} ({p['formula']})" for p in first["solids"]) + ".")
        order = _order(steps, use_step_grams=False)
        if len(order) > 1:
            out.append("Order of crystallization: "
                       + ", then ".join(f"{o['l'].lower()} at {o['t']:.0f} °C" for o in order) + ".")
        frac = (f" ({last['liquid_grams']:.1f} g of melt left from 100 g)."
                if st["inputs"].get("mode") == "fractional" else ".")
        out.append(f"At {last['t']:.0f} °C the system is {last['melt_pct']:.1f}% melt by mass{frac} The melt has "
                   f"{last['liquid'].get('SiO2', 0):.1f} wt% SiO2 and Mg# {last['mg_number']:.0f}.")
    if first and "clinopyroxene" not in (st["inputs"].get("phases") or []):
        out.append("Clinopyroxene was not allowed in this run.")
    rn = runaway_note(steps)
    if rn:
        out.append(rn)
    return out


def bowen(steps) -> dict:
    """Bowen's reaction series with what this run grew: {column: [(title, text, formed)]}."""
    first: dict[str, dict] = {}
    for s in steps:
        for p in s["solids"]:
            if _g(p) > 1e-6 and p["name"] not in first:
                first[p["name"]] = {"t": s["t"], "label": p["label"]}

    def box(title, names):
        seen = sorted((n for n in names if n in first), key=lambda n: -first[n]["t"])
        if not names:
            return (title, "not in Magma's model yet", False)
        if not seen:
            return (title, "did not form in this run", False)
        return (title, "; ".join(f"{first[n]['label'].lower()} from {first[n]['t']:.0f} °C" for n in seen), True)

    plag = [(s["t"], p) for s in steps for p in s["solids"]
            if p["name"] == "feldspar" and _g(p) > 1e-6 and "An" in (p.get("numbers") or {})]
    if not plag:
        cont = [box("Ca-rich plagioclase", ["feldspar"]), ("Na-rich plagioclase", "did not form in this run", False)]
    else:
        (ta, a), (tb, b) = plag[0], plag[-1]
        an_b = b["numbers"]["An"]
        cont = [("Ca-rich plagioclase", f"An{a['numbers']['An']:.0f} from {ta:.0f} °C", True),
                ("Na-rich plagioclase", f"An{an_b:.0f} at {tb:.0f} °C" if an_b < 50 else f"still An{an_b:.0f} at {tb:.0f} °C",
                 an_b < 50)]
    others = [n for n in ("spinel", "titanomagnetite", "rhombohedral", "garnet", "water") if n in first]
    return {"discontinuous": [box(t, n) for t, n in BOWEN_DISCONTINUOUS], "continuous": cont,
            "late": [box(t, n) for t, n in BOWEN_LATE],
            "also": ", ".join(f"{first[n]['label'].lower()} from {first[n]['t']:.0f} °C" for n in others)}
