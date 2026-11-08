"""Chemical formulas with real subscripts, as the desktop app's common.js writes them.

SiO2, Fe2O3, Na0.21Ca0.79Al1.79Si2.21O8, Fe2+, and end-member shorthand (Fo85, An60, Ilm69) get
subscripts; other words (QFM, Mg#, NNO+1, KLB-1, °C) are left alone, so the same function serves prose.
"""
from __future__ import annotations

import html
import re

FORMULA = re.compile(r"^(?:[A-Z][a-z]?(?:[1-9][+-])?(?:\d+(?:\.\d+)?)?|\(OH\)(?:\d+)?)+$")
PART = re.compile(r"([A-Z][a-z]?|\(OH\))([1-9][+-])?(\d+(?:\.\d+)?)?")
END_MEMBER = re.compile(r"^([A-Z][a-z]{2,3})(\d+(?:\.\d+)?)$")
WORD = re.compile(r"^([(\[“\"']*)(.*?)([)\].,;:!?”\"']*)$")
SUB = str.maketrans("0123456789+-", "₀₁₂₃₄₅₆₇₈₉₊₋")


def parts(text: str) -> list[tuple[str, str]]:
    """(text, kind) pieces; kind is '', 'sub' or 'sup'."""
    out: list[tuple[str, str]] = []

    def plain(t):
        if t:
            out.append((t, ""))
    for word in re.split(r"(\s+)", str(text)):
        m = WORD.match(word)
        pre, core, post = m.group(1), m.group(2), m.group(3)
        em = END_MEMBER.match(core)
        if em:
            plain(pre), plain(em.group(1)), out.append((em.group(2), "sub")), plain(post)
            continue
        if (not core or not FORMULA.match(core) or not re.search(r"\d|[+-]$", core)
                or re.match(r"^[A-Z][a-z]?[+-]?$", core)):
            plain(word)
            continue
        plain(pre)
        for el, charge, n in PART.findall(core):
            plain(el)
            if charge:
                out.append((charge, "sup"))
            if n:
                out.append((n, "sub"))
        plain(post)
    return out


def chem_html(text: str, escape: bool = False) -> str:
    """Formulas as <sub>/<sup> HTML, for st.markdown(unsafe_allow_html=True) and Plotly labels."""
    e = html.escape if escape else (lambda t: t)
    return "".join(f"<{k}>{e(t)}</{k}>" if k else e(t) for t, k in parts(text))


def chem_text(text: str) -> str:
    """Formulas with Unicode subscript digits, for places that take plain text only (labels, table headers).
    Only whole-number subscripts change; a formula with decimal amounts is left as it is."""
    out = []
    for word in re.split(r"(\s+)", str(text)):
        ps = parts(word)
        if any(k == "sub" and "." in t for t, k in ps):
            out.append(word)
        else:
            out.append("".join(t.translate(SUB) if k == "sub" else t for t, k in ps))
    return "".join(out)
