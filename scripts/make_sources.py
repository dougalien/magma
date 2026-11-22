"""Write docs/REFERENCES.md and the app's sources page (src/magma/app/static/sources.html) from one list.

Run from the repository root:  python scripts/make_sources.py
"""
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (authors, year, title, where, doi or "", used for)
SECTIONS = [
    ("Liquid and method (MELTS engine)", [
        ("Ghiorso, M.S., Sack, R.O.", "1995", "Chemical mass transfer in magmatic processes IV. A revised and internally consistent thermodynamic model for the interpolation and extrapolation of liquid-solid equilibria in magmatic systems at elevated temperatures and pressures",
         "Contributions to Mineralogy and Petrology 119, 197-212", "10.1007/BF00307281",
         "The silicate liquid: interaction parameters W (Table 4), liquid end members and standard states (Table A1 and appendix), anorthite and hematite adjustments."),
        ("Ghiorso, M.S.", "1985", "Chemical mass transfer in magmatic processes I. Thermodynamic relations and numerical algorithms",
         "Contributions to Mineralogy and Petrology 90, 107-120", "10.1007/BF00378254", "The equilibrium algorithm (equilibrium.py)."),
        ("Ghiorso, M.S.", "1994", "Algorithms for the estimation of phase stability in heterogeneous thermodynamic systems",
         "Geochimica et Cosmochimica Acta 58, 5489-5501", "10.1016/0016-7037(94)90245-3", "Adding and dropping phases by their affinities (equilibrium.py)."),
        ("Gualda, G.A.R., Ghiorso, M.S., Lemons, R.V., Carley, T.L.", "2012", "Rhyolite-MELTS: a modified calibration of MELTS optimized for silica-rich, fluid-bearing magmatic systems",
         "Journal of Petrology 53, 875-890", "10.1093/petrology/egr080", "The sanidine enthalpy adjustment (+3400 J) and the Early Bishop Tuff example."),
        ("Richet, P., Bottinga, Y., Denielou, L., Petitet, J.P., Tequi, C.", "1982", "Thermodynamic properties of quartz, cristobalite and amorphous SiO2: drop calorimetry measurements between 1000 and 1800 K and a review from 0 to 2000 K",
         "Geochimica et Cosmochimica Acta 46, 2639-2658", "10.1016/0016-7037(82)90383-0", "The amorphous-silica reference state of the SiO2 liquid component."),
        ("Robie, R.A., Hemingway, B.S., Fisher, J.R.", "1978", "Thermodynamic properties of minerals and related substances at 298.15 K and 1 bar (10^5 pascals) pressure and at higher temperatures",
         "U.S. Geological Survey Bulletin 1452", "10.3133/b1452", "The water term G_Robie in the dissolved-H2O standard state, as in Ghiorso and Sack (1995)."),
        ("Robie, R.A., Hemingway, B.S.", "1995", "Thermodynamic properties of minerals and related substances at 298.15 K and 1 bar (10^5 pascals) pressure and at higher temperatures",
         "U.S. Geological Survey Bulletin 2131", "10.3133/b2131", "The H and S of water at 298.15 K that put IAPWS-95 water on the minerals' energy scale (water_fluid.py)."),
    ]),
    ("Minerals (MELTS engine)", [
        ("Berman, R.G.", "1988", "Internally-consistent thermodynamic data for minerals in the system Na2O-K2O-CaO-MgO-FeO-Fe2O3-Al2O3-SiO2-TiO2-H2O-CO2",
         "Journal of Petrology 29, 445-522", "10.1093/petrology/29.2.445", "End-member H, S, heat capacity, volume, lambda and disorder terms (data/melts_minerals.tsv, minerals.py)."),
        ("Berman, R.G.", "1990", "Mixing properties of Ca-Mg-Fe-Mn garnets", "American Mineralogist 75, 328-344", "", "Garnet end members and mixing."),
        ("Berman, R.G., Koziol, A.M.", "1991", "Ternary excess properties of grossular-pyrope-almandine garnet and their influence in geothermobarometry",
         "American Mineralogist 76, 1223-1231", "", "Ternary garnet terms."),
        ("Sack, R.O., Ghiorso, M.S.", "1989", "Importance of considerations of mixing properties in establishing an internally consistent thermodynamic database: thermochemistry of minerals in the system Mg2SiO4-Fe2SiO4-SiO2",
         "Contributions to Mineralogy and Petrology 102, 41-68", "10.1007/BF01160190", "Olivine."),
        ("Hirschmann, M.", "1991", "Thermodynamics of multicomponent olivines and the solution properties of (Ni,Mg,Fe)2SiO4 and (Ca,Mg,Fe)2SiO4 olivines",
         "American Mineralogist 76, 1232-1248", "", "Olivine."),
        ("Elkins, L.T., Grove, T.L.", "1990", "Ternary feldspar experiments and thermodynamic models", "American Mineralogist 75, 544-559", "",
         "Plagioclase and alkali feldspar (Table 4)."),
        ("Salje, E.", "1985", "Thermodynamics of sodium feldspar I: order parameter treatment and strain induced coupling effects",
         "Physics and Chemistry of Minerals 12, 93-98", "10.1007/BF01046833", "The Landau ordering term of albite."),
        ("Carpenter, M.A.", "1988", "Thermochemistry of aluminium/silicon ordering in feldspar minerals",
         "In: Salje, E.K.H. (ed.) Physical Properties and Thermodynamic Behaviour of Minerals. NATO ASI Series C 225, Reidel, Dordrecht, 265-323", "",
         "The metastable C-1 anorthite adjustment, as in Ghiorso and Sack (1995)."),
        ("Sack, R.O., Ghiorso, M.S.", "1994", "Thermodynamics of multicomponent pyroxenes: III. Calibration of Fe2+(Mg)-1, TiAl2(MgSi2)-1, TiFe3+2(MgSi2)-1, AlFe3+(MgSi)-1, NaAl(CaMg)-1, Al2(MgSi)-1 and Ca(Mg)-1 exchange reactions between pyroxenes and silicate melts",
         "Contributions to Mineralogy and Petrology 118, 271-296", "10.1007/BF00306648",
         "Augite, pigeonite and orthopyroxene (cited in the code as 1994c, its letter in Ghiorso and Sack 1995)."),
        ("Sack, R.O., Ghiorso, M.S.", "1991", "An internally consistent model for the thermodynamic properties of Fe-Mg-titanomagnetite-aluminate spinels",
         "Contributions to Mineralogy and Petrology 106, 474-505", "10.1007/BF00321989", "Spinel, Cr-spinel to titanomagnetite."),
        ("Sack, R.O., Ghiorso, M.S.", "1991", "Erratum: An internally consistent model for the thermodynamic properties of Fe-Mg-titanomagnetite-aluminate spinels",
         "Contributions to Mineralogy and Petrology 107, 415", "10.1007/BF00325108", "The corrections to the spinel model."),
        ("Ghiorso, M.S., Evans, B.W.", "2008", "Thermodynamics of rhombohedral oxide solid solutions and a revision of the Fe-Ti two-oxide geothermometer and oxygen-barometer",
         "American Journal of Science 308, 957-1039", "10.2475/09.2008.01", "Ilmenite-hematite."),
    ]),
    ("Holland and Powell engine", [
        ("Holland, T.J.B., Powell, R.", "2011", "An improved and extended internally consistent thermodynamic dataset for phases of petrological interest, involving a new equation of state for solids",
         "Journal of Metamorphic Geology 29, 333-383", "10.1111/j.1525-1314.2010.00923.x",
         "End-member dataset ds636 (data/hp_ds636.tsv): heat capacity, modified Tait equation of state, Landau terms."),
        ("Holland, T.J.B., Powell, R.", "1996", "Thermodynamics of order-disorder in minerals: I. Symmetric formalism applied to minerals of fixed composition",
         "American Mineralogist 81, 1413-1424", "10.2138/am-1996-11-1214", "Bragg-Williams order-disorder of end members."),
        ("Holland, T.J.B., Powell, R.", "2003", "Activity-composition relations for phases in petrological calculations: an asymmetric multicomponent formulation",
         "Contributions to Mineralogy and Petrology 145, 492-501", "10.1007/s00410-003-0464-z", "The asymmetric (van Laar) excess terms."),
        ("White, R.W., Powell, R., Holland, T.J.B., Johnson, T.E., Green, E.C.R.", "2014", "New mineral activity-composition relations for thermodynamic calculations in metapelitic systems",
         "Journal of Metamorphic Geology 32, 261-286", "10.1111/jmg.12071", "Biotite, white mica, cordierite, chlorite."),
        ("Green, E.C.R., White, R.W., Diener, J.F.A., Powell, R., Holland, T.J.B., Palin, R.M.", "2016", "Activity-composition relations for the calculation of partial melting equilibria in metabasic rocks",
         "Journal of Metamorphic Geology 34, 845-869", "10.1111/jmg.12211", "Clinoamphibole."),
        ("Holland, T.J.B., Green, E.C.R., Powell, R.", "2018", "Melting of peridotites through to granites: a simple thermodynamic model in the system KNCFMASHTOCr",
         "Journal of Petrology 59, 881-900", "10.1093/petrology/egy048", "The igneous set: silicate liquid, aqueous fluid and the solid solutions."),
        ("Green, E.C.R., Holland, T.J.B., Powell, R., Weller, O.M., Riel, N.", "2025", "Corrigendum to: Melting of peridotites through to granites: a simple thermodynamic model in the system KNCFMASHTOCr",
         "Journal of Petrology 66, egae079", "10.1093/petrology/egae079",
         "Corrects the melt model of Holland et al. (2018). Magma follows MAGEMin 2.1.0; whether that release includes the correction has not been checked."),
        ("Pitzer, K.S., Sterner, S.M.", "1994", "Equations of state valid continuously from zero to extreme pressures for H2O and CO2",
         "Journal of Chemical Physics 101, 3111-3116", "10.1063/1.467624", "H2O and CO2 end members."),
    ]),
    ("Iron and oxygen", [
        ("Kress, V.C., Carmichael, I.S.E.", "1991", "The compressibility of silicate liquids containing Fe2O3 and the effect of composition, temperature, oxygen fugacity and pressure on their redox states",
         "Contributions to Mineralogy and Petrology 108, 82-92", "10.1007/BF00307328", "Melt Fe3+/Fe2+ at a given fO2 (eq. 7; redox.py)."),
        ("Frost, B.R.", "1991", "Introduction to oxygen fugacity and its petrologic importance",
         "In: Lindsley, D.H. (ed.) Oxide Minerals: Petrologic and Magnetic Significance. Reviews in Mineralogy 25, Mineralogical Society of America, 1-10",
         "10.1515/9781501508684-004", "The QFM, NNO, IW and MH buffers (Table 1; redox.py)."),
    ]),
    ("Water", [
        ("Wagner, W., Pruss, A.", "2002", "The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use",
         "Journal of Physical and Chemical Reference Data 31, 387-535", "10.1063/1.1461829", "Pure water and the water fluid (thermo/water.py, water_fluid.py)."),
        ("Haar, L., Gallagher, J.S., Kell, G.S.", "1984", "NBS/NRC Steam Tables: Thermodynamic and Transport Properties and Computer Programs for Vapor and Liquid States of Water in SI Units",
         "Hemisphere, Washington, DC", "", "The water equation of state MELTS uses; named for comparison only."),
        ("Delany, J.M., Helgeson, H.C.", "1978", "Calculation of the thermodynamic consequences of dehydration in subducting oceanic crust to 100 kb and >800 C",
         "American Journal of Science 278, 638-686", "10.2475/ajs.278.5.638", "The MELTS high-pressure water correction; named for comparison only."),
    ]),
    ("Compositions and names", [
        ("Gale, A., Dalton, C.A., Langmuir, C.H., Su, Y., Schilling, J.-G.", "2013", "The mean composition of ocean ridge basalts",
         "Geochemistry, Geophysics, Geosystems 14, 489-518", "10.1029/2012GC004334", "The ALL MORB preset (Table 1, arithmetic mean)."),
        ("Rudnick, R.L., Gao, S.", "2003", "Composition of the continental crust",
         "In: Rudnick, R.L. (ed.) The Crust. Treatise on Geochemistry 3, Elsevier, 1-64", "10.1016/B0-08-043751-6/03016-4",
         "The upper, lower and bulk continental crust presets."),
        ("Takahashi, E.", "1986", "Melting of a dry peridotite KLB-1 up to 14 GPa: implications on the origin of peridotitic upper mantle",
         "Journal of Geophysical Research 91(B9), 9367-9382", "10.1029/JB091iB09p09367", "The KLB-1 rock."),
        ("Davis, F.A., Tangeman, J.A., Tenner, T.J., Hirschmann, M.M.", "2009", "The composition of KLB-1 peridotite",
         "American Mineralogist 94, 176-180", "10.2138/am.2009.2984", "The KLB-1 preset (Table 1, recommended composition)."),
        ("Le Bas, M.J., Le Maitre, R.W., Streckeisen, A., Zanettin, B.", "1986", "A chemical classification of volcanic rocks based on the total alkali-silica diagram",
         "Journal of Petrology 27, 745-750", "10.1093/petrology/27.3.745", "The TAS names (classify.py)."),
        ("Wieser, M.E., Berglund, M.", "2009", "Atomic weights of the elements 2007 (IUPAC Technical Report)",
         "Pure and Applied Chemistry 81, 2131-2156", "10.1351/PAC-REP-09-08-03", "Atomic weights (composition.py)."),
    ]),
    ("Checks", [
        ("Hirose, K., Kushiro, I.", "1993", "Partial melting of dry peridotites at high pressures: determination of compositions of melts segregated from peridotite using aggregates of diamond",
         "Earth and Planetary Science Letters 114, 477-489", "10.1016/0012-821X(93)90077-M", "Comparison for KLB-1 melting (clinopyroxene out)."),
    ]),
    ("Software", [
        ("ENKI project (M.S. Ghiorso and others)", "n.d.", "ThermoEngine: the MELTS thermodynamic models and Python interface",
         "https://gitlab.com/ENKI-portal/ThermoEngine (AGPL-3)", "",
         "Used to check numbers, for parameter values that are hard to read in the scanned papers (CpxBerman.m, OpxBerman.m, SpinelBerman.m, RhombohedralBerman.m, LiquidMelts.m, LiquidMeltsH2O.m, BiotiteBerman.m, HornblendeBerman.m), and for the MELTS example basalt and Early Bishop Tuff compositions. Numbers only; no ThermoEngine code is in Magma."),
        ("Riel, N., Kaus, B.J.P., Green, E.C.R., Berlie, N.", "2022", "MAGEMin, an efficient Gibbs energy minimizer: application to igneous systems",
         "Geochemistry, Geophysics, Geosystems 23, e2022GC010427", "10.1029/2022GC010427",
         "The Holland and Powell parameter values as distributed in MAGEMin (GPL-3), and the reference equilibria in the tests. Numbers only; Magma's code is its own."),
    ]),
]


def md() -> str:
    out = ["# References", "",
           "Every paper, data set and program Magma's models and presets come from. The code and data",
           "files cite these by author and year. The same list is in the app (Sources link).", ""]
    for title, refs in SECTIONS:
        out += [f"## {title}", ""]
        for a, y, t, w, doi, use in refs:
            link = f" [doi:{doi}](<https://doi.org/{doi}>)" if doi else ""
            out.append(f"- {a} ({y}). {t}. *{w}*.{link}")
            out.append(f"  Used for: {use}")
        out.append("")
    return "\n".join(out)


def html_list() -> str:
    parts = []
    for title, refs in SECTIONS:
        items = []
        for a, y, t, w, doi, use in refs:
            link = (f' <a href="https://doi.org/{escape(doi)}" target="_blank" rel="noopener">doi:{escape(doi)}</a>'
                    if doi else "")
            items.append(f"<li>{escape(a)} ({y}). {escape(t)}. <i>{escape(w)}</i>.{link}"
                         f'<div class="muted small">Used for: {escape(use)}</div></li>')
        parts.append(f"<h2>{escape(title)}</h2><ol>{''.join(items)}</ol>")
    return "".join(parts)


CITE = ('<h2>How to cite</h2><p>If you use Magma in a study or a class, please cite it (the repository\'s '
        'CITATION.cff gives the reference), and cite the models and data it used for your result, listed below.</p>')

if __name__ == "__main__":
    (ROOT / "docs" / "REFERENCES.md").write_text(md(), encoding="utf-8")
    page = ROOT / "src" / "magma" / "app" / "static" / "sources.html"
    html = page.read_text(encoding="utf-8")
    head, rest = html.split('<main class="sources">', 1)
    intro = ("<p>Magma is an original implementation of the MELTS models and of the Holland and Powell igneous "
             "models, from the published equations and parameters. These are the papers, data and programs its "
             "models and starting compositions come from.</p>")
    page.write_text(head + '<main class="sources">\n' + intro + "\n" + CITE + "\n" + html_list() + "\n</main>\n</body>\n</html>\n",
                    encoding="utf-8")
    print("wrote docs/REFERENCES.md and sources.html")
