from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

DARK = RGBColor(0x2B, 0x2D, 0x42)   # basalt
GARNET = RGBColor(0x7A, 0x1F, 0x2B)
EMBER = RGBColor(0xD9, 0x77, 0x5C)
INK = RGBColor(0x1F, 0x21, 0x2B)
MUTED = RGBColor(0x6B, 0x6E, 0x7E)
ROCK = RGBColor(0xE9, 0xE4, 0xDC)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
HEAD, BODY = "Cambria", "Calibri"

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]


def box(slide, x, y, w, h, text, size=16, bold=False, color=INK, font=BODY,
        align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, italic=False, spacing=4):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    lines = text if isinstance(text, list) else [text]
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(spacing)
        runs = line if isinstance(line, list) else [(line, {})]
        for rt, opt in runs:
            r = p.add_run()
            r.text = rt
            f = r.font
            f.name = opt.get("font", font)
            f.size = Pt(opt.get("size", size))
            f.bold = opt.get("bold", bold)
            f.italic = opt.get("italic", italic)
            f.color.rgb = opt.get("color", color)
    return tb


def rect(slide, x, y, w, h, fill, line=None, shape=MSO_SHAPE.RECTANGLE):
    s = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line
        s.line.width = Pt(1)
    s.shadow.inherit = False
    return s


import json
EQ_SIZES = json.load(open("eq/sizes.json"))


def eq(slide, name, x, y, scale=2.0):
    """ Place a LaTeX-rendered equation (600 dpi PNG) at native 10pt * scale """
    w, h = EQ_SIZES[name]
    slide.shapes.add_picture(f"eq/{name}.png", Inches(x), Inches(y), width=Inches(w / 600 * scale))
    return h / 600 * scale


from pptx.oxml.ns import qn
from lxml import etree


def arrow(slide, x1, y1, x2, y2, color=MUTED, both=False, width=1.75):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1), Inches(x2), Inches(y2))
    c.line.color.rgb = color
    c.line.width = Pt(width)
    ln = c.line._get_or_add_ln()
    etree.SubElement(ln, qn("a:tailEnd")).set("type", "triangle")
    if both:
        etree.SubElement(ln, qn("a:headEnd")).set("type", "triangle")
    return c


def node(slide, x, y, w, h, head, body, fill=WHITE, line=DARK, head_color=DARK, body_color=INK, dashed=False):
    r = rect(slide, x, y, w, h, fill, line=line, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    r.adjustments[0] = 0.08
    r.line.width = Pt(1.25)
    if dashed:
        from pptx.enum.dml import MSO_LINE_DASH_STYLE
        r.line.dash_style = MSO_LINE_DASH_STYLE.DASH
    box(slide, x + 0.15, y + 0.08, w - 0.3, 0.3, head, size=13, bold=True, color=head_color)
    box(slide, x + 0.15, y + 0.38, w - 0.3, h - 0.45, body, size=10.5, color=body_color, spacing=1)
    return r


def tag(slide, x, y, w, text, color=MUTED, align=PP_ALIGN.LEFT):
    box(slide, x, y, w, 0.3, text, size=9.5, italic=True, color=color, align=align)


def title(slide, text, dark=False):
    box(slide, 0.6, 0.5, 12.1, 0.8, text, size=30, bold=True, font=HEAD,
        color=WHITE if dark else DARK)


def footer(slide, dark=False):
    box(slide, 0.6, 7.0, 8, 0.3, "pydfnworks.thermal  ·  dfnWorks, Los Alamos National Laboratory",
        size=10, color=RGBColor(0xC8, 0xC9, 0xD4) if dark else MUTED)



import csv
DATA = "../../examples/thermal_dfn_dataset/output/dataset"  # run the example first
rows = list(csv.DictReader(open(f"{DATA}/manifest.csv")))
by_sp = {}
for r in rows:
    by_sp.setdefault(float(r["spacing"]), []).append(r)

# ---------------------------------------------------------------- 1 title
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = "A DFN-based heat-extraction dataset: dfnWorks generates networks with parallel hydraulic fractures at a chosen spacing, graph flow and transport give every flow path, and an analytic fracture-rock exchange kernel turns the paths into produced temperature and heat rate. Each sample is one HDF5 file for the ML team."
rect(s, 0, 0, 13.333, 7.5, DARK)
rect(s, 0.6, 2.2, 0.12, 2.1, EMBER)
box(s, 0.95, 2.1, 11, 1.6, "A DFN heat-extraction dataset: produced temperature from fracture networks", size=40, bold=True, font=HEAD, color=WHITE)
box(s, 0.95, 3.75, 11, 0.6, "Parallel hydraulic fractures at a chosen spacing, graph flow and transport, and an analytic fracture-rock exchange kernel, one HDF5 sample at a time", size=18, color=RGBColor(0xE8, 0xA8, 0x7C), italic=True)
box(s, 0.95, 5.6, 11, 0.8, ["pydfnworks.thermal  ·  branch feature/thermal-mpf-surrogate (hymanjd/dfnWorks)", "Jeffrey Hyman, EES-16, Los Alamos National Laboratory  ·  November 2026"], size=14, color=RGBColor(0xC8, 0xC9, 0xD4))

# ---------------------------------------------------------------- 2 pipeline
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """The pipeline, left to right. dfnGen builds a network with two ingredients: parallel rectangular hydraulic fractures placed by the dataset helper at spacing D, spanning the block, and two stochastic natural fracture sets. Graph flow solves steady flow from the inflow face to the outflow face, so water enters every hydraulic fracture and flows along it, with the natural network cross-connecting the planes. Graph particle tracking, purely advective, records for every particle its advective travel time and the exchange integral beta, the sum over its edges of twice the residence time over the aperture. That pair is all the thermal model needs: the slab kernel with half-width D over 2 converts each path into a cooling curve, and the flux-weighted sum over particles is the cooled fraction of produced water. From that, T_out of t and Q dot of t. write_sample stores the inputs, the curves and the particle pairs in one HDF5 file and appends a manifest row. Nothing in the chain is meshed and no PFLOTRAN run is needed, so a sample costs seconds."""
title(s, "Pipeline: from a DFN to one dataset sample")
LIGHT = RGBColor(0xF4, 0xEF, 0xE8); TINT = RGBColor(0xFB, 0xF1, 0xEC)
for x, w, text in ((0.6, 3.0, "dfnWorks, existing"), (4.4, 3.6, "pydfnworks.thermal, new"), (9.2, 3.5, "Sample for the ML team")):
    box(s, x, 1.35, w, 0.3, text, size=12, bold=True, color=GARNET)
ax, aw, ah = 0.6, 3.0, 0.95
node(s, ax, 1.75, aw, ah, "dfnGen", "parallel hydraulic fractures at spacing D (user rectangles) + natural TPL sets; seed", fill=LIGHT)
node(s, ax, 3.0, aw, ah, "dfnGraph flow", "run_graph_flow, inflow to outflow face, at the drop that gives the target m_dot: flow graph G, dQ", fill=LIGHT)
node(s, ax, 4.25, aw, ah, "dfnGraph transport", "run_graph_transport, advective: per particle tau and beta = sum 2 tau_i / b_i", fill=LIGHT)
for y in (2.7, 3.95):
    arrow(s, ax + aw / 2, y, ax + aw / 2, y + 0.3)
bx, bw = 4.4, 3.6
node(s, bx, 1.75, bw, 0.95, "Inputs", "spacing D, seed, permeabilities, target m_dot, T_res, T_inj, rock and fluid properties, t_eval", fill=WHITE)
node(s, bx, 3.0, bw, 0.95, "Dataset helpers", "add_hydraulic_fractures: floor(L/D) rectangles spanning the block; graph_flow_for_rate: the drop that delivers m_dot", fill=TINT, line=GARNET, head_color=GARNET)
node(s, bx, 4.25, bw, 1.25, "heat_extraction()", "slab kernel, B = D/2, per particle; F(t) = flux-weighted cooled fraction; T_out(t), Q(t) with m_dot from the graph flow", fill=TINT, line=GARNET, head_color=GARNET)
node(s, bx, 5.8, bw, 0.8, "write_sample(), append_manifest()", "HDF5 + manifest.csv", fill=TINT, line=GARNET, head_color=GARNET)
arrow(s, bx + bw / 2, 2.7, bx + bw / 2, 3.0)
arrow(s, bx + bw / 2, 5.5, bx + bw / 2, 5.8)
arrow(s, bx, 3.47, ax + aw, 2.3)
tag(s, ax + aw + 0.05, 2.55, 0.8, "rectangles")
arrow(s, ax + aw, 4.72, bx, 4.72)
tag(s, ax + aw + 0.05, 4.42, 0.8, "tau, beta")
arrow(s, ax + aw, 3.47, bx - 0.02, 4.3)
tag(s, ax + aw + 0.05, 3.55, 0.8, "G")
cx, cw = 9.2, 3.5
node(s, cx, 1.75, cw, 2.5, "spacing_<D>_seed_<s>.h5", ["inputs/: spacing, seed, k, pressures, T, properties, fracture count, p32, dQ (attrs)", "curves/: t_eval, F, T_out, Q", "particles/: advect_time, beta", "attrs: m_dot, flow_rate, n_exit, n_stuck"], fill=LIGHT)
node(s, cx, 4.5, cw, 1.1, "manifest.csv", "one row per sample: inputs and headline outputs (pressure drop, Q(0), Q(30 yr), T_out(30 yr))", fill=LIGHT)
arrow(s, bx + bw, 6.2, cx + cw / 2, 6.2)
arrow(s, cx + cw / 2, 6.2, cx + cw / 2, 5.6)
box(s, 0.6, 6.75, 12.1, 0.3, "No meshing, PFLOTRAN or FEHM: a sample takes seconds, so spacings, seeds, permeabilities and intensities can be swept freely.", size=11, italic=True, color=MUTED)
footer(s)

# ---------------------------------------------------------------- 3 physics
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """The physics is heat treated as a solute. Dividing the fluid energy balance by the fluid heat capacity puts the thermal problem in the solute form with two substitutions: the matrix porosity becomes the rock-to-fluid volumetric heat capacity ratio, about one half for granite and water, and the matrix diffusivity becomes the rock thermal diffusivity. Graph particle tracking already accumulates, for every particle, the exchange integral beta, twice the residence time over the aperture summed along the path; it does this whether or not TDRW is on. For rock slabs of one half-width B everywhere, the thermal retention of a whole path has a closed Laplace transform, the compound-Poisson limit of the finite-slab TDRW model, and that is the Tang, Frind and Sudicky kernel, the same mathematics as Gringarten's parallel-fracture solution. Rescaled, it is the slab drawdown T_wD evaluated at the particle's own time scale and its own beta_G, so every particle contributes a known cooling curve starting at its advective arrival time. The flux-weighted sum over particles is the cooled fraction F of t, and the produced temperature and heat rate follow, with the mass flow rate the graph flow actually delivers. Two consequences: the result is deterministic given the particle set, with no extra Monte Carlo noise from thermal sampling, and the flow partition among fractures, which the old equal-flow surrogate assumed away, comes straight from the DFN."""
title(s, "Heat as a solute: each flow path carries its own cooling curve")
cols = ((0.6, (("Thermal to TDRW parameter mapping", "map", 1.5),
               ("Exchange integral the tracker accumulates along a path", "betaP", 1.5),
               ("Thermal retention of a path between slabs of half-width B (Laplace domain)", "kern", 1.5))),
        (6.9, (("Cooled fraction of produced water, flux-weighted over particles", "Fsum", 1.5),
               ("Per-particle time scale and slab parameter", "scales", 1.5),
               ("Produced temperature", "Tout1", 1.5),
               ("Heat extraction rate", "Qdot", 1.5))))
for x, items in cols:
    y = 1.5
    for head, name, scale in items:
        box(s, x, y, 5.9, 0.3, head, size=13, bold=True, color=GARNET)
        y += 0.34 + eq(s, name, x + 0.2, y + 0.32, scale) + 0.2
box(s, 0.6, 6.3, 12.0, 0.5, "B -> inf recovers Lauwerier; tau_i, b_i are edge residence time and aperture; w_p are flux weights (equal for flux-weighted release); m_dot = rho_f times the graph flow rate.", size=11, color=MUTED)
footer(s)

# ---------------------------------------------------------------- 4 kernel + verification
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """The slab drawdown kernel is inverted numerically with Gaver-Stehfest, 18 terms, vectorized over particles and times, using the routine the dfnWorks TDRW models already share. The curves show the family: close spacing, small beta_G, holds the full rate and then collapses once the slab is exhausted; wide spacing follows the single-fracture erfc decline. The verification has six parts and all pass. The erfc limit is matched exactly and the piston limit closes the heat balance to two parts in ten thousand. A finite-difference solution of the slab conduction and fracture advection equations, which never sees the Laplace transform, converges to the inverted kernel on grid refinement. A single synthetic path with infinite rock reproduces Lauwerier's closed form to machine precision, and with finite slabs reproduces the kernel. Twenty thousand draws from the graph-transport infinite-matrix TDRW sampler give an empirical CDF within one percent of the analytic breakthrough for the same edge, which ties the new kernel to the existing particle code. Superposition and the HDF5 round trip are checked as well."""
title(s, "Slab kernel: Stehfest inversion, verified six ways")
s.shapes.add_picture("drawdown_curves.png", Inches(0.6), Inches(1.5), width=Inches(6.2))
box(s, 0.6, 5.15, 6.2, 0.3, "Gaver-Stehfest, 18 terms, vectorized over particles and times", size=13, bold=True)
h = eq(s, "slab", 0.7, 5.5, 1.5)
eq(s, "steh", 0.7, 5.5 + h + 0.12, 1.35)
eq(s, "erfc", 4.4, 5.55, 1.35)
rows_v = [("Check", "Reference", "Result"),
          ("beta_G -> inf; beta_G = 0.3", "erfc; piston heat balance", "exact; 1.0002"),
          ("beta_G = 1", "finite-difference PDE solution", "9e-3 -> 5e-3 on refinement"),
          ("single path, infinite rock", "Lauwerier closed form", "1e-16"),
          ("single path, B = 5 m", "slab kernel", "exact"),
          ("one edge, 20000 draws", "TDRW infinite-matrix sampler (Monte Carlo)", "0.008"),
          ("two weighted paths; HDF5", "superposition; round trip", "pass")]
tx, ty = 7.3, 1.55
tbl = s.shapes.add_table(len(rows_v), 3, Inches(tx), Inches(ty), Inches(5.45), Inches(3.4)).table
tbl.columns[0].width, tbl.columns[1].width, tbl.columns[2].width = Inches(1.75), Inches(2.1), Inches(1.6)
for i, row in enumerate(rows_v):
    for j, val in enumerate(row):
        cell = tbl.cell(i, j); cell.text = ""
        r = cell.text_frame.paragraphs[0].add_run(); r.text = val
        r.font.size = Pt(11 if i else 12); r.font.name = BODY
        r.font.bold = (i == 0); r.font.color.rgb = WHITE if i == 0 else INK
        cell.fill.solid(); cell.fill.fore_color.rgb = GARNET if i == 0 else (WHITE if i % 2 else ROCK)
        cell.margin_left = cell.margin_right = Inches(0.06); cell.margin_top = cell.margin_bottom = Inches(0.04)
box(s, tx, 5.3, 5.45, 1.3, [[("The TDRW check matters most: ", {"bold": True}), ("the analytic kernel reproduces what the particle tracker's own matrix-diffusion sampler would do, without sampling. The finite-slab TDRW models remain unusable at thermal diffusivity (10^7 events per particle); the kernel replaces them.", {})]], size=12)
footer(s)

# ---------------------------------------------------------------- 5 dataset status
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """Where the example stands. The driver runs end to end: a 400 by 200 by 200 meter block, hydraulic fractures normal to y spanning the block, two natural sets, flow left to right at a target 10 kilograms per second with the pressure drop solved for, spacings 10 to 80 meters and two seeds. Six of eight samples wrote their HDF5 files; the table shows the fracture counts and the pressure drops, which scale as expected with fewer fractures. But the produced temperature stays at the reservoir value for the full 30 years in every sample, and that is a defect, not a result. The particle exchange integral beta is about fifty times the parallel-plate value: the pipe-network edges carry flow through an area of aperture times roughly 40 meters, not the full 200 meter fracture height, and detours through micron-aperture natural fractures add more, so the thermal time scale comes out at millions of years. The fix is to evaluate beta at the fracture level from the per-fracture flow partition, or to correct the edge widths, and then to re-verify on a one-fracture DFN against Lauwerier. The second defect is connectivity at 80 meter spacing, where dfnGen's largest-cluster removal left three fractures and graph flow had no inlet. Both are recorded in Documentation/thermal_status.md with the numbers."""
title(s, "Example status: pipeline runs, two defects to fix")
box(s, 0.6, 1.5, 6.0, 2.1, [
    [("Design. ", {"bold": True}), ("400 x 200 x 200 m block; hydraulic fractures normal to y spanning x and z, k = 1e-8 m^2; two natural TPL sets 20 to 60 m, p32 0.02 each, k = 1e-12 m^2; flow left to right at a target 10 kg/s (pressure drop solved for); 10^4 particles; 30-year window.", {})],
    [("Sweep. ", {"bold": True}), ("D = 10, 20, 40, 80 m x seeds 1, 2. Six of eight samples written; the 80 m seed-1 sample failed in dfnGen.", {})],
], size=12, spacing=6)
rows_t = [("D [m]", "n_hf", "fractures", "dp [kPa]", "n_exit", "T_out(30 yr) [C]")]
for d, rs in sorted(by_sp.items()):
    rows_t.append((f"{d:g}", rs[0]["n_hydraulic"], "/".join(r["num_fractures"] for r in rs),
                   "/".join(f"{float(r['pressure_drop_Pa'])/1e3:.1f}" for r in rs),
                   "/".join(r["n_exit"] for r in rs),
                   "/".join(f"{float(r['T_out_30yr']):.0f}" for r in rs)))
tbl = s.shapes.add_table(len(rows_t), 6, Inches(0.6), Inches(3.75), Inches(6.0), Inches(1.6)).table
for j, w in enumerate((0.7, 0.7, 1.1, 1.1, 1.2, 1.2)):
    tbl.columns[j].width = Inches(w)
for i, row in enumerate(rows_t):
    for j, val in enumerate(row):
        cell = tbl.cell(i, j); cell.text = ""
        r = cell.text_frame.paragraphs[0].add_run(); r.text = val
        r.font.size = Pt(10); r.font.name = BODY; r.font.bold = (i == 0)
        r.font.color.rgb = WHITE if i == 0 else INK
        cell.fill.solid(); cell.fill.fore_color.rgb = GARNET if i == 0 else (WHITE if i % 2 else ROCK)
        cell.margin_left = cell.margin_right = Inches(0.05); cell.margin_top = cell.margin_bottom = Inches(0.03)
box(s, 0.6, 5.45, 6.0, 0.4, "seed 1 / seed 2 per cell; dp is the pressure drop that delivers 10 kg/s", size=10, italic=True, color=MUTED)
x = 7.0
rect(s, x, 1.5, 5.7, 5.1, TINT, line=GARNET, shape=MSO_SHAPE.ROUNDED_RECTANGLE).adjustments[0] = 0.04
box(s, x + 0.3, 1.65, 5.1, 0.35, "Defect 1: exchange integral 50x too large", size=14, bold=True, color=GARNET)
box(s, x + 0.3, 2.05, 5.1, 1.9, [
    [("Parallel-plate estimate at D = 10 m: beta = 2LH/Q_f = 3.2e8 s/m, slab piston time ~26 yr. Stored particles: median beta 1.6e10 s/m, median beta_G 4e-4, median t_s 3e6 yr, so no cooling in 30 yr.", {})],
    [("Cause: pipe-network edges carry flow through b x ~40 m (median edge area 1.4e-2 m^2), not the 200 m fracture height; micron-aperture natural detours add more.", {})],
    [("Fix: fracture-level beta from the per-fracture flow partition along frac_seq, or corrected edge widths; re-verify on a one-fracture DFN against Lauwerier.", {"bold": True})],
], size=11, spacing=4)
box(s, x + 0.3, 4.15, 5.1, 0.35, "Defect 2: connectivity at wide spacing", size=14, bold=True, color=GARNET)
box(s, x + 0.3, 4.55, 5.1, 1.0, [
    [("At D = 80 m dfnGen's largest-cluster removal kept 3 of 133 fractures and graph flow found no inlet nodes. Guarantee connectivity (natural intensity or size), or catch and skip the sample.", {})],
], size=11, spacing=4)
box(s, x + 0.3, 5.65, 5.1, 0.8, "Not a cause: the Stehfest kernel reproduces the step at t_D = beta_G within 0.1 down to beta_G = 0.02. Details and numbers: Documentation/thermal_status.md.", size=11, italic=True, color=MUTED)
footer(s)

# ---------------------------------------------------------------- 6 limits + next
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """What the model does and does not do. It takes the DFN's flow partition and connectivity at face value and drops the equal-flow assumption of a parallel-plate surrogate. It applies one slab half-width, D over 2, to every edge, natural fractures included, which overstates the rock a natural fracture can draw on only where it runs parallel to the hydraulic planes; B is an argument, so it can be swept. Particles that never exit are dropped and counted in n_stuck. The boundary condition is face to face rather than a well doublet; the well package can replace it once the mesh-free well intersection is wired to the graph. For the ML team the natural next moves are more seeds and spacings, sweeps of natural intensity and hydraulic permeability, and storing the per-fracture flow partition alongside the curves so the features are richer than the headline numbers."""
rect(s, 0, 0, 13.333, 7.5, DARK)
title(s, "Scope, limits, and what to sweep next", dark=True)
cards = [("What the DFN adds", None, "Flow partition among fractures and connectivity through the natural network come from graph flow; no equal-flow assumption. The result is deterministic given the particles."),
         ("One slab width everywhere", "w_B", "is applied to every edge, natural fractures included. B is an argument of heat_extraction and can be swept or set per study."),
         ("Boundaries and losses", None, "Flow is face to face, not a well doublet (the well package can replace it). Particles that never exit are dropped and counted in n_stuck.")]
for i, (head, nm, text) in enumerate(cards):
    x = 0.6 + i * 4.1
    rect(s, x, 1.7, 3.85, 3.6, RGBColor(0x3A, 0x3D, 0x55))
    box(s, x + 0.3, 1.95, 3.25, 1.0, head, size=18, bold=True, font=HEAD, color=RGBColor(0xE8, 0xA8, 0x7C))
    y = 3.0
    if nm:
        y += eq(s, nm, x + 0.3, y, 1.6) + 0.15
    box(s, x + 0.3, y + 0.1, 3.25, 2.0, text, size=13, color=WHITE)
box(s, 0.6, 5.7, 12.1, 0.9, "Next for the dataset: more seeds and spacings, sweeps of natural-set intensity and hydraulic permeability, and the per-fracture flow partition stored with each sample as additional features.", size=14, italic=True, color=RGBColor(0xC8, 0xC9, 0xD4))
footer(s, dark=True)

prs.save("thermal_dfn_dataset_method.pptx")
print("saved")
