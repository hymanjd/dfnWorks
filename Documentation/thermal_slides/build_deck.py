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


def title(slide, text, dark=False):
    box(slide, 0.6, 0.5, 12.1, 0.8, text, size=30, bold=True, font=HEAD,
        color=WHITE if dark else DARK)


def footer(slide, dark=False):
    box(slide, 0.6, 7.0, 8, 0.3, "pydfnworks.thermal  ·  dfnWorks, Los Alamos National Laboratory",
        size=10, color=RGBColor(0xC8, 0xC9, 0xD4) if dark else MUTED)


# ---------------------------------------------------------------- 1 title
s = prs.slides.add_slide(BLANK)
rect(s, 0, 0, 13.333, 7.5, DARK)
rect(s, 0.6, 2.2, 0.12, 2.1, EMBER)
box(s, 0.95, 2.1, 11, 1.6, "Fracture spacing optimization with an analytic heat-extraction surrogate",
    size=40, bold=True, font=HEAD, color=WHITE)
box(s, 0.95, 3.75, 11, 0.6, "The Gringarten multiple-parallel-fracture model as a 100-microsecond objective inside a dfnWorks workflow",
    size=18, color=RGBColor(0xE8, 0xA8, 0x7C), italic=True)
s.notes_slide.notes_text_frame.text = "A new pydfnworks subpackage for heat extraction from fractured rock, built to sit inside a hydraulic-fracture spacing optimizer. Six slides: the interface, the model, the equations, the numerics and verification, the workflow on a FORGE-type network, and what the DFN-based backend has to add."
box(s, 0.95, 5.6, 11, 0.8, ["pydfnworks.thermal  ·  branch feature/thermal-mpf-surrogate",
                            "Jeffrey Hyman, EES-16, Los Alamos National Laboratory  ·  November 2026"],
    size=14, color=RGBColor(0xC8, 0xC9, 0xD4))

# ---------------------------------------------------------------- 2 problem + model sketch
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """The question is how to space hydraulic fractures for heat extraction at a FORGE-type site. The answer has to come from an optimizer, and an optimizer needs an objective it can call thousands of times. The design separates the objective from the physics behind it: one function signature, heat_extraction, with two backends. The first is the Gringarten multiple-parallel-fracture model, analytic in the Laplace domain and inverted numerically, which costs about 100 microseconds per call and is complete and verified. The second is the dfnWorks graph-transport backend, which will use the DFN flow partition and connectivity; its interface and parameter mapping are in place and the kernel is deferred. The sketch on the right is the surrogate's world: n identical fractures a distance D apart, water flowing up each one, heat conducting in from the rock slab on either side. Equal flow in every fracture is the assumption the DFN backend will eventually relax."""
title(s, "One interface, two backends")
box(s, 0.6, 1.55, 5.9, 1.4, [
    [("The question. ", {"bold": True}), ("How far apart should hydraulic fractures be spaced to maximize heat extracted from a FORGE-type enhanced geothermal reservoir?", {})],
    [("The design. ", {"bold": True}), ("An optimizer calls one function, and two backends implement it:", {})],
], size=15, spacing=10)
eq(s, "sig", 0.8, 2.95, 1.45)
box(s, 0.6, 3.45, 5.9, 3.3, [
    [("Phase 1  ", {"bold": True, "color": GARNET}), ("Gringarten et al. (1975) multiple parallel fractures, inverted from the Laplace domain. Done, verified, ~100 μs per call.", {})],
    [("Phase 2  ", {"bold": True, "color": GARNET}), ("dfnWorks graph transport with thermal TDRW parameters. Interface and parameter mapping in place; kernel deferred.", {})],
    [("Why a surrogate first. ", {"bold": True}), ("The optimizer needs thousands of evaluations; a particle run costs seconds to minutes. The surrogate sets the interface and the baseline the DFN backend must beat.", {})],
], size=15, spacing=10)

# sketch: n parallel fractures with slabs
x0, y0 = 7.2, 1.7
box(s, x0, y0 - 0.1, 5.5, 0.4, "n parallel fractures, spacing D, in an infinite rock mass", size=12, color=MUTED, italic=True)
rect(s, x0, y0 + 0.45, 5.5, 3.6, ROCK)
nf = 4
pitch = 5.5 / nf
for i in range(nf):
    fx = x0 + pitch * (i + 0.5)
    rect(s, fx - 0.04, y0 + 0.45, 0.08, 3.6, GARNET)
    # flow arrow up the fracture
    a = s.shapes.add_shape(MSO_SHAPE.UP_ARROW, Inches(fx - 0.17), Inches(y0 + 1.5), Inches(0.34), Inches(1.4))
    a.fill.solid(); a.fill.fore_color.rgb = EMBER; a.line.fill.background(); a.shadow.inherit = False
    if i < nf - 1:
        # conduction arrows from the slab into the fracture
        for yy in (y0 + 1.0, y0 + 3.3):
            c = s.shapes.add_shape(MSO_SHAPE.LEFT_RIGHT_ARROW, Inches(fx + 0.15), Inches(yy), Inches(pitch - 0.3), Inches(0.22))
            c.fill.solid(); c.fill.fore_color.rgb = MUTED; c.line.fill.background(); c.shadow.inherit = False
box(s, x0, y0 + 4.1, 1.6, 0.35, "T_inj in", size=12, bold=True, color=EMBER)
box(s, x0 + 3.9, y0 + 0.05, 1.6, 0.35, "T_out", size=12, bold=True, color=EMBER, align=PP_ALIGN.RIGHT)
box(s, x0 + pitch * 1.5 - 0.6, y0 + 4.1, 1.2, 0.35, "← D →", size=12, color=INK, align=PP_ALIGN.CENTER)
box(s, x0 + 5.55, y0 + 2.0, 0.6, 0.5, "H", size=12, color=INK)
box(s, x0, y0 + 4.55, 5.5, 0.9, "Each fracture draws on a slab of half-thickness D/2 by conduction normal to its plane; the slab mid-plane is insulated by symmetry. Equal flow in every fracture.",
    size=11, color=MUTED)
footer(s)

# ---------------------------------------------------------------- 3 equations
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """Everything in the surrogate follows from Equation A17 of Gringarten, Witherspoon and Ohnishi. In the Laplace domain the outlet drawdown is an exponential of root s times a hyperbolic tangent; the tanh is the finite slab. Two groups carry all the physics. The spacing parameter beta is the slab half-width in units of the conduction length set by the flow rate, so it is proportional to spacing times flow per fracture over conductivity and fracture area. The time scale t_s makes the dimensionless time. The heat rate is just mass flow times heat capacity times the available temperature drop, times one minus the drawdown. The two limits on the right are what the tests check. When beta is large the tanh goes to one and the inverse is the classic erfc solution for a single fracture between semi-infinite rock. When beta is small the response is a step at t_D equals beta, which in dimensional terms is the piston time: the heat in the slabs divided by the heat-capacity flow of the water. That piston time is the quantity that makes the optimization answer monotone, which comes up on the last slide."""
title(s, "The surrogate: Laplace solution and dimensionless groups")
rect(s, 0.6, 1.5, 12.1, 1.25, ROCK)
box(s, 0.9, 1.58, 11.5, 0.3, "Outlet drawdown, Laplace domain (Gringarten et al. 1975, Eq. A17)", size=13, color=MUTED)
eq(s, "a17", 0.9, 1.95, 1.8)
y = 3.05
for label, name in (("Spacing parameter", "beta"), ("Time scale", "ts"), ("Heat rate", "Q")):
    box(s, 0.6, y, 6.0, 0.3, label, size=14, bold=True, color=GARNET)
    y += 0.34 + eq(s, name, 0.7, y + 0.32, 1.7) + 0.12
h = eq(s, "Qf", 0.7, y + 0.08, 1.4)
box(s, 0.6, y + 0.1 + h + 0.05, 6.0, 0.6, "Q_f is the flow per fracture; H, w the fracture height (flow direction) and width; k, (ρc)_r rock conductivity and heat capacity.", size=12, color=MUTED)
x = 7.0
box(s, x, 3.05, 5.7, 0.3, "Two limits anchor the model", size=14, bold=True, color=GARNET)
box(s, x, 3.5, 5.7, 0.3, [[("β → ∞, wide spacing: ", {"bold": True}), ("semi-infinite single fracture (Lauwerier 1955; Bödvarsson and Tsang 1982)", {})]], size=13)
h = eq(s, "erfc", x + 0.1, 3.95, 1.9)
box(s, x, 3.95 + h + 0.3, 5.7, 0.3, [[("β → 0, close spacing: ", {"bold": True}), ("a step at t_D = β, the piston displacement of the slab", {})]], size=13)
h2 = eq(s, "tp", x + 0.1, 3.95 + h + 0.7, 1.9)
box(s, x, 3.95 + h + 0.7 + h2 + 0.2, 5.7, 0.7, "the heat stored in the n slabs divided by the heat-capacity flow of the water. Both limits are closed-form checks on the numerical inversion.", size=12, color=MUTED)
footer(s)

# ---------------------------------------------------------------- 4 inversion + verification
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """The inversion is Gaver-Stehfest with 18 terms, the same routine the dfnWorks TDRW slab models already use, vectorized so that a full 30-year curve costs about 100 microseconds. GEOPHIRES-X inverts the same expression with mpmath at 15 digits and takes about a second, which is fine for a single reservoir run and useless inside an optimizer. The term count was chosen against a 40-digit de Hoog reference; 18 is uniformly at least as good as 16, and 20 or more loses accuracy in double precision. The curves show the family: close spacing, small beta, holds the full rate and then collapses; wide spacing follows the erfc decline from the start. The table is the verification. The erfc limit is matched to six parts in a million. The energy balance at close spacing closes to a tenth of a percent. The mpmath route agrees to one part in ten thousand. And because there is no closed form at finite beta, there is a finite-difference solution of the slab conduction and fracture advection equations that never touches the Laplace transform; it converges toward the inverted solution on grid refinement. Together those cover the whole beta range."""
title(s, "Inverted in double precision, verified four independent ways")
s.shapes.add_picture("drawdown_curves.png", Inches(0.6), Inches(1.6), width=Inches(6.4))
box(s, 0.6, 5.3, 6.4, 0.3, "Gaver-Stehfest inversion, 18 terms, vectorized over time", size=14, bold=True)
eq(s, "steh", 0.7, 5.65, 1.7)
box(s, 4.6, 5.65, 2.4, 1.0, "Shares the routine of the dfnWorks TDRW slab models. 361 points in ~100 μs; the GEOPHIRES-X mpmath route takes ~1 s.", size=11, color=MUTED)
rows = [("Check", "Reference", "Result"),
        ("Wide spacing, β = 3·10⁵", "erfc closed form", "max error 6·10⁻⁶"),
        ("Close spacing, β = 0.6 and 1.4", "energy balance, heat in slabs", "1.0001 and 0.9988"),
        ("β = 0.5, 1.5, 6", "GEOPHIRES-X mpmath Stehfest", "≤ 1·10⁻⁴"),
        ("Finite spacing, β = 1", "finite-difference solution of the PDEs", "9·10⁻³ → 5·10⁻³ on refinement"),
        ("Close vs wide, same flow", "physics", "never more heat; 75 % drawdown at 51 vs 123 d")]
tx, ty = 7.3, 1.65
tbl = s.shapes.add_table(len(rows), 3, Inches(tx), Inches(ty), Inches(5.45), Inches(3.3)).table
tbl.columns[0].width, tbl.columns[1].width, tbl.columns[2].width = Inches(1.75), Inches(1.9), Inches(1.8)
for i, row in enumerate(rows):
    for j, val in enumerate(row):
        cell = tbl.cell(i, j)
        cell.text = ""
        p = cell.text_frame.paragraphs[0]
        r = p.add_run(); r.text = val
        r.font.size = Pt(11 if i else 12); r.font.name = BODY
        r.font.bold = (i == 0); r.font.color.rgb = WHITE if i == 0 else INK
        cell.fill.solid(); cell.fill.fore_color.rgb = GARNET if i == 0 else (WHITE if i % 2 else ROCK)
        cell.margin_left = cell.margin_right = Inches(0.06)
        cell.margin_top = cell.margin_bottom = Inches(0.04)
box(s, tx, 5.2, 5.45, 1.4, [
    [("No closed form exists at finite β", {"bold": True}), (", so the general case is checked against a finite-difference solution of the slab conduction and fracture advection equations that never sees the Laplace transform. The erfc limit, the energy balance and the finite-difference run together cover the whole β range.", {})],
], size=13)
footer(s)

# ---------------------------------------------------------------- 5 workflow + result
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """The workflow is the point of putting this inside pydfnworks rather than in a standalone script. dfnWorks generates a FORGE-representative network, the well package places an injection lateral and a production lateral and reports which fractures each cuts, and from that the surrogate gets its geometry: the fracture count, the mean diameter as the fracture height, and the lateral length. The optimizer then evaluates candidate spacings L over n for integer stage counts. The figure is the run as committed. On the left, 30-year cumulative heat against spacing for stages filling the lateral: it is flat near the closest spacing and falls off above about 20 meters, so the best is the closest allowed, 7.7 meters and 78 stages, 11.6 petajoules. On the right, the heat rate at that spacing holds 15.6 megawatts for 18 years and then drops as the slabs run out, while the spacing the DFN itself implies, 150 meters and five fractures, declines from the first day and extracts about half as much. Swapping the DFN backend in is one argument to the optimizer."""
title(s, "Workflow: DFN geometry in, surrogate objective out")
steps = [("1", "Generate", "FORGE-representative DFN: 600 m granitoid block, a set parallel to the hydraulic fractures and a natural set"),
         ("2", "Wells", "Injection and production laterals placed with the well package; fractures each cuts are listed"),
         ("3", "Geometry", "n = fractures cut by the injector, H = their mean diameter, L = lateral length. Here n = 5, H = 219 m, L = 600 m"),
         ("4", "Optimize", "optimize_spacing() over candidates L/n, objective = 30-year cumulative heat, 20 kg/s, 226 °C reservoir, 40 °C injection")]
for i, (num, head, text) in enumerate(steps):
    y = 1.6 + i * 1.25
    c = rect(s, 0.6, y, 0.55, 0.55, GARNET, shape=MSO_SHAPE.OVAL)
    box(s, 0.6, y + 0.1, 0.55, 0.4, num, size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    box(s, 1.35, y - 0.02, 4.2, 0.4, head, size=16, bold=True, color=DARK)
    box(s, 1.35, y + 0.36, 4.2, 0.9, text, size=12, color=INK)
s.shapes.add_picture("forge_result.png", Inches(5.9), Inches(1.6), width=Inches(6.9))
box(s, 5.9, 4.45, 6.9, 2.3, [
    [("Read the result. ", {"bold": True}), ("Left: 30-year heat vs spacing for stages filling the lateral, best at the closest allowed spacing (7.7 m, 78 stages, 11.6 PJ). Right: the heat rate at that spacing holds 15.6 MW for 18 years and then collapses as the slabs deplete, while the DFN-implied 150 m spacing (5 fractures, 6.4 PJ) declines from day one.", {})],
    [("Same optimizer, other backend. ", {"bold": True}), ("Pointing heat_extraction at the DFN backend changes one argument; it raises NotImplementedError until the kernel lands.", {})],
], size=13, spacing=6)
footer(s)

# ---------------------------------------------------------------- 6 findings
s = prs.slides.add_slide(BLANK)
s.notes_slide.notes_text_frame.text = """Three things to take away. First, with the surrogate alone the objective is monotone in spacing, and that is physics, not a weakness of the optimizer: with stages filling a fixed lateral at fixed flow, n times D is the lateral length, so the piston time does not depend on spacing at all and the closest spacing always approximates the ideal piston front best. At fixed fracture count the widest spacing wins because it adds rock. An interior optimum needs something the surrogate does not contain, either a cost per stage or the DFN flow partition. Second, heat is a solute with two coefficient substitutions: matrix porosity becomes the rock-to-fluid heat capacity ratio, about one half, and matrix diffusivity becomes the rock thermal diffusivity. With those, run_graph_transport already runs a thermal problem; the infinite-matrix model handled ten thousand particles in two seconds. Third, the finite-slab sampler cannot be reused as is. It samples excursions event by event at a rate proportional to porosity times diffusivity, which is ten million times larger for heat, so a single particle took minutes. The Phase 2 kernel has to draw the whole slab retention time per edge. That, plus the arrival-time to heat-rate conversion, is the remaining work."""
rect(s, 0, 0, 13.333, 7.5, DARK)
title(s, "What the surrogate says, and what Phase 2 has to add", dark=True)
cards = [
    ("Objectives are monotone in spacing", "w_tpL",
     "At fixed flow and fixed rock volume the piston time does not depend on D, so the closest spacing always wins; at fixed fracture count the widest does. An interior optimum needs a cost per stage or the DFN flow partition."),
    ("Heat is a solute with new coefficients", ("w_phi", "w_D"),
     "turn run_graph_transport() into a thermal run. The infinite-matrix model handles 10⁴ particles in 2 s."),
    ("The slab sampler cannot be reused as is", "w_gam",
     "is the event-by-event trapping rate: 10⁷ times faster for heat than for a solute, minutes per particle. Phase 2 draws the whole slab retention time per edge from the finite-slab kernel instead."),
]
for i, (head, names, text) in enumerate(cards):
    x = 0.6 + i * 4.1
    rect(s, x, 1.7, 3.85, 3.9, RGBColor(0x3A, 0x3D, 0x55))
    box(s, x + 0.3, 1.95, 3.25, 1.0, head, size=18, bold=True, font=HEAD, color=RGBColor(0xE8, 0xA8, 0x7C))
    y = 3.0
    for nm in (names if isinstance(names, tuple) else (names,)):
        y += eq(s, nm, x + 0.3, y, 1.65) + 0.15
    box(s, x + 0.3, y + 0.1, 3.25, 1.8, text, size=13, color=WHITE)
box(s, 0.6, 5.95, 12.1, 0.5, "Next: per-edge finite-slab kernel and the arrival-time → Q̇(t) conversion (Shook 2001), then the same optimizer with the DFN backend.",
    size=14, italic=True, color=RGBColor(0xC8, 0xC9, 0xD4))
footer(s, dark=True)

prs.save("thermal_surrogate_method.pptx")
print("saved")
