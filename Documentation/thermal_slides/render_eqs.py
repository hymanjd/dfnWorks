"""Typeset each equation with LaTeX (same MacTeX engine LaTeXiT drives) to a transparent PNG."""
import subprocess, pathlib, json
EQS = {
 "map":   r"\phi_m \;\to\; \frac{(\rho c)_r}{(\rho c)_f},\qquad D_m \;\to\; \alpha = \frac{k}{(\rho c)_r}",
 "betaP": r"\beta_p = \sum_{i\in\text{path}} \frac{2\,\tau_i}{b_i}",
 "kern":  r"\widetilde{R}_p(s) = \exp\!\Big[-\phi_m\sqrt{D_m s}\;\beta_p\,\tanh\!\big(B\sqrt{s/D_m}\big)\Big]",
 "F":     r"F(t) = \sum_p w_p\, T_{wD}\!\left(\frac{t-\tau_p}{t_{s,p}};\,\beta_{G,p}\right),\qquad t_{s,p} = \big(\phi_m\sqrt{D_m}\,\beta_p\big)^2,\qquad \beta_{G,p} = \frac{B}{\phi_m D_m \beta_p}",
 "Tout":  r"T_{\mathrm{out}}(t) = T_{\mathrm{res}} - (T_{\mathrm{res}}-T_{\mathrm{inj}})\,F(t),\qquad \dot Q(t) = \dot m\, c_f\,(T_{\mathrm{res}}-T_{\mathrm{inj}})\,[1-F(t)]",
 "slab":  r"\bar T_{wD}(u) = \frac{1}{u}\exp\!\big[-\sqrt{u}\,\tanh(\beta_G\sqrt{u})\big]",
 "erfc":  r"T_{wD} \to \operatorname{erfc}\!\left(\frac{1}{2\sqrt{t_D}}\right)\ \ (\beta_G\to\infty)",
 "steh":  r"T_{wD}(t_D) \approx \frac{\ln 2}{t_D}\sum_{i=1}^{18} V_i\, \bar T_{wD}\!\left(\frac{i\ln 2}{t_D}\right)",
 "w_B":   r"B = D/2\ \text{everywhere, natural fractures included}",
}
out = pathlib.Path("eq"); out.mkdir(exist_ok=True)
sizes = {}
for name, body in EQS.items():
    (out / f"{name}.tex").write_text(
        "\\documentclass[preview,border=1.5pt]{standalone}\n\\usepackage{amsmath,amssymb,lmodern}\n"
        "\\begin{document}\n$\\displaystyle " + body + "$\n\\end{document}\n")
    subprocess.run(["latex", "-interaction=nonstopmode", f"{name}.tex"], cwd=out, check=True,
                   stdout=subprocess.DEVNULL)
    fg = "rgb 1 1 1" if name.startswith("w_") else "rgb 0.12 0.13 0.17"
    subprocess.run(["dvipng", "-q", "-D", "600", "-T", "tight", "-bg", "Transparent", "-fg", fg,
                    "-o", f"{name}.png", f"{name}.dvi"], cwd=out, check=True, stdout=subprocess.DEVNULL)
    from PIL import Image
    sizes[name] = Image.open(out / f"{name}.png").size
json.dump(sizes, open(out / "sizes.json", "w"))
print({k: (round(w/600, 2), round(h/600, 2)) for k, (w, h) in sizes.items()})
