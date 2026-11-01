"""Typeset each equation with LaTeX (same MacTeX engine LaTeXiT drives) to a transparent PNG."""
import subprocess, pathlib, json
EQS = {
 "a17":   r"\bar T_{wD}(s) = \frac{1}{s}\exp\!\left[-\sqrt{s}\,\tanh\!\big(\beta\sqrt{s}\big)\right],\qquad T_{wD} = \frac{T_{\mathrm{res}}-T_{\mathrm{out}}}{T_{\mathrm{res}}-T_{\mathrm{inj}}},\qquad t_D = \frac{t}{t_s}",
 "beta":  r"\beta = \frac{(\rho c)_f\, Q_f\, D}{4\,k\,w\,H}",
 "ts":    r"t_s = \frac{4\,k\,(\rho c)_r\,(wH)^2}{\big[(\rho c)_f\,Q_f\big]^2}",
 "Q":     r"\dot Q(t) = \dot m\, c_f\,\big(T_{\mathrm{res}}-T_{\mathrm{inj}}\big)\big[1 - T_{wD}(t/t_s;\beta)\big]",
 "Qf":    r"Q_f = \dot m/(\rho_f\, n)",
 "erfc":  r"T_{wD} = \operatorname{erfc}\!\left(\frac{1}{2\sqrt{t_D}}\right)",
 "tp":    r"t_{\mathrm{piston}} = \beta\, t_s = \frac{(\rho c)_r\, n\, D\, w\, H}{(\rho c)_f\, \dot m/\rho_f}",
 "steh":  r"T_{wD}(t_D) \approx \frac{\ln 2}{t_D}\sum_{i=1}^{18} V_i\, \bar T_{wD}\!\left(\frac{i\ln 2}{t_D}\right)",
 "sig":   r"\texttt{heat\_extraction}(D,\ \dot m,\ T_{\mathrm{res}},\ T_{\mathrm{inj}},\ \text{rock},\ \text{fluid},\ n,\ H,\ t)\ \to\ \dot Q(t)",
 # white versions for the dark closing slide
 "w_tpL": r"t_p = \frac{(\rho c)_r\, L\, w\, H}{(\rho c)_f\, \dot m/\rho_f}",
 "w_phi": r"\phi_m \;\to\; \frac{(\rho c)_r}{(\rho c)_f} = 0.52",
 "w_D":   r"D_m \;\to\; \frac{k}{(\rho c)_r} = 1.4\times 10^{-6}\ \mathrm{m^2/s}",
 "w_gam": r"\gamma = \frac{2\,\phi_m D_m}{b\,\varepsilon\,B}",
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
