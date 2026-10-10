#!/usr/bin/env python3
"""Assemble paper/qrh.tex from paper/src/*.tex, keeping only cited bibliography entries."""
import re
from pathlib import Path
src = Path(__file__).parent / "src"
parts = ["00_preamble", "01_intro", "02_outline", "03_reduction", "04_poisson", "05_descent",
         "06_reflection", "07_transfer", "08_appA", "09_appB", "11_lean"]
body = {p: (src / f"{p}.tex").read_text() for p in parts}
front = r"""\begin{document}
\maketitle
\begin{abstract}
We give an account of a proof that the Riemann zeta function has no zeros in the
half-plane $\Rea s>23/24$. The argument is that of OpenAI's \emph{The Quasi-Riemann Hypothesis},
which proves the half-plane $\Rea s>11/12$ for all Dirichlet $L$-functions and all
finite-order Hecke $L$-functions over $\Q(\sqrt{-3})$, specialized here to $\zeta$ and to a fixed value of its parameter. A
M\"obius sum over the ideals of $\mathbb Z[\omega]$ is embedded in a sextic-twisted family. Poisson
summation converts the family's mean square into one with cubic Gauss-sum coefficients;
Kubota's cubic theta function converts the sextic twists into quadratic ones; the quadratic
large sieve then bounds the completed sums, and a finite descent removes the cube factors.
The argument is matched to a Lean~4 formalization: every result of Section~3 is a Lean theorem
with the same statement, and the main results of Sections~4--7 are theorems of the formalization,
proved with somewhat different bookkeeping.
\end{abstract}

\tableofcontents

"""
appendix_start = "\\appendix\n"
text = body["00_preamble"] + front
for p in parts[1:]:
    t = body[p]
    if p == "08_appA":
        t = t.replace("\\appendix\n", "")
        text += appendix_start
    text += t
cites = set()
for m in re.finditer(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", text):
    cites.update(k.strip() for k in m.group(1).split(","))
bib = (src / "10_bib.tex").read_text()
entries = re.split(r"(?=\\bibitem\{)", bib)
head, items = entries[0], entries[1:]
kept = []
for it in items:
    key = re.match(r"\\bibitem\{([^}]*)\}", it).group(1)
    if key in cites:
        kept.append(it.replace("\\end{thebibliography}\n\\end{document}\n", "").rstrip() + "\n\n")
missing = cites - {re.match(r"\\bibitem\{([^}]*)\}", it).group(1) for it in items}
text += head + "".join(kept) + "\\end{thebibliography}\n\\end{document}\n"
(Path(__file__).parent / "qrh.tex").write_text(text)
print(f"cited {len(cites)}, kept {len(kept)} entries; missing: {sorted(missing)}")
