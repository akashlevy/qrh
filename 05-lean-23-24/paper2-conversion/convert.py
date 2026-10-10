#!/usr/bin/env python3
"""Convert paper2.tex (The Quasi-Riemann Hypothesis, 11/12) into the 23/24 paper.

Every change is an explicit replacement that must match exactly once, so the result is the
original text with a reviewable list of edits:
  * only zeta (the Hecke character nu is trivial),
  * the fixed parameter vartheta = 1/10, giving 23/24,
  * extraction from the rows u = b^6 (no prime ideal theorem), as in the Lean proof,
  * the top level (Section 3) split into the three propositions of the Lean spine,
  * parts the argument does not need are cut (consequences, prior work, the proofs of
    Appendices A and B, which are cited from the original),
  * remarks and an appendix recording the correspondence with the Lean formalization.

Usage: python3 paper23/convert.py SRC_PAPER2_TEX OUT_TEX
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rep(text: str, old: str, new: str, count: int = 1) -> str:
    n = text.count(old)
    if n != count:
        sys.exit(f"expected {count} occurrence(s), found {n}:\n{old[:300]}")
    return text.replace(old, new)


def between(text: str, start: str, end: str) -> tuple[int, int]:
    i = text.index(start)
    j = text.index(end, i)
    return i, j


def splice(text: str, start: str, end: str, new: str) -> str:
    """Replace text[start:end) (start marker included, end marker kept) by new."""
    i, j = between(text, start, end)
    return text[:i] + new + text[j:]


def piece(name: str) -> str:
    return (HERE / "pieces" / name).read_text()


def main() -> int:
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    t = src.read_text()

    # ---------------------------------------------------------------- front matter
    t = rep(t, r"\hypersetup{pdftitle={The Quasi-Riemann Hypothesis},pdfauthor={OpenAI},"
               r"pdfsubject={A cubic theta argument with a fixed ray class character}}",
            r"\hypersetup{pdftitle={The Quasi-Riemann Hypothesis: the half-plane Re s > 23/24 for Dirichlet L-functions},"
            r"pdfsubject={A cubic theta argument, matched to a Lean formalization}}")
    t = rep(t, "\\title{The Quasi-Riemann Hypothesis}\n\\author{OpenAI}\n\\date{October 5, 2026}",
            "\\title{The Quasi-Riemann Hypothesis:\\\\ a zero-free half-plane $\\Rea s>23/24$ for Dirichlet $L$-functions}\n"
            "\\author{Adapted from OpenAI, \\emph{The Quasi-Riemann Hypothesis} (October 5, 2026)}\n"
            "\\date{October 9, 2026}")
    t = splice(t, "We establish the quasi-Riemann hypothesis by proving", "\\end{abstract}", piece("abstract.tex"))

    # ---------------------------------------------------------------- Section 1
    t = rep(t, "Every finite-order Hecke $L$-function over $K$ has no zeros in the half-plane $\\Rea s>11/12$. "
               "Thus every Dirichlet $L$-function, including $\\zeta(s)$, has no zeros for $\\Rea s>11/12$.",
            "Every Dirichlet $L$-function, including $\\zeta(s)$, has no zeros in the half-plane $\\Rea s>23/24$.")
    t = rep(t, "Theorem~\\ref{thm:main} is a natural intermediate step toward the stronger zero-free region "
               "$\\Rea s>7/8$ established in \\cite[Theorem~1.1]{OpenAI26}. For simplicity, we have isolated "
               "Theorem~\\ref{thm:main} and its proof here.\n",
            piece("relation.tex"))
    t = rep(t, " \\ll x^{11/12}\\log x,\n", " \\ll x^{23/24}\\log x,\n")
    t = rep(t, "\\subsection{Organization}", piece("formal.tex") + "\\subsection{Organization}")
    t = rep(t, "Appendix~\\ref{sec:arithmetic} supplies the arithmetic identities and\n"
               "the detailed theta calculation. Appendix~\\ref{app:weights} records the smooth\n"
               "separation lemma used throughout.",
            "Appendix~\\ref{sec:arithmetic} supplies the arithmetic identities and\n"
            "the detailed theta calculation. Appendix~\\ref{app:weights} records the smooth\n"
            "separation lemma used throughout, and Appendix~\\ref{app:lean} the correspondence with\n"
            "the formal proof.")

    # ---------------------------------------------------------------- Section 2 (outline)
    t = rep(t, "Fix a finite-order Hecke character $\\nu$ of $K$. Let $L_K(s, \\nu)$ be the corresponding "
               "Hecke $L$-function.",
            "Fix a Dirichlet character $\\chi$, and let $\\nu=\\chi\\circ\\NK$ be the corresponding finite-order "
            "Hecke character of $K$. Its Hecke $L$-function is $L_K(s,\\nu)=L(s,\\chi)L(s,\\chi\\chi_{-3})$, where "
            "$\\chi_{-3}$ is the nontrivial character modulo $3$.")
    t = rep(t, "Section~\\ref{sec:reduction} gives the complete argument needed here.\n",
            "Section~\\ref{sec:reduction} gives the complete argument needed here, for a variant of\n"
            "$A_1(D)$ with an additional weight $\\varphi(\\mathfrak n)/\\NK(\\mathfrak n)$ (see Step~2).\n")
    t = rep(t, "Fix $0<\\vartheta\\le1/10$. The crucial estimate, Proposition~\\ref{thm:ms}, is that",
            "Fix $\\vartheta=1/10$. The crucial estimate, Proposition~\\ref{thm:ms}, is that")
    t = splice(t, "The mean-square estimate \\eqref{eq:intro-ms} gives the desired power saving",
               "\\subsection{Step 3: Poisson summation}", piece("outline-extraction.tex"))

    # ---------------------------------------------------------------- Section 3
    t = splice(t, "\\section{From the mean-square estimate to the zero-free region}\\label{sec:reduction}",
               "\\section{Poisson summation and the dual mean square}", piece("reduction.tex"))

    # ---------------------------------------------------------------- Sections 4-7: vartheta = 1/10
    i = t.index("\\section{Poisson summation and the dual mean square}")
    j = t.index("\\appendix")
    body = t[i:j]
    body = rep(body, "Fix $0<\\vartheta\\le1/10$ and put $H=D^{1+\\vartheta}$.",
               "Let $\\vartheta=1/10$ and put $H=D^{1+\\vartheta}$.", count=2)
    # remarks on the formal bookkeeping
    body = rep(body, "\\section[Iteration of the dual mean-square estimate]",
               piece("remark-dual.tex") + "\n\\section[Iteration of the dual mean-square estimate]")
    body = rep(body, "\\subsection{The original mean square and the exponent $11/12$}",
               piece("remark-descent.tex") + "\n\\subsection{The original mean square}")
    t = t[:i] + body + t[j:]

    # ---------------------------------------------------------------- appendices
    # Appendices A and B are kept unchanged (they never use nu); Appendix C is new.
    t = rep(t, "\\begin{thebibliography}", piece("lean.tex") + "\n\\begin{thebibliography}")

    # ---------------------------------------------------------------- bibliography
    t = rep(t, "\\bibitem{OpenAI26}\n", piece("bib-openai26b.tex") + "\\bibitem{OpenAI26}\n")
    t = prune_bibliography(t)

    out.write_text(t)
    return 0


def prune_bibliography(t: str) -> str:
    i = t.index("\\begin{thebibliography}")
    j = t.index("\\end{thebibliography}")
    body, bib = t[:i], t[i:j]
    cites = set()
    for m in re.finditer(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", body):
        cites.update(k.strip() for k in m.group(1).split(","))
    head, *items = re.split(r"(?=\\bibitem\{)", bib)
    keys = [re.match(r"\\bibitem\{([^}]*)\}", it).group(1) for it in items]
    missing = cites - set(keys)
    if missing:
        sys.exit(f"cited but not in the bibliography: {sorted(missing)}")
    kept = [it for k, it in zip(keys, items) if k in cites]
    print(f"bibliography: {len(kept)} of {len(items)} entries cited")
    return body + head + "".join(kept) + t[j:]


if __name__ == "__main__":
    sys.exit(main())
