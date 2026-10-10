#!/usr/bin/env python3
"""Check that paper numbers cited in the Lean spine docstrings match the compiled paper.

Every docstring mention of a LaTeX label in backticks (e.g. `prop:poisson-reduction`) must be
preceded, within the same docstring line or the line before, by the number LaTeX assigned to that
label (read from paper/qrh.aux). Run after `tectonic -X compile qrh.tex --keep-intermediates`.
"""
import re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
aux = (ROOT.parent / "paper" / "qrh.aux").read_text()  # repo layout: ../paper (was ROOT/paper)
num = dict(re.findall(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", aux))
bad = 0
for f in sorted((ROOT / "QRH" / "Spine").glob("*.lean")) + [ROOT / "QRH.lean"]:
    lines = f.read_text().split("\n")
    for i, l in enumerate(lines):
        for m in re.finditer(r"`([a-z]+:[A-Za-z0-9-]+)`", l):
            label = m.group(1)
            if label not in num:
                print(f"{f.name}:{i+1}: unknown label {label}"); bad += 1; continue
            context = (lines[i - 1] if i else "") + " " + l[:m.start()]
            cited = re.findall(r"(?:Proposition|Lemma|Theorem|Remark|[Ee]q\.)\s*\(?([0-9A-Z]+\.[0-9]+)\)?", context)
            cited += re.findall(r"\(([0-9A-Z]+\.[0-9]+)\)", context)  # "paper (3.1)" style
            if num[label] not in cited:
                print(f"{f.name}:{i+1}: {label} is {num[label]} in the paper, docstring cites {cited[-1] if cited else 'nothing'}")
                bad += 1
print("numbering OK" if not bad else f"{bad} numbering problems")
sys.exit(1 if bad else 0)
