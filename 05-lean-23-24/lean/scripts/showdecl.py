#!/usr/bin/env python3
"""Print the source of declarations by (short or qualified) name from the QRH library.
Usage: python3 scripts/showdecl.py NAME [NAME...] [--full]   (default: signature up to ':=')"""
import re, sys
from pathlib import Path
full = "--full" in sys.argv
names = [a for a in sys.argv[1:] if a != "--full"]
KW = re.compile(r"^\s*(?:@\[[^\]]*\]\s*)*(?:(?:private|protected|noncomputable|nonrec)\s+)*(theorem|lemma|def|abbrev|structure|class|instance|inductive)\s+([^\s(:{\[]+)")
for p in sorted(Path("QRH").rglob("*.lean")):
    L = p.read_text().split("\n")
    for i, l in enumerate(L):
        m = KW.match(l)
        if not m: continue
        nm = m.group(2)
        if any(nm == n or nm.endswith("." + n) or n.endswith("." + nm) and n.split(".")[-1] == nm.split(".")[-1] for n in names):
            out = [l]; j = i + 1
            while j < len(L) and (L[j].startswith((" ", "\t")) or L[j] == ""):
                if not full and ":=" in out[-1]: break
                out.append(L[j]); j += 1
            if not full:
                k = next((t for t, x in enumerate(out) if ":=" in x), len(out) - 1)
                out = out[:k + 1]
            print(f"--- {p}:{i+1}"); print("\n".join(out).rstrip())
