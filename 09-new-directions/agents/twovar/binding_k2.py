"""Constraint values of the k=2 cubic probe at its GRH-floor vertex (MODEL)."""
import numpy as np
c = 1/3
def energy(M, lp=0.0):
    f, l = M, 1 + lp; lR = 2 * (M + lp) - (1 + lp)
    nodes = {"direct cubic LS": max(0, l - f, 2/3*(l+f) - f), "Poisson": max(0, l - 2/3*f),
             "reflect (cB=2)+cubic LS": max(0, lR - f, 2/3*(lR+f) - f), "reflect+Poisson": max(0, lR - 2/3*f)}
    third = (M + 4*lp - 2) / 3
    return nodes, third
for (l, y) in [((0.0, 0.0), (0.75, 0.75)), ((0.679, 0.071), (0.679, 0.071))]:
    l = np.array(l); y = np.array(y); M = l.sum() + y.sum()
    base = 0.5 + sum((0.5 - c) + c * li for li in l)
    gram = sum(max(0, c * (yi - li)) for li, yi in zip(l, y))
    nodes, third = energy(M)
    e = max(min(nodes.values()), third, 0)
    print(f"geometry l={l}, y={y}, M={M:.3f}")
    print(f"  low base (signal/low) = {base:.4f};  + Gram exceptional/2 = {base + gram/2:.4f};  energy excess e = {e:.4f} "
          f"(nodes { {k: round(v,3) for k,v in nodes.items()} }, third {third:.3f})")
    print(f"  LOW total = {base + gram/2 + e/2:.4f}")
    gen = 0.5 + sum(0.5*(1 - li - yi) + (0.5 - c)*(1 - li) for li, yi in zip(l, y))
    p1 = 0.5 + 0.5*(1 - l[1] - y[1]) + (0.5 - c)*(1 - l[1])
    p2 = 0.5 + 0.5*(1 - l[0] - y[0]) + (0.5 - c)*(1 - l[0])
    print(f"  GRH rows: generic {gen:.4f}; copy in var1 only {p1:.4f}; copy in var2 only {p2:.4f}")
