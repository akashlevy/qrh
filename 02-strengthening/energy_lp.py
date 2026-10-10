"""
Max of Lemma reflected-energy's E_ref - M' over all admissible block data,
for the all-marked completed row of Lemma probe-row-norm at rescaled level d.
Variables: O, DH, za, v, lb, el, A0, N0, B0, S0  (H = M'-O-DH, N* = 1+ell')
E_ref = O/2 + max(H, v+lb [-w]) - S0 - B0 + za - u - lb - 2el/3 - (Td - y)/2
u = min(v, za, (v+za)/3), Td = 2H+2A0+2za-N*-N0-3B0, y=v+3lb+el <= Td
Branch enumeration -> LPs.
"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog

def emax(Mp, lp, w=0.0, mode="soft", fix_O0=False):
    # var order: O, DH, za, v, lb, el, A0, N0, B0, S0
    n = 10
    Ns = 1 + lp
    best = (-9, None)
    # H = Mp - O - DH ; Td = 2H + 2A0 + 2za - Ns - N0 - 3B0
    def Hvec():
        c = np.zeros(n); c[0] = -1; c[1] = -1; return c, Mp
    def Tdvec():
        c = np.zeros(n); c[0] = -2; c[1] = -2; c[2] = 2; c[6] = 2; c[7] = -1; c[8] = -3
        return c, 2 * Mp - Ns
    for ubr in ["v", "za", "avg"]:
        for wb in ["H", "dual"]:
            # objective: O/2 + W - S0 - B0 + za - u - lb - 2el/3 - (Td - y)/2
            c = np.zeros(n); const = 0.0
            c[0] += 0.5
            if wb == "H":
                hc, h0 = Hvec(); c += hc; const += h0
            else:
                c[3] += 1; c[4] += 1; const -= (w if mode == "soft" else 0.0)
            c[9] -= 1; c[8] -= 1; c[2] += 1
            if ubr == "v": c[3] -= 1
            elif ubr == "za": c[2] -= 1
            else: c[3] -= 1/3; c[2] -= 1/3
            c[4] -= 1; c[5] -= 2/3
            tc, t0 = Tdvec(); c -= tc / 2; const -= t0 / 2
            c[3] += 0.5; c[4] += 1.5; c[5] += 0.5
            A = []; b = []
            # H >= 0 : O + DH <= Mp
            r = np.zeros(n); r[0] = 1; r[1] = 1; A.append(r); b.append(Mp)
            # za <= lp
            r = np.zeros(n); r[2] = 1; A.append(r); b.append(lp)
            # y <= Td :  v + 3lb + el - Td <= 0
            r = np.zeros(n); r[3] = 1; r[4] = 3; r[5] = 1; r -= tc; A.append(r); b.append(t0)
            # 2A0 <= O
            r = np.zeros(n); r[6] = 2; r[0] = -1; A.append(r); b.append(0)
            # N0 + B0 <= A0 ; S0 <= A0
            r = np.zeros(n); r[7] = 1; r[8] = 1; r[6] = -1; A.append(r); b.append(0)
            r = np.zeros(n); r[9] = 1; r[6] = -1; A.append(r); b.append(0)
            # u branch validity: u = chosen value is <= others  (we maximize -u, so must ensure chosen is the min)
            uv = {"v": np.eye(n)[3], "za": np.eye(n)[2], "avg": (np.eye(n)[3] + np.eye(n)[2]) / 3}
            for other in uv:
                if other != ubr:
                    A.append(uv[ubr] - uv[other]); b.append(0)
            # W branch validity
            hc, h0 = Hvec()
            dual = np.zeros(n); dual[3] = 1; dual[4] = 1
            if mode == "soft":
                # W = max(H, v+lb-w)
                if wb == "H": A.append(dual - hc); b.append(h0 + w)
                else: A.append(hc - dual); b.append(-h0 - w)
            else:
                # hard: W = H if Td <= H + w, else max(H, v+lb)
                pass
            if fix_O0:
                r = np.zeros(n); r[0] = 1; A.append(r); b.append(0)
            if mode == "hard":
                # sub-branches: (Td <= H+w -> W=H) or (Td >= H+w -> W=max(H,v+lb))
                for sub in ["inside", "outside"]:
                    A2 = list(A); b2 = list(b)
                    if sub == "inside":
                        if wb != "H": continue
                        A2.append(tc - hc); b2.append(h0 + w - t0)
                    else:
                        A2.append(hc - tc); b2.append(t0 - h0 - w)
                        if wb == "H": A2.append(dual - hc); b2.append(h0)
                        else: A2.append(hc - dual); b2.append(-h0)
                    res = linprog(-c, A_ub=np.array(A2), b_ub=np.array(b2), bounds=[(0, None)] * n, method="highs")
                    if res.status == 0:
                        val = -res.fun + const
                        if val > best[0]: best = (val, (ubr, wb, sub, np.round(res.x, 5)))
                continue
            res = linprog(-c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * n, method="highs")
            if res.status == 0:
                val = -res.fun + const
                if val > best[0]: best = (val, (ubr, wb, np.round(res.x, 5)))
    return best[0] - Mp, best[1]

if __name__ == "__main__":
    w = float(sys.argv[1]); lx = eval(sys.argv[2]); ly = eval(sys.argv[3]); ell = eval(sys.argv[4])
    M = lx + ly
    print("geometry M=%.5f ell=%.5f M+ell-1=%.5f w=%.3f" % (M, ell, M + ell - 1, w))
    for d in np.linspace(0, ell, 7):
        Mp, lp = M - 2 * d, ell - d
        script_e = max(0, (1 + 2 * lp - 2 * Mp) / 3, Mp + lp - 1 - w)
        paper_e = max(0, (1 + 3 * lp - 2 * Mp) / 4) if 2 * Mp + lp >= 1 else 0
        e_soft, arg = emax(Mp, lp, w, "soft")
        e_hard, argh = emax(Mp, lp, w, "hard")
        e_O0, _ = emax(Mp, lp, w, "soft", fix_O0=True)
        print("d=%.4f  script e=%.5f  LP(soft W)=%.5f  LP(hard W)=%.5f  LP(O=0)=%.5f  paper(1+3l'-2M')/4=%.5f  -d+e/2: script %.5f LP %.5f"
              % (d, script_e, e_soft, e_hard, e_O0, paper_e, -d + script_e / 2, -d + max(0, e_soft) / 2))
        if d == 0: print("   argmax soft:", arg)
