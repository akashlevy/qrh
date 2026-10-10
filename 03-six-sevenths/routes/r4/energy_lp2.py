"""
Independent re-derivation of the low-side energy loss of Lemma probe-row-norm
(7/8 paper, Lemma reflected-energy, eq. old-eq:4.6), as an LP over block data.

E_ref = cnt + max(H, v+lb [-w]) - S0 - B0 + za - u - lb - (Td - y)/2        (e_lambda = 0)
   u  = min(v, za, (v+za)/3)
   Td = 2H + 2A0 + 2za - Ns - N0 - 3B0,  Ns = 1 + l'   (completed length)
   y  = v + 3lb <= Td        (retained dyads; others are negligible)
constraints:  H <= M' - O,  A0 = Anz + Az,  O >= 2Anz + 6Az  (powerful primes: v_p>=2, j=0 actives have v_p>=6),
              N0 + B0 + S4 <= Anz (j=4 primes split), S0 = S4 + Az,  za <= l'
count of frozen powerful parts:
     paper:    cnt = O/2                   (every powerful ideal is x^2 y^3)
     refined:  cnt = O/6 + 2*Anz/3         (sixth-power part counted at O0/6, the rest has radical Anz)
Returns max(E_ref) - M'  (the energy loss e(M',l')).
"""
import itertools, numpy as np
from scipy.optimize import linprog

names = "O H Anz Az N0 B0 S4 za v lb".split()
ix = {k: i for i, k in enumerate(names)}
n = len(names)

def vec(**kw):
    r = np.zeros(n)
    for k, c in kw.items():
        r[ix[k]] += c
    return r

def emax(Mp, lp, w=0.0, count="paper", extra=None, ret_arg=False):
    Ns = 1 + lp
    # Td = 2H + 2A0 + 2za - Ns - N0 - 3B0 , A0 = Anz + Az
    Td = vec(H=2, Anz=2, Az=2, za=2, N0=-1, B0=-3); Td0 = -Ns
    y = vec(v=1, lb=3)
    best = (-99, None)
    for ub in ("v", "za", "avg"):
        for wb in ("H", "dual"):
            c = np.zeros(n); const = 0.0
            if count == "paper":
                c += vec(O=0.5)
            elif count == "refined":
                c += vec(O=1/6, Anz=2/3)
            elif count == "sixth":   # hypothetical: only sixth powers counted (not a theorem)
                c += vec(O=1/6)
            if wb == "H":
                c += vec(H=1)
            else:
                c += vec(v=1, lb=1); const -= w
            c += vec(S4=-1, Az=-1, B0=-1, za=1, lb=-1)          # -S0 - B0 + za - lb
            uv = {"v": vec(v=1), "za": vec(za=1), "avg": vec(v=1/3, za=1/3)}
            c -= uv[ub]
            c -= 0.5 * (Td - y); const -= 0.5 * Td0
            A, b = [], []
            def le(row, rhs):
                A.append(row); b.append(rhs)
            le(vec(H=1, O=1), Mp)                              # H <= M' - O
            le(vec(Anz=2, Az=6, O=-1), 0)                      # O >= 2Anz + 6Az
            le(vec(N0=1, B0=1, S4=1, Anz=-1), 0)               # j=4 primes are among Anz
            le(vec(za=1), lp)
            le(y - Td, Td0)                                    # y <= Td  <=>  y - Td_lin <= Td0
            for o in uv:
                if o != ub:
                    le(uv[ub] - uv[o], 0)
            if wb == "H":
                le(vec(v=1, lb=1) - vec(H=1), w)
            else:
                le(vec(H=1) - vec(v=1, lb=1), -w)
            if extra:
                for row, rhs in extra:
                    le(row, rhs)
            res = linprog(-c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * n, method="highs")
            if res.status == 0:
                val = -res.fun + const
                if val > best[0]:
                    best = (val, (ub, wb, dict(zip(names, np.round(res.x, 5)))))
    if ret_arg:
        return best[0] - Mp, best[1]
    return best[0] - Mp

if __name__ == "__main__":
    print("check against closed forms: paper count -> max(0,(1+3l'-2M')/4, M'+l'-1)")
    for (Mp, lp) in [(5/6, 1/6), (0.8, 0.2), (0.75, 0.25), (6/7, 2/7), (0.9, 0.2), (0.7, 0.3), (0.85, 0.25)]:
        for cnt in ("paper", "refined", "sixth"):
            e, arg = emax(Mp, lp, 0.0, cnt, ret_arg=True)
            print(f"M'={Mp:.4f} l'={lp:.4f} count={cnt:8s}: e={e:+.5f}   caseA={(1+3*lp-2*Mp)/4:+.5f} "
                  f"O=0 form={(1+2*lp-2*Mp)/3:+.5f} wall={Mp+lp-1:+.5f}   arg={arg}")
