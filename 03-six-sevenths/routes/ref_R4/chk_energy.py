import sys, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/wf67/r4')
from energy_lp3 import emax
def closed4(Mp, lp, w, count):
    a = (1 + 3*lp - 2*Mp)/4 if count == "paper" else (1 + 2*lp - 2*Mp)/3
    return max(0.0, a, Mp + lp - 1 - w, (Mp + 4*lp - 2)/3 - w)
def closed3(Mp, lp, w, count):
    a = (1 + 3*lp - 2*Mp)/4 if count == "paper" else (1 + 2*lp - 2*Mp)/3
    return max(0.0, a, Mp + lp - 1 - w)
rng = np.random.default_rng(11)
for count in ("paper","refined"):
    md4=md3=0; arg4=arg3=None; ninf=0; nover=0
    for _ in range(3000):
        Mp = rng.uniform(0.3, 1.3); lp = rng.uniform(0, 0.7); w = rng.choice([0, 0.05, 0.1, 1/7, 0.2, 2/7, 0.3, 4/7])
        e = emax(Mp, lp, w, count)
        if e < -50: ninf+=1; continue
        c4 = closed4(Mp, lp, w, count); c3 = closed3(Mp,lp,w,count)
        if e - c4 > 1e-7: nover+=1
        if abs(e-c4)>md4: md4=abs(e-c4); arg4=(Mp,lp,w,e,c4)
        if abs(e-c3)>md3: md3=abs(e-c3); arg3=(Mp,lp,w,e,c3)
    print(count, "infeasible:",ninf, "LP>closed4 count:", nover, "max|LP-closed4|=",md4, arg4)
    print(count, "max|LP-closed3|=",md3, arg3)
