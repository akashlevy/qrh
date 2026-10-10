"""Patterson bias check: S(N) = sum_{n sf, odd, N<N(n)<=2N} c(n) w(N(n)/N) for c = beta = gamma_2 conj(alpha),
gamma_2, gamma_2 alpha, and the Moebius partner a = mu conj G.  Random size ~ (sum w^2)^{1/2};
a Patterson pole gives ~ count^{5/6}."""
import numpy as np, time, json
from edge_lib import *
out = {}
fields, pnorm, plist = setup(130000)
for N0 in (2000, 4000, 8000, 16000, 32000, 64000):
    t0 = time.time()
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    w = C.w
    g2 = C.beta * C.alpha
    res = {}
    for name, c in (("beta", C.beta), ("gamma2", g2), ("gamma2*alpha", g2 * C.alpha), ("mob_partner", C.amob), ("mu", C.mu)):
        S = np.sum(c * w); rnd = np.sqrt(np.sum(w ** 2))
        res[name] = (abs(S) / rnd, float(np.angle(S)))
    out[N0] = dict(count=C.n, ratio=res)
    print(N0, C.n, "  ".join(f"{k}: |S|/sqrt(sum w^2)={v[0]:.2f}" for k, v in res.items()), f"  (count^(1/3)={C.n**(1/3):.1f}) {time.time()-t0:.0f}s", flush=True)
json.dump(out, open("bias_check.json", "w"), indent=1)
