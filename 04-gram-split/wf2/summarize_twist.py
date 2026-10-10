import numpy as np
edges = [7, 32, 128, 512, 2048, 8192, 24001]
for f in (3, 12):
    print("\n=== f = %d : mean Rat over nu (s.e. over nu) ===" % f)
    print("%-8s %-14s %-6s %-26s %-14s %-14s %-8s %-8s %-8s %-8s" % ("Q", "N(nu) bin", "#nu", "N^2/Q | N^3/Q | N^4/Q", "rho", "pi", "1", "gamma2", "zeta", "rand"))
    for Q in (1500, 3000, 6000):
        d = np.load("twist_Q%d.npz" % Q)
        N = d["nu_N"]
        for b in range(len(edges) - 1):
            s = (N >= edges[b]) & (N < edges[b + 1])
            if not s.sum(): continue
            mid = np.exp(np.log(N[s]).mean())
            line = "%-8d [%5d,%5d)  %-6d %-26s " % (Q, edges[b], edges[b + 1], s.sum(), "%.2g | %.2g | %.2g" % (mid ** 2 / Q, mid ** 3 / Q, mid ** 4 / Q))
            for nm in ("rho", "pi"):
                r = d["%s_f%d" % (nm, f)][s]; line += "%-14s " % ("%.2f(%.2f)" % (r.mean(), r.std() / np.sqrt(len(r))))
            for nm in ("1", "gamma2", "zeta", "rand"):
                line += "%-8s " % ("%.2f" % d["%s_f%d" % (nm, f)][s].mean())
            print(line)
    # pooled over Q, N(nu) < 128
    for nm in ("rho", "pi", "1", "conj(alpha)^2", "gamma2", "zeta", "rand"):
        v = np.concatenate([np.load("twist_Q%d.npz" % Q)["%s_f%d" % (nm, f)][np.load("twist_Q%d.npz" % Q)["nu_N"] < 128] for Q in (1500, 3000, 6000)])
        print("  pooled N(nu)<128, all Q: %-14s %.2f (s.e. %.2f, n=%d)" % (nm, v.mean(), v.std() / np.sqrt(len(v)), len(v)))
