import json, sys, numpy as np
for fn in sys.argv[1:]:
    for line in open(fn):
        r = json.loads(line); th = r["theta"]; co = r.get("control", {})
        print(f"q={r['q']:2d} k={r['kappa']} l={r['lam']} nu={r['nu']} D={r['D']} win={int(r['window'])} rows={th['rows']:6d}{'s' if r['subset'] else ' '} "
              f"rho_th={th['rho_mean']:+.4f} (rms {th['rho_rms']:.4f}, min {min(th['rho']):+.4f} max {max(th['rho']):+.4f}) "
              f"ctrl={np.mean(co['rho']):+.4f} rms {co['rho_rms']:.4f} | m4 {np.nanmean(th['m4']):.2f}/{np.nanmean(co['m4']):.2f} | "
              f"cross {th['cross_rms'] or 0:.4f}/{co['cross_rms'] or 0:.4f} 1/sqrt(rows)={1/np.sqrt(th['rows']):.4f} | Ncoh {th['N_coherent']:+.4f} Nrnd {th['N_random']:+.4f} | op {th['opnorm2_over_diag']:.3f}/edge {th['rmt_edge']:.3f} | {r['secs']:.0f}s")
