"""Aggregate h2 runs into runs/summary.txt (tables quoted in H2.md)."""
import json, glob, os
import numpy as np

D = "/home/user/qrh/09-new-directions/agents/h2/runs"
out = []
P = out.append

P("== payoff: theta vs (H2) saving theta2 = 1 - Hoff  ((H1)-average assumed, bias piece kept; model_h2.py) ==")
rows = []
for f in glob.glob(f"{D}/hoff_*.json"):
    r = json.load(open(f)); rows.append((r["cfg"]["Hoff"], r["theta"], r["geom"]))
rows.sort()
for h, th, g in rows:
    best = min(th, 0.8749568)
    M = g["lx"] + g["ly"]
    P(f"  theta2 = {1-h:.2f}  (Hoff {h:.2f}): offset-branch {th:.6f}  -> theta = {best:.6f}   geometry M={M:.3f} ell={g['ell']:.3f}  M+ell-1={M+g['ell']-1:.3f}  ell<=M/2: {g['ell']<=M/2}")
hs = np.array([r[0] for r in rows]); ts = np.array([r[1] for r in rows])
a, b = np.polyfit(hs[hs <= 0.13], ts[hs <= 0.13], 1)
P(f"  linear fit theta ~ {b:.6f} + {a:.4f} Hoff ;  break-even with 0.874957 at Hoff = {(0.8749568-b)/a:.4f};  27/31 at Hoff = {(27/31-b)/a:.4f}")

for N0 in (20000, 40000):
    f = f"{D}/struct_N{N0}.json"
    if not os.path.exists(f): continue
    r = json.load(open(f))
    P(f"\n== structure hunt N0 = {N0} (cols {r['ncols']}, rows {r['nrows']}) : pure off-diagonal O = (G - Delta)/sqrt(G11 G22) ==")
    P("  block      K      rows  L'   rmsO th/rnd      rmsC th   Delta   ineqfail  indef  <z2_O> th/rnd  zmax_O th(arg)            cubic1 th/rnd  conj O th/rnd  slopeN th   Tcoh(p78)/sd")
    for b in r["blocks"]:
        t, q = b["theta"], b["random"]
        sd = t["T_rand_mean_sd_max"][1]
        P(f"  {b['label']:9s} {b['K']:<6d} {b['rows']:<5d} {b['marks']:<4d} {t['rms_O']:.4f}/{q['rms_O']:.4f}  {t['rms_C']:.4f}  {t['rms_Delta_rel']:.4f}  "
          f"{t['frac_ineq_fails']:.2f}      {t['frac_2x2_indefinite']:.2f}   {t['z2_mean_O']:.2f}/{q['z2_mean_O']:.2f}      {t['z_max_O']:.2f}({t['z_argmax_O']:12s})   "
          f"{(t['cubic1_O'] or 0):.2f}/{(q['cubic1_O'] or 0):.2f}      {(t['conjpair_O'] or 0):.2f}/{(q['conjpair_O'] or 0):.2f}      {(t['slopeN_O'] or 0):+.3f}     {t['T_coh']['p78']/sd:+.2f}")
    # aggregate z^2 per candidate over the L0 blocks (disjoint mark sets; windows ~independent after removing Delta)
    blocks = [b for b in r["blocks"] if b["label"] != "all"]
    keys = list(blocks[0]["theta"]["z_O"].keys())
    zt = np.array([[b["theta"]["z_O"][k] ** 2 for k in keys] for b in blocks])
    zr = np.array([[b["random"]["z_O"][k] ** 2 for k in keys] for b in blocks])
    P(f"  phase tests (O): mean z^2 over {len(blocks)} blocks x {len(keys)} phases: theta {zt.mean():.3f}, random {zr.mean():.3f} (null 1);"
      f" largest per-phase mean z^2: theta {zt.mean(0).max():.2f} ({keys[int(zt.mean(0).argmax())]}), random {zr.mean(0).max():.2f} ({keys[int(zr.mean(0).argmax())]})")
    cub = [b["theta"]["cubic1_O"] for b in blocks if b["theta"]["cubic1_O"]]
    cubr = [b["random"]["cubic1_O"] for b in blocks if b["random"]["cubic1_O"]]
    cj = [b["theta"]["conjpair_O"] for b in blocks if b["theta"]["conjpair_O"]]
    cjr = [b["random"]["conjpair_O"] for b in blocks if b["random"]["conjpair_O"]]
    P(f"  pair classes (O): (P1/P2)_3 = 1 vs all: theta {np.mean(cub):.3f}, random {np.mean(cubr):.3f};  conjugate pairs: theta {np.mean(cj):.2f}, random {np.mean(cjr):.2f}")
    byL = {}
    for b in blocks:
        byL.setdefault(b["label"], []).append((b["K"], b["theta"]["rms_O"], b["random"]["rms_O"], b["rows"]))
    P("  L-dependence of rms O at fixed K (theta / random):")
    Ks = sorted({b["K"] for b in blocks})
    for K in Ks:
        P(f"    K={K:<6d}: " + "  ".join(f"{lab}: {th:.4f}/{rn:.4f}" for lab, v in byL.items() for (k, th, rn, nr) in v if k == K)
          + f"   (1/sqrt(rows) = {1/np.sqrt([nr for v in byL.values() for (k, th, rn, nr) in v if k == K][0]):.4f})")
    for p in r["persistence"]:
        for w in p["windows"]:
            P(f"  persistence K={w['K']}: sector coherence ratio th {w['theta_sector_coherence_ratio']:.3f} rnd {w['random_sector_coherence_ratio']:.3f} (null 1, deterministic 12);"
              f" mean inter-sector corr th {w['theta_sector_corr_mean']:+.4f} rnd {w['random_sector_corr_mean']:+.4f}; column-halves corr {w['theta_halves_corr'][0]:+.4f}")
        P("  dyadic-window corr (O): theta " + " ".join(f"{d['K']}:{d['corr_re']:+.3f}" for d in p["dyadic"]["theta"]) +
          " | random " + " ".join(f"{d['K']}:{d['corr_re']:+.3f}" for d in p["dyadic"]["random"]))
    for tw in r["twisted"]:
        for nm in ("theta", "random"):
            T = np.array(tw[nm])
            P(f"  cocycle twists chi_k(P)^a chi_P(k)^b, rms of normalised G (K={tw['K']}, {nm}):")
            for a in range(6):
                P("     a=%d: " % a + " ".join(f"{T[a, bb]:.4f}" for bb in range(6)))

G = f"{D}/struct_N20000_G.json"
for f in sorted(glob.glob(f"{D}/dft_*.json")):
    r = json.load(open(f))
    P(f"\n== all-character spectrum in k mod P1P2 ({os.path.basename(f)}; off-diagonal part, per-row n-diagonal removed) ==")
    for nm in ("theta", "random"):
        s = r["summary_" + nm]
        P(f"  {nm:7s}: pairs {s['pairs']}, |O12|^2/mean|F|^2 = {s['triv_over_mean_mean']:.3f} (null 1), rank of trivial char {s['triv_rank_mean']:.3f} (null 0.5),"
          f" max/mean {s['max_over_mean_mean']:.2f} (iid {s['expected_max_iid']:.2f}), argmax at sextic char {s['argmax_sextic_frac']:.3f}, sextic points mean {s['sextic_mean']:.3f}")

for f in sorted(glob.glob(f"{D}/coh_*.json")):
    r = json.load(open(f))
    P(f"\n== coherent vectors ({os.path.basename(f)}): t = T(c)/sd_random(T), {r['summary']['theta']['one']['n']} independent blocks; null mean t^2 = 1 +- {np.sqrt(2/r['summary']['theta']['one']['n']):.2f} ==")
    for k in r["summary"]["theta"]:
        a, b = r["summary"]["theta"][k], r["summary"]["random"][k]
        P(f"  {k:9s} theta: mean t {a['mean_t']:+.3f} mean t^2 {a['mean_t2']:.3f}   | random columns: mean t {b['mean_t']:+.3f} mean t^2 {b['mean_t2']:.3f}")

f = f"{D}/conjfrac.json"
if os.path.exists(f):
    r = json.load(open(f))
    P("\n== conjugation symmetry explains the T(c) excess along structured (c_conjP = conj c_P) vectors: E t^2 = 1 + f ==")
    for N0, v in r.items():
        P(f"  N0={N0}: <f> = {v['f_mean']:.3f}, predicted <t^2> = {v['pred_t2']:.3f}, observed theta {v['obs_t2_theta']:.3f}, random columns {v['obs_t2_random']:.3f};"
          f" per-block fit t^2 = {v['fit_t2_vs_f'][0]:.2f} + {v['fit_t2_vs_f'][1]:.2f} f (prediction 1 + 1 f)")

f = f"{D}/bias_N20000.json"
if os.path.exists(f):
    r = json.load(open(f))
    P(f"\n== Patterson bias of slot sums X_c(n) = sum_P c_P chi_P(n)^-2: max |corr| with conj(gamma_2(n)) alpha(n)^j dictionary (null {r['null']:.4f}) ==")
    for row in r["rows"]:
        P(f"  L0={row['L0']:<5d} L'={row['marks']:<4d} c={row['vec']:9s} max|r| {row['best_abs']:.4f} at {row['best']}")

for f in sorted(glob.glob(f"{D}/kernel_*.json")):
    r = json.load(open(f)); c = r["flat_counterexample"]
    P(f"\n== idea 1 kernel test ({os.path.basename(f)}): N={r['N']}, K={r['K']}, K'={r['Kp']} rows ==")
    P(f"  rank deficiency of Gamma {r['nker']} (N-K' = {r['N_minus_Kp']}); Q eigenvalues [{r['Q_min']:.1f}, {r['Q_max']:.1f}], {100*r['Q_frac_negative']:.1f}% negative, "
      f"median negative {r['Q_median_negative']:.1f} vs mean diagonal {r['mean_diag']:.1f}; mark mass in ker Gamma {np.mean(r['mark_mass_in_kernel']):.3f} (1-K'/N = {r['predicted_kernel_mass']:.3f})")
    P(f"  flat certificate: rms | |x2|-|x1| | / rms|x1| = {c['rms_rel_modulus_dev']:.1e};  O11 = O22 = {c['O11']:.1f}, |Delta12| = {c['Delta12_abs']:.1e}, |O12| = {c['O12_abs']:.1f} = G11;"
      f"  |O12|/sqrt(O11 O22) = {c['ratio_O12_over_sqrt_O11O22']:.1f},  |O12|/Delta11 = {c['ratio_O12_over_diag']:.2f}")

open(f"{D}/summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
