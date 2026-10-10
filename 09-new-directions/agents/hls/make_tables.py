"""hls: markdown rows for HLS.md from runs/scan*.json"""
import json, sys
for fn in sys.argv[1:]:
    for r in json.load(open(fn)):
        c = r["cfg"]
        hb = f"{r['quadref']['lam_over_sum']:.3f}" if "quadref" in r else "—"
        print(f"| {c['N0']} | {c['K']} | {c['L0']} | {r['nc']} | {r['nk']}×{r['nP']} = {r['nrows']} | "
              f"{r['lam_over_sum']:.3f} | {r['lam'][0]/r['lam_prod']:.4f} | {r['lam'][0]/r['rmt']:.3f} | "
              f"{r['seq_bound']/r['lam'][0]:.1f} | {hb} |")
