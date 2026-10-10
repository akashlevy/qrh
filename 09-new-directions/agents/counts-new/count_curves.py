"""Count curves R(delta, x) for the paper pieces and each NEW piece (counts_lab), at x = 0 and x = 1/2.
Prints the first delta where each new piece lowers R by > 1e-6, and R at the binding bins."""
import numpy as np, json, sys
import counts_lab as CL
import model_sa as MS
CL.install()

PAPER = dict(MS.DEFAULT)
DIAG = dict(PAPER, caseA="refined", delta0=0.001, mk=1.0, w2scale=10.0, cP=1.0)
NEW = {
    "GM+copies (heur.)": dict(gm=True),
    "Huxley+copies (heur.)": dict(hux=True),
    "HM, GLH off-diag (cond.)": dict(hm=True),
    "GM+Hux+HM": dict(gm=True, hux=True, hm=True),
    "Moebius 4th moment diag. H>=D^2 (hyp.)": dict(hk=2),
    "Moebius 6th moment diag. H>=D^3 (hyp.)": dict(hk=3),
    "plain 6th moment, m<=1/2 (beyond edge)": dict(p6=True, m6=0.5),
    "plain 6th moment, m<=1/3 (within edge)": dict(p6=True, m6=1 / 3),
    "plain 8th moment, m<=1/2": dict(p2k=4, mk2=0.5),
}
ds = np.round(np.arange(0.20, 0.76, 0.01), 4)

def curve(base, extra, x, S):
    cfg = dict(base); cfg.update(CL.NEWDEF); cfg.update(extra)
    return MS.R_of(ds, np.full_like(ds, x), S, cfg)

out = {}
for bname, base, S in [("paper", PAPER, 0.2055), ("DIAG", DIAG, 0.2308)]:
    for x in (0.0, 0.5):
        r0 = curve(base, {}, x, S)
        out[f"{bname} x={x}"] = dict(base=r0.tolist())
        print(f"\n=== {bname}, x = {x}  (S = {S})")
        print("delta: " + " ".join(f"{d:.2f}" for d in ds[::5]))
        print("R    : " + " ".join(f"{v:.4f}" for v in r0[::5]))
        for name, extra in NEW.items():
            rn = curve(base, extra, x, S)
            gain = r0 - rn
            k = np.where(gain > 1e-6)[0]
            first = ds[k[0]] if len(k) else None
            at = {d: float(gain[np.argmin(abs(ds - d))]) for d in (0.27, 0.33, 0.39, 0.42)}
            print(f"  {name:42s} first gain at delta = {first}; gain at .27/.33/.39/.42 = "
                  + "/".join(f"{at[d]:.4f}" for d in (0.27, 0.33, 0.39, 0.42)))
            out[f"{bname} x={x}"][name] = rn.tolist()
json.dump(dict(deltas=ds.tolist(), curves=out), open("count_curves.json", "w"))
# analytic check: R0 closed form at x = 0 (paper)
print("\nR0 closed form check (x=0):", [round(CL.R0(d), 5) for d in (0.27, 0.33, 0.39)])
