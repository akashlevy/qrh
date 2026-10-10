"""
numerics: window / sector statistics and fits for the energy runs (data/energy_N*.npz).

Window K = N0/F: rows K/2 < N(k) <= K with V = sin^2(2 pi (N(k)/K - 1/2)); principal row excluded.
rho = sum V (|S_k|^2 - D_k) / sum V D_k   (D_k the exact n-diagonal).
Sector statistics: rows split into S angular sectors of arg(k) (independent sub-families); for each sector
rho_s; rms_S = sqrt(mean rho_s^2), mean_S.  Single-mark and per-mark (Gram) energies use rows with (P,k)=1.
Fits (OLS):  log rms_12 = a + b log(rows per sector) + c log F [+ d log L].
"""
import sys, os, json, glob, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SECTORS = (6, 12, 24, 48)


def load(path):
    Z = np.load(path, allow_pickle=False)
    meta = json.loads(str(Z["meta"]))
    d = {k: Z[k] for k in Z.files if k != "meta"}
    d["meta"] = meta
    d["names"] = [str(x) for x in d["names"]]
    d["fac"] = [set(int(x) for x in s.split(";")) if s else set() for s in d["rowfac"]]
    return d


def row_masks(d):
    """name -> boolean row mask (rows allowed), for single marks and Gram per-mark columns."""
    m = {}
    meta = d["meta"]
    for s in meta["singles"]:
        ok = np.array([s["idx"] not in f for f in d["fac"]])
        m[f"bP{s['N']}"] = ok; m[f"rP{s['N']}"] = ok
    for j, g in enumerate(meta["gram"]):
        m[f"g{j}"] = np.array([g["idx"] not in f for f in d["fac"]])
    return m


def window_stats(d, Fs=None, sectors=SECTORS):
    meta = d["meta"]; N0 = meta["N0"]
    rowN = d["rowN"].astype(float); S = d["S"]; D = d["D"]; names = d["names"]
    E = np.abs(S) ** 2 - D
    masks = row_masks(d)
    Mk = np.ones((len(rowN), len(names)))
    for j, n in enumerate(names):
        if n in masks: Mk[:, j] = masks[n]
    if Fs is None:
        Fs = [2 ** (j / 2) for j in range(0, 40)]
    secs = {ns: (np.floor((d["rowarg"] + np.pi) / (2 * np.pi) * ns).astype(int)) % ns for ns in sectors}
    out = []
    for F in Fs:
        K = N0 / F
        x = rowN / K
        inwin = (x > 0.5) & (x <= 1) & (rowN > 1)
        nrow = int(inwin.sum())
        if nrow < 6 * 12: break
        V = np.where(inwin, np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)
        Vm = V[:, None] * Mk
        num = (Vm * E).sum(0); den = (Vm * D).sum(0)
        w = dict(F=F, K=K, nrow=nrow, rho=dict(zip(names, (num / den).tolist())), sec={})
        for ns, sec in secs.items():
            nums = np.zeros((ns, len(names))); dens = np.zeros((ns, len(names)))
            np.add.at(nums, sec, Vm * E); np.add.at(dens, sec, Vm * D)
            cnt = np.bincount(sec[inwin], minlength=ns)
            ok = cnt > 0
            with np.errstate(invalid="ignore", divide="ignore"):
                r = nums[ok] / dens[ok]
            w["sec"][ns] = dict(rows_per_sector=float(nrow / ns), min_rows=int(cnt.min()),
                                rms=dict(zip(names, np.sqrt(np.nanmean(r ** 2, 0)).tolist())),
                                mean=dict(zip(names, np.nanmean(r, 0).tolist())),
                                maxabs=dict(zip(names, np.nanmax(np.abs(r), 0).tolist())))
        out.append(w)
    return out


def ols(X, y):
    X = np.column_stack([np.ones(len(y))] + list(X))
    beta, res, rk, sv = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    dof = max(1, len(y) - X.shape[1])
    s2 = float(r @ r) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    return beta, np.sqrt(np.diag(cov)), cov, math.sqrt(s2)


def collect(allstats, sel, ns=12, Fmin=1.0, minrps=8, pool=None):
    """rows: (N0, F, rows/sector, rms) for the names selected by sel(name) (pooled in quadrature over names)."""
    pts = []
    for N0, W in allstats.items():
        for w in W:
            if w["F"] < Fmin - 1e-9: continue
            s = w["sec"][ns]
            if s["rows_per_sector"] < minrps: continue
            vals = [s["rms"][n] for n in s["rms"] if sel(n)]
            if not vals: continue
            pts.append((N0, w["F"], s["rows_per_sector"], float(np.sqrt(np.mean(np.square(vals)))), len(vals)))
    return pts


def fit_rows_F(pts, extra=None):
    y = np.log([p[3] for p in pts])
    X = [np.log([p[2] for p in pts]), np.log([p[1] for p in pts])]
    if extra is not None: X.append(np.log(extra))
    return ols(X, y)


def main(files=None, out="analysis"):
    files = files or sorted(glob.glob(os.path.join(HERE, "data", "energy_N*.npz")),
                            key=lambda f: int(os.path.basename(f)[8:-4]))
    allstats = {}; metas = {}
    for f in files:
        d = load(f)
        N0 = d["meta"]["N0"]; metas[N0] = d["meta"]
        allstats[N0] = window_stats(d)
        print(f"{os.path.basename(f)}: cols {d['meta']['cols']} rows {d['meta']['rows']}, windows {len(allstats[N0])}", flush=True)
    json.dump(dict(stats={str(k): v for k, v in allstats.items()}, metas={str(k): v for k, v in metas.items()}),
              open(os.path.join(HERE, out + "_windows.json"), "w"))
    return allstats, metas


if __name__ == "__main__":
    main()
