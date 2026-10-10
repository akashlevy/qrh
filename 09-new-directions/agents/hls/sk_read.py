"""reader for sextic_kubota.c output"""
import numpy as np
def read(fn):
    with open(fn, "rb") as f:
        T, NP = np.frombuffer(f.read(8), np.int32)
        tp = np.frombuffer(f.read(16 * T), np.int64).reshape(T, 2)
        prec = np.dtype([("N", np.int64), ("a", np.int64), ("b", np.int64), ("g", np.float64, (5, 2))])
        pr = np.frombuffer(f.read(prec.itemsize * NP), prec)
        rec = np.dtype([("N", np.int32), ("mu", np.int32), ("g", np.float64, (5, 2)), ("e", np.int8, (T,))])
        data = np.frombuffer(f.read(), rec)
    g = data["g"][..., 0] + 1j * data["g"][..., 1]          # (n, 5): gamma_m, m = 1..5
    pg = pr["g"][..., 0] + 1j * pr["g"][..., 1]
    return dict(T=int(T), tp=tp, primes=pr, pg=pg, N=data["N"].astype(np.int64), mu=data["mu"], g=g, e=data["e"])
