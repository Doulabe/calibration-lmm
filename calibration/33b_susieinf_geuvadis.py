# ==============================================================================
# 33b_susieinf_geuvadis.py
# SuSiE-inf (official fine-mapping-inf code, settings of its driver run_fine_mapping.py,
# as in 31_susieinf_selection.py) on the GEUVADIS application. FINEMAP-inf was dropped here:
# on windows with very strong leads its shotgun search exceeded 2 GB and several hours per window.
# 8 cis windows, EUR (n = 358) and YRI (n = 87-89), in-sample LD.
# Input: results/geuvadis_comp/data/<pop>_<symbol>_{Z,y,pos}.txt (33a_export_geuvadis.R)
# Output: results/geuvadis_comp/susieinf.jsonl (one line per window; positions are 1-based indices)
# Usage (from the repository root): python3 calibration/33b_susieinf_geuvadis.py
# ==============================================================================
import os, sys, json, glob, time
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool

ROOT = os.path.abspath(".")
FMI = os.environ.get("FINEMAPINF_DIR", os.path.join(ROOT, "comparators", "fine-mapping-inf-master"))
sys.path.insert(0, os.path.join(FMI, "susieinf")); sys.path.insert(0, os.path.join(FMI, "finemapinf"))
import susieinf, finemapinf
RESDIR = os.environ.get("GEU_RESDIR", "geuvadis_comp")       # geuvadis_comp_cov: data residualized on the covariates (33f)
DATA = os.path.join(ROOT, "calibration", "results", RESDIR, "data")
OUT = os.path.join(ROOT, "calibration", "results", RESDIR, "susieinf.jsonl")

def one(pre):
    Z = np.loadtxt(pre + "_Z.txt"); y = np.loadtxt(pre + "_y.txt"); pos = np.loadtxt(pre + "_pos.txt", dtype=int)
    n = Z.shape[0]
    X = Z - Z.mean(0); X = X / np.sqrt((X ** 2).mean(0))
    y = (y - y.mean()) / np.sqrt(((y - y.mean()) ** 2).mean()); z = X.T @ y / np.sqrt(n)
    Dsq, V = np.linalg.eigh(X.T @ X); Dsq = np.maximum(Dsq, 0)
    t0 = time.time()
    s = susieinf.susie(z, 1.0, n, 10, V=V, Dsq=Dsq, method="moments", maxiter=100, PIP_tol=1e-3, verbose=False)
    conv = bool(s.get("converged", True))
    if not conv:
        s = susieinf.susie(z, 1.0, n, 10, V=V, Dsq=Dsq, est_tausq=False, method="moments", maxiter=100, PIP_tol=1e-3, verbose=False)
    cs = susieinf.cred(s["PIP"], coverage=0.95, purity=0.5, V=V, Dsq=Dsq, n=n)
    pip = 1 - np.prod(1 - s["PIP"], axis=1)
    so = np.argsort(-pip)[:5]
    return {"window": os.path.basename(pre), "n": int(n), "m": int(Z.shape[1]), "susieinf_converged": conv,
            "tausq": float(s["tausq"]), "sigmasq": float(s["sigmasq"]),
            "susieinf_cs": [[int(i) + 1 for i in c] for c in cs],
            "susieinf_top": [int(c[int(np.argmax(pip[list(c)]))]) + 1 for c in cs],
            "susieinf_top5": [int(i) + 1 for i in so], "susieinf_top5_pip": [float(pip[i]) for i in so],
            "secs": time.time() - t0}

if __name__ == "__main__":
    pres = sorted({f[:-6] for f in glob.glob(os.path.join(DATA, "*_Z.txt"))})
    with Pool(1, maxtasksperchild=1) as pool, open(OUT, "w") as fh:
        for r in pool.imap(one, pres):
            fh.write(json.dumps(r) + "\n"); fh.flush(); print(r["window"], "%.0f s" % r["secs"], r["susieinf_top"], flush=True)
