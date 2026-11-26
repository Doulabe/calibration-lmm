# ==============================================================================
# 31_susieinf_selection.py
# SuSiE-inf and FINEMAP-inf (Cui et al. 2024, fine-mapping-inf, official code in
# comparators/fine-mapping-inf-master) on the datasets of the selection-level
# study exported by 29_export_selection_datasets.R.
# In-sample LD: columns standardized so that X'X/n has unit diagonal, y standardized
# (meansq = 1), z = X'y / sqrt(n); V, Dsq = eigendecomposition of X'X, once per window.
# Settings are those of the official driver run_fine_mapping.py: SuSiE-inf with L = 10,
# method of moments, 95% credible sets with purity 0.5 (falls back to tau^2 = 0 when it
# does not converge, as the driver does); FINEMAP-inf with SuSiE-inf's sigma^2 and tau^2,
# S = [0.05^2], sched_sss = [10000], nmodels = 5000.
# Two passes:
#   --susie    SuSiE-inf on all 1,800 datasets
#   --finemap  FINEMAP-inf on the 1,000 main datasets (200 global nulls + 800 with causals);
#              one fresh worker per dataset, because the hash tables of the C++ extension
#              are not released between calls (about 0.4 GB per call at m = 2,699).
# Output: results/comparators/susieinf/<gene>_<pass>.jsonl, one line per dataset
# Usage (from the repository root): python3 calibration/31_susieinf_selection.py --susie|--finemap [workers]
# ==============================================================================
import os, sys, json, time, re
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import pandas as pd
from multiprocessing import Pool

ROOT = os.path.abspath(".")
FMI = os.environ.get("FINEMAPINF_DIR", os.path.join(ROOT, "comparators", "fine-mapping-inf-master"))
sys.path.insert(0, os.path.join(FMI, "susieinf")); sys.path.insert(0, os.path.join(FMI, "finemapinf"))
import susieinf, finemapinf

CROOT = os.environ.get("COMP_ROOT", "comparators")           # comparators_purenull: y = noise (29b)
DATA = os.path.join(ROOT, "calibration", "results", CROOT, "data")
OUT = os.path.join(ROOT, "calibration", "results", CROOT, "susieinf")
GENES = ["ERAP2", "RPS26", "SLFN5", "PEX6-region"]
G = {}

def load(gene):
    """Standardized genotypes, phenotypes and the eigendecomposition of X'X (cached in an .npz)."""
    cache = os.path.join(OUT, gene + "_eigen.npz")
    Y = np.loadtxt(os.path.join(DATA, gene + "_Y.txt"))
    meta = pd.read_csv(os.path.join(DATA, gene + "_meta.csv"))
    if os.path.exists(cache):
        e = np.load(cache); X, V, Dsq = e["X"], e["V"], e["Dsq"]
    else:
        Z = np.loadtxt(os.path.join(DATA, gene + "_Z.txt"))
        X = Z - Z.mean(0); X = X / np.sqrt((X ** 2).mean(0))
        Dsq, V = np.linalg.eigh(X.T @ X); Dsq = np.maximum(Dsq, 0)
        np.savez(cache, X=X, V=V, Dsq=Dsq)
    return dict(X=X, Y=Y, meta=meta, n=X.shape[0], V=V, Dsq=Dsq)

def init_worker(gene):                               # spawn start method (macOS): load in each worker
    G[gene] = load(gene)

def zscores(d, k):
    y = d["Y"][:, k]; y = (y - y.mean()) / np.sqrt(((y - y.mean()) ** 2).mean())
    return d["X"].T @ y / np.sqrt(d["n"])

def fit_susieinf(d, z):
    n = d["n"]
    s = susieinf.susie(z, 1.0, n, 10, V=d["V"], Dsq=d["Dsq"], est_ssq=True, ssq=None, ssq_range=(0, 1),
                       method="moments", maxiter=100, PIP_tol=1e-3, verbose=False)
    conv = bool(s.get("converged", True))
    if not conv:                                     # driver fallback: tau^2 = 0
        s = susieinf.susie(z, 1.0, n, 10, V=d["V"], Dsq=d["Dsq"], est_tausq=False, est_ssq=True, ssq=None,
                           ssq_range=(0, 1), method="moments", maxiter=100, PIP_tol=1e-3, verbose=False)
    return s, conv

def run_susie(args):
    gene, k = args; d = G[gene]; z = zscores(d, k); t0 = time.time()
    s, conv = fit_susieinf(d, z)
    cs = susieinf.cred(s["PIP"], coverage=0.95, purity=0.5, V=d["V"], Dsq=d["Dsq"], n=d["n"])
    pip = 1 - np.prod(1 - s["PIP"], axis=1)
    return {"id": d["meta"]["id"][k], "cs": [[int(i) + 1 for i in c] for c in cs],
            "top": [int(c[int(np.argmax(pip[list(c)]))]) + 1 for c in cs],
            "pip_idx": [int(i) + 1 for i in np.nonzero(pip > 0.01)[0]], "pip": [float(v) for v in pip[pip > 0.01]],
            "tausq": float(s["tausq"]), "sigmasq": float(s["sigmasq"]), "converged": conv, "secs": time.time() - t0}

def run_finemap(args):
    gene, k = args
    if gene not in G: G[gene] = load(gene)
    d = G[gene]; z = zscores(d, k); t0 = time.time()
    s, conv = fit_susieinf(d, z)
    f = finemapinf.finemap(z, 1.0, d["n"], 10, V=d["V"], Dsq=d["Dsq"], sigmasq=s["sigmasq"], tausq=s["tausq"],
                           S=[0.05 ** 2], sched_sss=[10000], n_conv_sss=100, nmodels=5000,
                           prob_tol_sss=1e-3, seed=21, verbose=0)
    fp = np.asarray(f["PIP"])
    return {"id": d["meta"]["id"][k], "pip_idx": [int(i) + 1 for i in np.nonzero(fp > 0.01)[0]],
            "pip": [float(v) for v in fp[fp > 0.01]], "tausq": float(s["tausq"]), "secs": time.time() - t0}

def is_main(i):                                      # b <= 25, or a dataset with causals
    return bool(re.search(r"_K[13]_", i)) or int(i[-3:]) <= 25

if __name__ == "__main__":
    mode = "finemap" if "--finemap" in sys.argv else "susie"
    nums = [a for a in sys.argv[1:] if a.isdigit()]
    workers = int(nums[0]) if nums else 2
    os.makedirs(OUT, exist_ok=True)
    for gene in GENES:
        fn = os.path.join(OUT, "%s_%s.jsonl" % (gene, mode))
        done = set()
        if os.path.exists(fn):
            with open(fn) as fh: done = {json.loads(l)["id"] for l in fh if l.strip()}
        meta = pd.read_csv(os.path.join(DATA, gene + "_meta.csv"))
        ks = [k for k in range(len(meta)) if meta["id"][k] not in done and (mode == "susie" or is_main(meta["id"][k]))]
        print(gene, mode, "to do:", len(ks), flush=True)
        if not ks: continue
        if mode == "susie": load(gene)             # make sure the eigen cache exists before the workers start
        t0 = time.time()
        if mode == "susie":
            pool = Pool(workers, initializer=init_worker, initargs=(gene,))
        else:
            load(gene); pool = Pool(workers, maxtasksperchild=1)
        with open(fn, "a") as fh:
            for res in pool.imap_unordered(run_susie if mode == "susie" else run_finemap, [(gene, k) for k in ks],
                                           chunksize=1):
                fh.write(json.dumps(res) + "\n"); fh.flush()
        pool.close(); pool.join()
        print(gene, mode, "%d datasets in %.0f s" % (len(ks), time.time() - t0), flush=True)
