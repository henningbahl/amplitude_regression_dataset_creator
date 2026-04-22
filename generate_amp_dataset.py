#!/usr/bin/env python
"""Generate q q > Z + n g events with MadGraph, compute amplitudes via matrix2py.

Workflow follows https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/FAQ-General-4 :
  1. `output madevent` to generate LHE events with the requested cuts.
  2. `output standalone` of the same process to build matrix2py.so.
  3. Parse LHE, optionally apply post-cuts on m_inv and costheta, evaluate |M|^2.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import warnings
import xml.etree.ElementTree as ET

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MG_PATH_FILE = os.path.join(HERE, ".mg_path")
MG_BIN = None  # resolved in main() / run()


def _load_mg_bin():
    if os.path.isfile(MG_PATH_FILE):
        with open(MG_PATH_FILE) as f:
            p = f.read().strip()
            if p:
                return p
    return os.path.join(HERE, "MG", "bin", "mg5_aMC")


def _save_mg_bin(path):
    with open(MG_PATH_FILE, "w") as f:
        f.write(path.strip() + "\n")


def _resolve_mg_bin(cli_path):
    global MG_BIN
    if cli_path:
        MG_BIN = os.path.abspath(os.path.expanduser(cli_path))
        _save_mg_bin(MG_BIN)
        print(f"[generate] saved MadGraph path to {MG_PATH_FILE}: {MG_BIN}")
    else:
        MG_BIN = _load_mg_bin()
    if not os.path.isfile(MG_BIN):
        raise FileNotFoundError(
            f"MadGraph binary not found at {MG_BIN!r}. "
            f"Pass --mg-path /path/to/mg5_aMC to set it."
        )
    return MG_BIN


def _proc_line(process):
    p = process.strip()
    if not p.lower().startswith("generate "):
        p = "generate " + p
    return p


def _count_final_state(process):
    """Count final-state particles in a MadGraph process string."""
    after = process.split(">", 1)[1]
    n = 0
    for tok in after.split():
        if "=" in tok or tok[0] in "/$@,":
            break
        n += 1
    if n == 0:
        raise ValueError(f"Could not parse final-state particles from {process!r}")
    return n


def _slug(process):
    s = process.split(">", 1)[1] if ">" in process else process
    s = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")
    return s or "proc"


def _write_mg_script(path, mode, out_dir, process, group=False):
    with open(path, "w") as f:
        f.write("import model sm\n")
        if not group:
            f.write("set group_subprocesses False\n")
        f.write("define q = u d s c u~ d~ s~ c~\n")
        f.write(_proc_line(process) + "\n")
        extra = " --prefix=int" if mode == "standalone" else ""
        f.write(f"output {mode} {out_dir}{extra}\n")


def _log_call(cmd, log, cwd=None):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(cmd)}  (cwd={cwd or os.getcwd()})\n")
        f.flush()
        subprocess.check_call(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)


def _run_mg(script, log):
    _log_call([MG_BIN, "-f", script], log)


def _patch_run_card(run_card, n_events, dR=0.3):
    with open(run_card) as f:
        txt = f.read()

    def sub(key, value):
        nonlocal txt
        txt = re.sub(
            rf"^(\s*)\S+(\s*=\s*{re.escape(key)}\s*!.*)$",
            lambda m: f"{m.group(1)}{value}{m.group(2)}",
            txt,
            count=1,
            flags=re.MULTILINE,
        )

    sub("nevents", int(n_events))
    # max eta = -1 -> disable eta cut on jets
    sub("etaj", -1.0)
    sub("etab", -1.0)
    # delta R cut between jets
    sub("drjj", dR)
    # loosen ptj so the run actually proceeds
    sub("sde_strategy", 1)

    with open(run_card, "w") as f:
        f.write(txt)


_DUMMY_CUTS_TMPL = r"""      logical FUNCTION dummy_cuts(P)
      IMPLICIT NONE
      include 'genps.inc'
      include 'nexternal.inc'
      REAL*8 P(0:3,nexternal)
      integer i
      double precision psum0,psum1,psum2,psum3,m2,mi,ct
      double precision dt0,dt1,dt2,dt3,du0,du1,du2,du3,tval,uval
      LOGICAL  IS_A_J(NEXTERNAL),IS_A_L(NEXTERNAL)
      LOGICAL  IS_A_B(NEXTERNAL),IS_A_A(NEXTERNAL),IS_A_ONIUM(NEXTERNAL)
      LOGICAL  IS_A_NU(NEXTERNAL),IS_HEAVY(NEXTERNAL)
      logical  do_cuts(nexternal)
      COMMON /TO_SPECISA/IS_A_J,IS_A_A,IS_A_L,IS_A_B,IS_A_NU,IS_HEAVY,
     . IS_A_ONIUM, do_cuts

      dummy_cuts=.true.
      psum0=0d0
      psum1=0d0
      psum2=0d0
      psum3=0d0
      do i=3,nexternal
         psum0=psum0+P(0,i)
         psum1=psum1+P(1,i)
         psum2=psum2+P(2,i)
         psum3=psum3+P(3,i)
      enddo
      m2 = psum0*psum0-psum1*psum1-psum2*psum2-psum3*psum3
      if (m2.lt.0d0) m2=0d0
      mi = sqrt(m2)
      if (mi.lt.__MINV_LO__d0) then
         dummy_cuts=.false.
         return
      endif
      if (mi.gt.__MINV_HI__d0) then
         dummy_cuts=.false.
         return
      endif
__CT_BLOCK__
      return
      end
"""

_CT_BLOCK = r"""      if (nexternal.eq.4) then
c        Lorentz-invariant cos(theta*) = (u - t)/(u + t)
c        with t = (p1 - p3)^2, u = (p1 - p4)^2.
         dt0 = P(0,1) - P(0,3)
         dt1 = P(1,1) - P(1,3)
         dt2 = P(2,1) - P(2,3)
         dt3 = P(3,1) - P(3,3)
         tval = dt0*dt0 - dt1*dt1 - dt2*dt2 - dt3*dt3
         du0 = P(0,1) - P(0,4)
         du1 = P(1,1) - P(1,4)
         du2 = P(2,1) - P(2,4)
         du3 = P(3,1) - P(3,4)
         uval = du0*du0 - du1*du1 - du2*du2 - du3*du3
         ct = (uval - tval)/(uval + tval)
         if (ct.lt.__CT_LO__d0) then
            dummy_cuts=.false.
            return
         endif
         if (ct.gt.__CT_HI__d0) then
            dummy_cuts=.false.
            return
         endif
      endif
"""


def _patch_dummy_cuts(me_dir, m_inv_range, costheta_range):
    # Always rewrite dummy_cuts so reused directories don't retain stale cuts.
    lo, hi = m_inv_range or (None, None)
    lo = 0.0 if lo is None else float(lo)
    hi = 1.0e12 if hi is None else float(hi)
    body = _DUMMY_CUTS_TMPL.replace("__MINV_LO__", repr(lo)).replace(
        "__MINV_HI__", repr(hi)
    )
    if costheta_range is not None:
        clo, chi = costheta_range
        clo = -1.0 if clo is None else float(clo)
        chi = 1.0 if chi is None else float(chi)
        ct_block = _CT_BLOCK.replace("__CT_LO__", repr(clo)).replace(
            "__CT_HI__", repr(chi)
        )
    else:
        ct_block = ""
    body = body.replace("__CT_BLOCK__", ct_block)

    sub = os.path.join(me_dir, "SubProcesses")
    targets = [os.path.join(sub, "dummy_fct.f")]
    for d in os.listdir(sub):
        f = os.path.join(sub, d, "dummy_fct.f")
        if os.path.isfile(f):
            targets.append(f)
    for tgt in targets:
        with open(tgt) as fh:
            txt = fh.read()
        # replace from `logical FUNCTION dummy_cuts` up to its terminating `      end`
        new = re.sub(
            r"      logical FUNCTION dummy_cuts\(P\).*?\n      end\n",
            body,
            txt,
            count=1,
            flags=re.DOTALL,
        )
        with open(tgt, "w") as fh:
            fh.write(new)


def _generate_events(out_dir, n_events, dR, log):
    run_card = os.path.join(out_dir, "Cards", "run_card.dat")
    _patch_run_card(run_card, n_events, dR=dR)
    # disable MadAnalysis5
    ma5 = os.path.join(out_dir, "Cards", "madanalysis5_parton_card.dat")
    if os.path.exists(ma5):
        os.remove(ma5)
    ma5h = os.path.join(out_dir, "Cards", "madanalysis5_hadron_card.dat")
    if os.path.exists(ma5h):
        os.remove(ma5h)
    launch = os.path.join(HERE, "_launch.mg5")
    with open(launch, "w") as f:
        f.write(f"launch {out_dir}\n0\n")
    _run_mg(launch, log)
    os.remove(launch)
    # find produced LHE
    evt_dir = os.path.join(out_dir, "Events", "run_01")
    for name in os.listdir(evt_dir):
        if (
            name.endswith("unweighted_events.lhe.gz")
            or name == "unweighted_events.lhe.gz"
        ):
            lhe = os.path.join(evt_dir, name)
            subprocess.check_call(["gunzip", "-f", lhe])
            return lhe[:-3]
        if name.endswith(".lhe"):
            return os.path.join(evt_dir, name)
    raise RuntimeError("No LHE file found in " + evt_dir)


def _build_allmatrix2py(sa_dir, log):
    sub = os.path.join(sa_dir, "SubProcesses")
    env = os.environ.copy()
    env["MENUM"] = "2"
    # ensure f2py matches the current Python interpreter (not whatever 'f2py3' resolves to)
    env["F2PY"] = f"{sys.executable} -m numpy.f2py"
    # meson f2py compiles in a tmpdir; tell gfortran where to find coupl.inc etc.
    p_dirs = [
        d
        for d in os.listdir(sub)
        if d.startswith("P") and os.path.isdir(os.path.join(sub, d))
    ]
    inc_flags = " ".join([f"-I{os.path.join(sub, d)}" for d in p_dirs] + [f"-I{sub}"])
    env["FFLAGS"] = (env.get("FFLAGS", "") + " -fPIC " + inc_flags).strip()
    # remove stale non-PIC matrix.o / libmatrix.a from a previous build
    # (do NOT touch libdhelas.a / libmodel.a in standalone/lib)
    for d in p_dirs + [""]:
        full = os.path.join(sub, d) if d else sub
        for fn in os.listdir(full):
            if fn.endswith(".o"):
                os.remove(os.path.join(full, fn))
    libmatrix = os.path.join(sa_dir, "lib", "libmatrix.a")
    if os.path.exists(libmatrix):
        os.remove(libmatrix)
    # remove stale .so built against a different Python
    for root, _, files in os.walk(sub):
        for fn in files:
            if fn.endswith(".so") and "2py" in fn:
                os.remove(os.path.join(root, fn))
    with open(log, "a") as f:
        f.write(f"\n$ F2PY='{env['F2PY']}' make all_matrix2py.so  (cwd={sub})\n")
        f.flush()
        abs_lib = os.path.abspath(os.path.join(sa_dir, "lib"))
        subprocess.check_call(
            [
                "make",
                "all_matrix2py.so",
                f"LINKLIBS_ME=-L{abs_lib} -ldhelas -lmodel",
                f"LINKLIBS_ALL=-L{abs_lib} -lmatrix -ldhelas -lmodel",
            ],
            cwd=sub,
            stdout=f,
            stderr=subprocess.STDOUT,
            env=env,
        )
    return sub


_WORKER = r"""
import sys, json, os
sys.path.insert(0, sys.argv[1])
import matrix2py
matrix2py.initialise(sys.argv[2])
data = json.load(sys.stdin)
out = []
for pdgs, p in data:
    # transpose to fortran ordering
    p2 = [[p[i][j] for i in range(len(p))] for j in range(len(p[0]))]
    out.append(matrix2py.smatrixhel(pdgs, -1, p2, 0.13, 0.0, -1))
json.dump(out, sys.stdout)
"""


def _eval_amplitudes(p_dirs, param_card, events_by_pdir):
    """events_by_pdir: dict[p_dir] -> list of (pdgs, mom_list). Returns same keys -> list of me2."""
    import json, tempfile

    results = {}
    worker_path = os.path.join(tempfile.gettempdir(), "_mg_amp_worker.py")
    with open(worker_path, "w") as f:
        f.write(_WORKER)
    for p_dir, evts in events_by_pdir.items():
        payload = json.dumps([(pdgs, mom.tolist()) for pdgs, mom in evts])
        proc = subprocess.run(
            [sys.executable, worker_path, p_dir, param_card],
            input=payload,
            capture_output=True,
            text=True,
            check=True,
        )
        results[p_dir] = json.loads(proc.stdout)
    return results


def _parse_lhe(lhe_path):
    """Yield (pdgs, momenta[N,4]) per event. Only final-state ordering preserved;
    matrix2py needs the full initial+final state in MG ordering — we return all
    particles in the order they appear in the LHE event block, which matches MG."""
    with open(lhe_path) as f:
        text = f.read()
    # Wrap so ET can parse
    if not text.lstrip().startswith("<"):
        raise RuntimeError("Unexpected LHE")
    root = ET.fromstring(text)
    for event in root.iter("event"):
        lines = [ln for ln in event.text.strip().splitlines() if ln.strip()]
        header = lines[0].split()
        # LHE event header: NUP IDPRUP XWGTUP SCALUP AQEDUP AQCDUP
        nup = int(header[0])
        alphas = float(header[5])
        pdgs, mom = [], []
        for ln in lines[1 : 1 + nup]:
            parts = ln.split()
            pdgs.append(int(parts[0]))
            px, py, pz, E = (
                float(parts[6]),
                float(parts[7]),
                float(parts[8]),
                float(parts[9]),
            )
            mom.append([E, px, py, pz])
        yield pdgs, np.array(mom), alphas


def _invert(p):
    return [[p[i][j] for i in range(len(p))] for j in range(len(p[0]))]


def _post_cuts(pdgs, mom, m_inv_range, costheta_range):
    # invariant mass of all final-state particles (idx >= 2 in MG convention)
    pf = mom[2:]
    P = pf.sum(axis=0)
    m2 = P[0] ** 2 - (P[1] ** 2 + P[2] ** 2 + P[3] ** 2)
    m_inv = np.sqrt(max(m2, 0.0))
    if m_inv_range is not None:
        lo, hi = m_inv_range
        if (lo is not None and m_inv < lo) or (hi is not None and m_inv > hi):
            return False
    if costheta_range is not None:
        if len(pf) != 2:
            raise ValueError("costheta cut requires 2 -> 2 scattering")

        # Lorentz-invariant cos(theta*) = (u - t)/(u + t)
        def _sq(p):
            return p[0] ** 2 - p[1] ** 2 - p[2] ** 2 - p[3] ** 2

        p1, p3, p4 = mom[0], mom[2], mom[3]
        tval = _sq(p1 - p3)
        uval = _sq(p1 - p4)
        ct = (uval - tval) / (uval + tval)
        lo, hi = costheta_range
        if (lo is not None and ct < lo) or (hi is not None and ct > hi):
            return False
    return True


def _plot_hexbin(momenta, me2, n_final, out_path):
    import matplotlib.pyplot as plt

    if n_final != 2:
        raise ValueError("hexbin plot only available for 2->2 processes")

    # Mandelstams (Lorentz invariant): s=(p1+p2)^2, t=(p1-p3)^2, u=(p1-p4)^2.
    # cos(theta*) in the partonic CM frame = (u - t)/(u + t).
    def _sq(p):
        return p[:, 0] ** 2 - p[:, 1] ** 2 - p[:, 2] ** 2 - p[:, 3] ** 2

    p1, p2, p3, p4 = momenta[:, 0], momenta[:, 1], momenta[:, 2], momenta[:, 3]
    s = _sq(p1 + p2)
    t = _sq(p1 - p3)
    u = _sq(p1 - p4)
    m_inv = np.sqrt(np.maximum(s, 0.0))
    ct = (u - t) / (u + t)
    fig, ax = plt.subplots(figsize=(6, 5))
    hb = ax.hexbin(
        m_inv,
        ct,
        C=me2,
        reduce_C_function=np.mean,
        gridsize=50,
        bins="log",
        cmap="viridis",
    )
    ax.set_xlabel(r"$m_{\mathrm{inv}}$ [GeV]")
    ax.set_ylabel(r"$\cos\theta^{*} = (u-t)/(u+t)$")
    fig.colorbar(hb, ax=ax, label=r"mean $|\mathcal{M}|^2$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote plot to {out_path}")


def run(
    process,
    n_events=10000,
    m_inv_range=None,
    costheta_range=None,
    workdir=None,
    dR=0.3,
    plot=False,
    mg_path=None,
    overwrite_alphas=None,
    overwrite_PDGs=None,
):
    _resolve_mg_bin(mg_path)
    n_final = _count_final_state(process)
    slug = _slug(process)
    workdir = workdir or os.path.join(HERE, "output", f"proc_{slug}")
    os.makedirs(workdir, exist_ok=True)
    me_dir = os.path.join(workdir, "madevent")
    sa_dir = os.path.join(workdir, "standalone")

    log = os.path.join(workdir, "madgraph.log")
    open(log, "w").close()
    print(f"[generate] process: {process!r} (n_final={n_final})")
    print(f"[generate] logging MadGraph output to {log}")

    me_script = os.path.join(workdir, "proc_me.mg5")
    sa_script = os.path.join(workdir, "proc_sa.mg5")
    # Always disable subprocess grouping: with grouping enabled the standalone
    # build only registers a single representative diagram per equivalence
    # class, and smatrixhel returns 0 for every event whose specific (pdgs)
    # tuple does not match the representative.
    group = False

    me_ready = os.path.isfile(os.path.join(me_dir, "Cards", "run_card.dat"))
    sa_ready = os.path.isfile(os.path.join(sa_dir, "Cards", "param_card.dat"))
    if not me_ready:
        _write_mg_script(me_script, "madevent", me_dir, process, group=group)
        print("[generate] generating madevent process directory ...")
        _run_mg(me_script, log)
    else:
        print(f"[generate] reusing existing madevent directory {me_dir}")
    if not sa_ready:
        _write_mg_script(sa_script, "standalone", sa_dir, process, group=group)
        print("[generate] generating standalone process directory ...")
        _run_mg(sa_script, log)
    else:
        print(f"[generate] reusing existing standalone directory {sa_dir}")

    # stale Events/run_01 would make MadGraph refuse to launch with the same tag
    run_dir = os.path.join(me_dir, "Events", "run_01")
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)

    if costheta_range is not None or (plot and n_final == 2):
        warnings.warn(
            "The scattering angle is computed as (u - t)/(u + t); this is only "
            "exact if at least one of the two final-state particles is massless. "
            "For two massive final states the value is a Lorentz-invariant proxy, "
            "not the true cos(theta*).",
            stacklevel=2,
        )
    print("[generate] patching dummy_cuts (m_inv/costheta) ...")
    _patch_dummy_cuts(me_dir, m_inv_range, costheta_range)
    print(f"[generate] launching event generation ({n_events} events) ...")
    lhe = _generate_events(me_dir, n_events, dR, log)
    print("[generate] building allmatrix2py.so ...")
    sub_dir = _build_allmatrix2py(sa_dir, log)

    sys.path.insert(0, sub_dir)
    import all_matrix2py  # noqa: E402

    # lha_read.f calls Fortran STOP if ident_card.dat isn't found under ./Cards.
    # Run initialise from sa_dir so the relative lookup succeeds.
    _cwd = os.getcwd()
    os.chdir(sa_dir)
    try:
        all_matrix2py.initialise(os.path.join(sa_dir, "Cards", "param_card.dat"))
    finally:
        os.chdir(_cwd)

    print("[generate] computing amplitudes ...")
    rows = []
    for pdgs, mom, alphas in _parse_lhe(lhe):
        if not _post_cuts(pdgs, mom, m_inv_range, costheta_range):
            continue
        p2 = _invert(mom.tolist())
        if overwrite_alphas is not None:
            alphas = overwrite_alphas
        if overwrite_PDGs is not None:
            pdgs = overwrite_PDGs
        me2 = all_matrix2py.smatrixhel(pdgs, -1, p2, alphas, 0.0, -1)
        rows.append(mom.flatten().tolist() + [alphas, me2])

    n_part = 2 + n_final
    # row layout: [momenta(4*n_part), alphas, me2]
    arr = np.array(rows, dtype=np.float64) if rows else np.empty((0, 4 * n_part + 2))
    momenta = arr[:, : 4 * n_part].reshape(-1, n_part, 4)
    alphas = arr[:, -2]
    me2 = arr[:, -1]
    out = os.path.join(workdir, f"proc_{slug}.npy")
    np.save(out, arr)
    print(f"Wrote {len(arr)} events (after cuts) to {out}  shape={arr.shape}")
    if plot:
        _plot_hexbin(
            momenta,
            me2,
            n_final,
            os.path.join(workdir, f"proc_{slug}_minv_costheta.png"),
        )
    return momenta, alphas, me2


def _parse_range(s):
    if s is None:
        return None
    lo, hi = s.split(",")
    return (float(lo) if lo else None, float(hi) if hi else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--process",
        type=str,
        required=True,
        help="MadGraph process string, e.g. 'q q > z g' or 'g g > t t~'",
    )
    ap.add_argument(
        "--n-events",
        type=int,
        default=10000,
        help="number of unweighted events to generate (default: 10000)",
    )
    ap.add_argument(
        "--mg-path",
        type=str,
        default=None,
        help=f"path to mg5_aMC binary; saved to {MG_PATH_FILE} for reuse",
    )
    ap.add_argument(
        "--m-inv",
        type=str,
        default=None,
        help="invariant-mass cut 'lo,hi' on all final-state particles",
    )
    ap.add_argument(
        "--costheta",
        type=str,
        default=None,
        help="costheta cut 'lo,hi' (only for 2->2 processes)",
    )
    ap.add_argument(
        "--overwrite_alphas",
        type=float,
        default=None,
        help="Overwrite alpha_s in lhe files",
    )
    ap.add_argument(
        "--overwrite_PDGs",
        type=str,
        default=None,
        help="Comma-separated list of PDG IDs to overwrite in the LHE files, e.g. '21,1,-1'",
    )
    ap.add_argument("--workdir", type=str, default=None)
    ap.add_argument(
        "--plot",
        action="store_true",
        help="hexbin plot of m_inv vs costheta colored by mean |M|^2 (2->2 only)",
    )
    args = ap.parse_args()
    run(
        args.process,
        n_events=args.n_events,
        m_inv_range=_parse_range(args.m_inv),
        costheta_range=_parse_range(args.costheta),
        workdir=args.workdir,
        plot=args.plot,
        mg_path=args.mg_path,
        overwrite_alphas=args.overwrite_alphas,
        overwrite_PDGs=(
            [int(x) for x in args.overwrite_PDGs.split(",")]
            if args.overwrite_PDGs
            else None
        ),
    )


if __name__ == "__main__":
    main()
