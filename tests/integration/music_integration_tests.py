#!/usr/bin/env python3
"""High-level integration tests for MUSIC.

Each test generates its own small initial state (no downloads), runs
MUSIChydro in a temporary directory, and checks a physical invariant or the
agreement of two code paths that must give the same answer, instead of
golden numbers.

  conservation_ideal    E, P and N_B are conserved in a 3D ideal TATB run
  conservation_viscous  the same with shear and bulk viscosity
  conservation_outflow  the same while the fluid leaves through the grid
                        faces: the diagnostic adds the time-integrated face
                        flux
  lumpy_ideal           E, Pz and N_B are conserved in a lumpy ideal TATB
                        event (the superluminal-momentum cap in the dilute
                        cells must keep T^{tau eta})
  lumpy_second_order    a lumpy TATB event with second-order shear terms
                        stays stable (phi7 must not diverge where eta -> 0)
  lumpy_viscous_Du      the lumpy event with shear and bulk viscosity
                        conserves E and Pz to O(dtau^2) with partial_tau u
                        from the ideal equations of motion in dilute cells
                        (default), and stays stable with the backward
                        difference everywhere (Du_from_ideal_EoM = 0)
  lumpy_diffusion       the same event with baryon diffusion keeps N_B
                        positive and conserved (kappa must not turn negative,
                        q must stay smooth in time, the J^tau correction
                        must not be skipped)
  freezeout_balance     E, Px, Pz through the Cornelius surface equal the
                        initial content of the hot region (ideal 3D)
  beastmode3            beastMode 3 (grid enlargement, dtau doubling,
                        coarsening) agrees with a large fixed grid
  rerun                 a 3D TATB run that reaches the grid edge is re-run
                        on a larger grid and agrees with a large-grid run; a
                        profile that cannot enlarge the grid stops; input
                        files survive the output clean-up
  rerun_2d              the same for a boost-invariant run
  rerun_source          the same for the TATB source (Initial_profile 112)
  threads               results do not depend on the number of OpenMP threads
  ipglasma_paths        the same T^{mu nu} given as initial condition
                        (profile 9) and as source term (profile 17) give the
                        same T, v_x, v_z fields
  boost_invariance      a 3D run of an eta-independent state matches 2+1D
  bjorken_rk_order      viscous Bjorken flow is second order in dtau (the
                        second RK stage must use its own tau and transport
                        coefficients)
  gubser_face_reconstruction
                        ideal Gubser flow with n_B against the analytic
                        solution with face_reconstruction = 0 and 1; mode 1
                        must be within 1e-3 and closer than mode 0
  tatb_source           the TATB source (112) deposits the same E, Pz and
                        N_B as the TATB initial condition (111), and both
                        match the model totals sqrt(s)/2 (N_A + N_B)
  memory_history        the in-memory hydro history reproduces the
                        evolution file (needs integration_hydro_memory.e)

usage (from anywhere):
  python3 music_integration_tests.py --music ./MUSIChydro \
      [--memory-driver build/src/integration_hydro_memory.e] [tests ...]
"""

import argparse
import glob
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

HBARC = 0.19733


class TestFailure(Exception):
    pass


class TestSkipped(Exception):
    pass


def check(condition, message):
    if not condition:
        raise TestFailure(message)


# ---------------------------------------------------------------------------
# running MUSIC
# ---------------------------------------------------------------------------

class Runner:
    def __init__(self, music, memory_driver, workroot, keep, threads):
        self.music = os.path.abspath(music)
        self.memory_driver = (os.path.abspath(memory_driver)
                              if memory_driver else None)
        self.workroot = workroot
        self.keep = keep
        self.threads = threads

    def workdir(self, name):
        path = os.path.join(self.workroot, name)
        shutil.rmtree(path, ignore_errors=True)
        os.makedirs(path)
        return path

    def run(self, workdir, params, input_name="music_input", threads=None,
            exe=None, expect_success=True, timeout=1800):
        """write the input file, run MUSIC in workdir, return the log"""
        write_input(os.path.join(workdir, input_name), params)
        env = dict(os.environ)
        nthreads = threads if threads is not None else self.threads
        if nthreads is not None:
            env["OMP_NUM_THREADS"] = str(nthreads)
        cmd = [exe or self.music, input_name]
        t0 = time.time()
        proc = subprocess.run(cmd, cwd=workdir, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=timeout)
        log = proc.stdout.decode("utf-8", errors="replace")
        logname = os.path.splitext(input_name)[0] + ".log"
        with open(os.path.join(workdir, logname), "w") as f:
            f.write(log)
        if expect_success:
            check(proc.returncode == 0,
                  "{} exited with {} in {} (see {})".format(
                      os.path.basename(cmd[0]), proc.returncode, workdir,
                      logname))
        return log, proc.returncode, time.time() - t0


def write_input(path, params):
    with open(path, "w") as f:
        for key, value in params.items():
            f.write("{} {}\n".format(key, value))
        f.write("EndOfData\n")


def freeze_out_params(e_fo, every=5):
    return {
        "Do_FreezeOut_Yes_1_No_0": 1,
        "Do_FreezeOut_lowtemp": 0,
        "freeze_out_method": 4,
        "freeze_surface_in_binary": 1,
        "average_surface_over_this_many_time_steps": every,
        "freeze_Ncell_x_step": 1,
        "freeze_Ncell_eta_step": 1,
        "freeze_eps_flag": 0,
        "N_freeze_out": 1,
        "eps_freeze_max": e_fo,
        "eps_freeze_min": e_fo,
        "use_eps_for_freeze_out": 1,
    }


def base_params(**kw):
    params = {
        "echo_level": 1,
        "mode": 2,
        "EOS_to_use": 0,
        "Viscosity_Flag_Yes_1_No_0": 0,
        "Include_Shear_Visc_Yes_1_No_0": 0,
        "Include_Bulk_Visc_Yes_1_No_0": 0,
        "Include_Rhob_Yes_1_No_0": 0,
        "turn_on_baryon_diffusion": 0,
        "output_evolution_data": 0,
        "Do_FreezeOut_Yes_1_No_0": 0,
        "beastMode": 0,
    }
    params.update(kw)
    return params


# ---------------------------------------------------------------------------
# initial states
# ---------------------------------------------------------------------------

def write_TATB(workdir, n, size, sigma=1.2, NA=15., NB=10., b=1.0):
    """Gaussian thickness functions on an n x n grid of side `size` fm"""
    x = np.linspace(-size / 2., size / 2., n)
    X, Y = np.meshgrid(x, x, indexing="ij")

    def gauss(xc, norm):
        return norm / (2 * math.pi * sigma**2) * np.exp(
            -((X - xc)**2 + Y**2) / (2 * sigma**2))

    np.savetxt(os.path.join(workdir, "TA.dat"), gauss(+b / 2., NA), fmt="%.10e")
    np.savetxt(os.path.join(workdir, "TB.dat"), gauss(-b / 2., NB), fmt="%.10e")


def tatb_params(n, size, neta, eta_size, duration, **kw):
    params = base_params(
        Initial_profile=111,
        Initial_TA_Distribution_Filename="TA.dat",
        Initial_TB_Distribution_Filename="TB.dat",
        ecm=19.6,
        Eta_plateau_size=1.0,
        Eta_fall_off=0.7,
        initial_eta_rhob_profile=2,
        eta_rhob_0=1.0,
        eta_rhob_width_1=0.3,
        eta_rhob_width_2=0.5,
        Include_Rhob_Yes_1_No_0=1,
        boost_invariant=0,
        Initial_time_tau_0=1.0,
        Total_evolution_time_tau=duration,
        Delta_Tau=0.02,
        Eta_grid_size=eta_size,
        Grid_size_in_eta=neta,
        X_grid_size_in_fm=size,
        Y_grid_size_in_fm=size,
        Grid_size_in_x=n,
        Grid_size_in_y=n,
    )
    params.update(kw)
    return params


def gaussian_fireball(x, y, eta, e0=20., R=1.5, eta0=0.3, sigma_eta=1.0,
                      x0=0.3, ux0=0.0):
    """e (GeV/fm^3) and u^x on a grid; eta may be None for 2+1D"""
    X, Y = np.meshgrid(x, y, indexing="ij")
    e2d = e0 * np.exp(-((X - x0)**2 + (Y / 1.2)**2) / (2 * R**2))
    ux = ux0 * np.ones_like(X)
    if eta is None:
        return e2d, ux
    env = np.exp(-(np.asarray(eta) - eta0)**2 / (2 * sigma_eta**2))
    return env[:, None, None] * e2d[None, :, :], np.broadcast_to(
        ux, (len(eta),) + X.shape)


def write_ipglasma_file(path, x, y, e, ux, uy, ueta, tau0, eta=None):
    """IP-Glasma text format (profiles 9, 17 and, with eta, 95)

    e in GeV/fm^3; u^eta in 1/fm (the readers multiply by tau0).
    2D: e[ix, iy]; 3D: e[ieta, ix, iy]. Rows: (ieta), ix, iy with iy fastest.
    """
    nx, ny = len(x), len(y)
    dx, dy = x[1] - x[0], y[1] - y[0]
    if eta is None:
        e, ux, uy, ueta = (np.asarray(a)[None] for a in (e, ux, uy, ueta))
        etas, neta, deta = [0.], 1, 0.
    else:
        etas, neta = list(eta), len(eta)
        deta = eta[1] - eta[0] if neta > 1 else 0.
    zeros = " ".join(["0"] * 11)
    with open(path, "w") as f:
        f.write("# dummy 1 etamax= {} xmax= {} ymax= {} deta= {} dx= {} "
                "dy= {}\n".format(neta, nx, ny, deta, dx, dy))
        for ieta in range(neta):
            for ix in range(nx):
                for iy in range(ny):
                    uxv = ux[ieta, ix, iy]
                    uyv = uy[ieta, ix, iy]
                    uev = ueta[ieta, ix, iy]
                    utau = math.sqrt(1. + uxv**2 + uyv**2 + (tau0 * uev)**2)
                    f.write("{:.6g} {:.8g} {:.8g} {:.10e} {:.12g} {:.12g} "
                            "{:.12g} {:.12g} {}\n".format(
                                etas[ieta], x[ix], y[iy], e[ieta, ix, iy],
                                utau, uxv, uyv, uev, zeros))


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------

def read_conservation(workdir):
    path = os.path.join(workdir, "global_conservation_laws.dat")
    check(os.path.exists(path), "global_conservation_laws.dat not written")
    data = np.loadtxt(path)
    check(data.ndim == 2 and len(data) > 2,
          "global_conservation_laws.dat has too few rows")
    return data[:, :6]  # tau, E, Px, Py, Pz, N_B (N_Q, N_S follow)


def read_surface(workdir, nfloat=36):
    files = sorted(glob.glob(os.path.join(workdir, "surface_eps_*.dat")))
    check(files, "no freeze-out surface written")
    chunks = [np.fromfile(f, dtype=np.float32) for f in files]
    data = np.concatenate(chunks).astype(np.float64)
    check(data.size % nfloat == 0, "surface file size is not a multiple of "
          "{} floats".format(nfloat))
    data = data.reshape(-1, nfloat)
    check(len(data) > 0, "empty freeze-out surface")
    check(np.all(np.isfinite(data)), "non-finite numbers in the surface")
    return data


def surface_fluxes(surf):
    """ideal-fluid fluxes through the surface, in GeV (E, Px, Py, Pz) and
    the entropy flux; for a boost-invariant surface they are per unit eta"""
    tau, eta = surf[:, 0], surf[:, 3]
    S = surf[:, 4:8]
    u0, ux, uy = surf[:, 8], surf[:, 9], surf[:, 10]
    ueta = surf[:, 11] / tau  # true u^eta
    e = surf[:, 12]
    T = surf[:, 13]
    w = surf[:, 17] * T  # e + P
    P = w - e
    u = np.stack([u0, ux, uy, ueta], axis=1)
    g = np.stack([np.ones_like(tau), -np.ones_like(tau), -np.ones_like(tau),
                  -1. / tau**2], axis=1)

    def row(unu, nu):
        # T^{mu nu} for fixed nu, as a function of mu
        V = w[:, None] * u * unu[:, None]
        V[:, nu] -= P * g[:, nu]
        return V

    def dot(V):
        return tau * np.sum(S * V, axis=1)

    Ttau, Tx, Ty = row(u0, 0), row(ux, 1), row(uy, 2)
    Teta = row(ueta, 3)
    ch, sh = np.cosh(eta), np.sinh(eta)
    E = np.sum(ch * dot(Ttau) + sh * tau * dot(Teta)) * HBARC
    Pz = np.sum(sh * dot(Ttau) + ch * tau * dot(Teta)) * HBARC
    Px = np.sum(dot(Tx)) * HBARC
    Py = np.sum(dot(Ty)) * HBARC
    s = w / T
    entropy = np.sum(s * dot(u))
    return {"E": E, "Px": Px, "Py": Py, "Pz": Pz, "S": entropy,
            "n": len(surf), "u_perp": np.sum(
                s * dot(u) * np.hypot(ux, uy)) / entropy}


def read_table(workdir, name):
    path = os.path.join(workdir, name)
    check(os.path.exists(path), name + " not written")
    return np.loadtxt(path)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def report(name, value, limit):
    print("    {:<46s} {:10.3e}  (limit {:.1e})".format(name, value, limit))
    check(value <= limit, "{} = {:.3e} exceeds {:.1e}".format(
        name, value, limit))


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

def _conservation(r, name, viscous):
    wd = r.workdir(name)
    write_TATB(wd, 41, 12.0, NA=15., NB=10., b=1.0)
    params = tatb_params(41, 12.0, 49, 12.0, 1.8)
    if viscous:
        params.update({
            "Viscosity_Flag_Yes_1_No_0": 1,
            "Include_Shear_Visc_Yes_1_No_0": 1,
            "Shear_to_S_ratio": 0.12,
            "Include_Bulk_Visc_Yes_1_No_0": 1,
            "Include_second_order_terms": 1,
        })
    r.run(wd, params)
    data = read_conservation(wd)
    tau, E, Px, Py, Pz, NB = data.T
    print("    tau = {:.3f} -> {:.3f} fm/c, E0 = {:.4g} GeV, N_B0 = {:.4g}, "
          "Pz0 = {:.4g} GeV".format(tau[0], tau[-1], E[0], NB[0], Pz[0]))
    check(abs(Pz[0]) > 1e-3 * E[0],
          "the initial state should carry net Pz (TA != TB)")
    # the TATB model puts sqrt(s)/2 per participant into the fluid
    report("|E0 / (sqrt(s)/2 (N_A + N_B)) - 1|", rel(E[0], 9.8 * 25.), 1e-3)
    report("|N_B0 / (N_A + N_B) - 1|", rel(NB[0], 25.), 1e-3)
    # the diagnostic reports the two-level total the scheme conserves; with
    # the deferred viscous correction the viscous run is conserved to
    # O(dtau^2) (2e-6 at dtau = 0.02)
    report("max |E/E0 - 1|", np.max(np.abs(E / E[0] - 1)), 2e-5)
    report("max |N_B/N_B0 - 1|", np.max(np.abs(NB / NB[0] - 1)), 2e-5)
    report("max |Pz - Pz0| / E0", np.max(np.abs(Pz - Pz[0])) / E[0], 2e-5)
    report("max |Px - Px0| / E0", np.max(np.abs(Px - Px[0])) / E[0], 1e-7)
    report("max |Py| / E0", np.max(np.abs(Py)) / E[0], 1e-8)


def test_conservation_ideal(r):
    _conservation(r, "conservation_ideal", viscous=False)


def test_conservation_viscous(r):
    _conservation(r, "conservation_viscous", viscous=True)


def test_conservation_outflow(r):
    # an 8 fm box: matter flows out through the x, y and eta faces
    wd = r.workdir("conservation_outflow")
    write_TATB(wd, 27, 8.0, NA=15., NB=10., b=1.0)
    r.run(wd, tatb_params(27, 8.0, 17, 4.0, 3.0))
    tau, E, Px, Py, Pz, NB = read_conservation(wd).T
    print("    tau = {:.3f} -> {:.3f} fm/c, E0 = {:.4g} GeV".format(
        tau[0], tau[-1], E[0]))
    report("max |E/E0 - 1|", np.max(np.abs(E / E[0] - 1)), 5e-5)
    report("max |N_B/N_B0 - 1|", np.max(np.abs(NB / NB[0] - 1)), 5e-5)
    report("max |Pz - Pz0| / E0", np.max(np.abs(Pz - Pz[0])) / E[0], 1e-5)
    report("max |Px - Px0| / E0", np.max(np.abs(Px - Px[0])) / E[0], 1e-6)


def _write_lumpy_TATB(wd, n, size):
    # 40 nucleons of width 0.5 fm per nucleus
    x = np.linspace(-size / 2., size / 2., n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    rng = np.random.default_rng(7)
    width = 0.5
    for name, xc in (("TA.dat", -1.0), ("TB.dat", 1.0)):
        T = np.zeros_like(X)
        nucleons = rng.normal(0., 2.2, size=(40, 2))
        nucleons[:, 0] += xc
        for px, py in nucleons:
            T += (np.exp(-((X - px)**2 + (Y - py)**2) / (2. * width**2))
                  / (2. * np.pi * width**2))
        np.savetxt(os.path.join(wd, name), T, fmt="%.10e")


def test_lumpy_ideal(r):
    # the dilute cells at the edges of the lumps are pushed past the light
    # cone every step and capped to |M| = 0.999 T^{tau tau}. Scaling
    # T^{tau eta} along with the transverse momentum removed sinh(eta)
    # tau T^{tau eta} from E in each of them: E fell by 3.2e-4 and Pz drifted
    # by 4e-5, independent of dtau. Now 1.6e-5 and 1e-5; the rest comes
    # from the few cells where T^{tau eta} alone is superluminal and falls
    # with deta.
    n, size = 41, 12.0
    wd = r.workdir("lumpy_ideal")
    _write_lumpy_TATB(wd, n, size)
    r.run(wd, tatb_params(n, size, 17, 10.0, 2.0))
    tau, E, Px, Py, Pz, NB = read_conservation(wd).T
    print("    tau = {:.3f} -> {:.3f} fm/c, E0 = {:.4g} GeV".format(
        tau[0], tau[-1], E[0]))
    report("max |E/E0 - 1|", np.max(np.abs(E / E[0] - 1)), 5e-5)
    report("max |Pz - Pz0| / E0", np.max(np.abs(Pz - Pz[0])) / E[0], 2e-5)
    report("max |N_B/N_B0 - 1|", np.max(np.abs(NB / NB[0] - 1)), 2e-5)


def test_lumpy_second_order(r):
    # 40 nucleons of width 0.5 fm per nucleus: near-vacuum cells between the
    # lumps have eta -> 0 while tau_pi is clamped. The phi7 pi pi term used
    # to diverge there, drive pi to the regulator cap (7.6 (e + P)) and
    # blow up the event within ~1 fm/c.
    n, size = 41, 12.0
    wd = r.workdir("lumpy_second_order")
    _write_lumpy_TATB(wd, n, size)
    params = tatb_params(n, size, 17, 10.0, 1.0)
    params.update({
        "Viscosity_Flag_Yes_1_No_0": 1,
        "Include_Shear_Visc_Yes_1_No_0": 1,
        "Shear_to_S_ratio": 0.12,
        "Include_Bulk_Visc_Yes_1_No_0": 1,
        "Include_second_order_terms": 1,
    })
    log, status, _ = r.run(wd, params, expect_success=False)
    eps_max = [float(l.split("eps_max =")[1].split()[0])
               for l in log.splitlines() if "eps_max =" in l]
    check(eps_max, "no eps_max in the log")
    check(status == 0, "the lumpy event did not finish (exit {}, eps_max "
          "reached {:.3g} GeV/fm^3)".format(status, max(eps_max)))
    # the event only expands: the maximum energy density must not grow
    report("max eps_max / initial eps_max - 1",
           max(eps_max) / eps_max[0] - 1., 1e-2)
    # 4e-6 at dtau = 0.02 (2e-4 with Du_from_ideal_EoM = 0, see
    # lumpy_viscous_Du)
    tau, E = read_conservation(wd)[:, :2].T
    report("max |E/E0 - 1|", np.max(np.abs(E / E[0] - 1)), 2e-5)


def test_lumpy_viscous_Du(r):
    # with partial_tau u from a backward time difference the viscous terms
    # are not smooth in time in the dilute, fast cells between the lumps;
    # their deferred conservation correction has to be skipped there
    # (smoothness guard 0.25) and E drifts by 2.8e-4 (O(dtau)). By default
    # partial_tau u comes from the ideal equations of motion in cells below
    # 0.003 GeV/fm^3 (blended up to 0.015), where the guard is relaxed to 4:
    # 5e-6.
    n, size = 41, 12.0
    for flag, limit_E, limit_Pz in ((None, 2e-5, 1e-5), (0, 1e-3, 1e-4)):
        name = "lumpy_viscous_Du" + ("" if flag is None else "_backward")
        wd = r.workdir(name)
        _write_lumpy_TATB(wd, n, size)
        params = tatb_params(n, size, 17, 10.0, 2.0)
        params.update({
            "Viscosity_Flag_Yes_1_No_0": 1,
            "Include_Shear_Visc_Yes_1_No_0": 1,
            "Shear_to_S_ratio": 0.12,
            "Include_Bulk_Visc_Yes_1_No_0": 1,
            "Include_second_order_terms": 1,
        })
        if flag is not None:
            params["Du_from_ideal_EoM"] = flag
        r.run(wd, params, input_name="music_input")
        tau, E, Px, Py, Pz, NB = read_conservation(wd).T
        print("    Du_from_ideal_EoM = {}: tau -> {:.3f} fm/c".format(
            "default" if flag is None else flag, tau[-1]))
        report("max |E/E0 - 1|", np.max(np.abs(E / E[0] - 1)), limit_E)
        report("max |Pz - Pz0| / E0", np.max(np.abs(Pz - Pz[0])) / E[0],
               limit_Pz)


def test_lumpy_diffusion(r):
    # cells between the lumps get n ~ 0.7 fm^-3 with e ~ 0.1 GeV/fm^3 (less
    # than a nucleon mass per baryon); there the kinetic-theory kappa used to
    # turn negative (anti-diffusion) and N_B blew up within ~0.5 fm/c. With
    # muB/T ~ 50 there, D(muB/T) from a backward time difference put
    # time-step-scale noise on q and N_B drifted by 4e-3, growing as
    # dtau -> 0; D(muB/T) now comes from the ideal equations of motion.
    n, size = 41, 12.0
    wd = r.workdir("lumpy_diffusion")
    _write_lumpy_TATB(wd, n, size)
    params = tatb_params(n, size, 17, 10.0, 2.0)
    params.update({
        "Viscosity_Flag_Yes_1_No_0": 1,
        "turn_on_baryon_diffusion": 1,
        "kappa_coefficient": 0.2,
    })
    log, status, _ = r.run(wd, params, expect_success=False)
    check(status == 0, "the lumpy event with diffusion did not finish "
          "(exit {}, see music_input.log)".format(status))
    tau, NB = read_conservation(wd)[:, [0, 5]].T
    print("    tau = {:.3f} -> {:.3f} fm/c, N_B0 = {:.4g}".format(
        tau[0], tau[-1], NB[0]))
    check(np.all(NB > 0.), "N_B went negative")
    # 8e-5 at dtau = 0.02. The J^tau part of the deferred viscous correction
    # has its own smoothness guard (4); with the shared guard (0.25) it was
    # skipped in the cold, baryon-rich cells and N_B drifted by 7e-4
    report("max |N_B/N_B0 - 1|", np.max(np.abs(NB / NB[0] - 1)), 3e-4)


def _fireball_3d(n=41, size=12.0, neta=41, eta_size=8.0):
    x = np.linspace(-size / 2., size / 2., n)
    eta = np.linspace(-eta_size / 2., eta_size / 2., neta)
    e, ux = gaussian_fireball(x, x, eta, ux0=0.05)
    return x, eta, e, ux


def _hot_region_content(x, eta, e, ux, tau0, e_fo):
    """E, Px, Pz (GeV) of the cells above e_fo at tau0 (ideal gas, P=e/3)"""
    dx, deta = x[1] - x[0], eta[1] - eta[0]
    hot = e > e_fo
    u0 = np.sqrt(1. + ux**2)
    w = 4. / 3. * e
    Ttt = w * u0**2 - e / 3.
    Ttx = w * u0 * ux
    ch = np.cosh(eta)[:, None, None]
    sh = np.sinh(eta)[:, None, None]
    vol = tau0 * dx * dx * deta
    return {"E": np.sum((ch * Ttt)[hot]) * vol,
            "Px": np.sum(Ttx[hot]) * vol,
            "Pz": np.sum((sh * Ttt)[hot]) * vol}


def test_freezeout_balance(r):
    wd = r.workdir("freezeout_balance")
    tau0, e_fo = 0.6, 0.3
    x, eta, e, ux = _fireball_3d()
    zero = np.zeros_like(e)
    write_ipglasma_file(os.path.join(wd, "ic.dat"), x, x, e, ux, zero, zero,
                        tau0, eta=eta)
    params = base_params(
        Initial_profile=95, Initial_Distribution_input_filename="ic.dat",
        boost_invariant=0, Initial_time_tau_0=tau0,
        Total_evolution_time_tau=20., Delta_Tau=0.02,
        X_grid_size_in_fm=x[-1] - x[0], Y_grid_size_in_fm=x[-1] - x[0],
        Grid_size_in_x=len(x), Grid_size_in_y=len(x),
        Eta_grid_size=eta[-1] - eta[0], Grid_size_in_eta=len(eta),
        Eta_plateau_size=100., **freeze_out_params(e_fo))
    log, _, _ = r.run(wd, params)
    check("All cells frozen out" in log, "the run did not freeze out fully")
    hot = _hot_region_content(x, eta, e, ux, tau0, e_fo)
    flux = surface_fluxes(read_surface(wd))
    print("    hot region at tau0: E = {:.4g} GeV, Px = {:.4g} GeV, "
          "Pz = {:.4g} GeV; {} surface cells".format(
              hot["E"], hot["Px"], hot["Pz"], flux["n"]))
    report("|E_surface / E_hot - 1|", rel(flux["E"], hot["E"]), 1e-3)
    report("|Px_surface / Px_hot - 1|", rel(flux["Px"], hot["Px"]), 1e-2)
    report("|Pz_surface / Pz_hot - 1|", rel(flux["Pz"], hot["Pz"]), 1e-3)
    report("|Py_surface| / E_hot", abs(flux["Py"]) / hot["E"], 1e-3)


def _fireball_2d_params(x, tau0, e_fo, **kw):
    params = base_params(
        Initial_profile=9, Initial_Distribution_input_filename="ic.dat",
        boost_invariant=1, Initial_time_tau_0=tau0,
        Total_evolution_time_tau=30., Delta_Tau=0.02,
        X_grid_size_in_fm=x[-1] - x[0], Y_grid_size_in_fm=x[-1] - x[0],
        Grid_size_in_x=len(x), Grid_size_in_y=len(x),
        Eta_grid_size=0., Grid_size_in_eta=1, Eta_plateau_size=100.,
        **freeze_out_params(e_fo, every=5))
    params.update(kw)
    return params


def woods_saxon_2d(x, e0=60., R=2.0, a=0.25):
    """flat-top fireball with a sharp edge (e in GeV/fm^3)"""
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X, Y / 1.1)
    return e0 / (1. + np.exp((r - R) / a))


def _write_2d_fireball(wd, n, size, tau0):
    x = np.linspace(-size / 2., size / 2., n)
    e = woods_saxon_2d(x)
    zero = np.zeros_like(e)
    write_ipglasma_file(os.path.join(wd, "ic.dat"), x, x, e, zero, zero,
                        zero, tau0)
    return x


def test_beastmode3(r):
    # dx = 0.19 fm, dtau = 0.02 fm/c: after 192 steps dx/(2 dtau) < 5, so
    # beastMode 3 doubles dtau and coarsens the grid. The fireball fits in
    # the 8.36 fm start grid at tau0 but outgrows it later.
    # beastMode > 0 limits the face states in the primitive variables
    # (face_reconstruction = 1), so the reference does too: this test is
    # about the grid changes, not the face scheme
    tau0, e_fo = 0.6, 0.3
    ref = r.workdir("beastmode3_reference")
    x = _write_2d_fireball(ref, 161, 30.4, tau0)
    log, _, t_ref = r.run(ref, _fireball_2d_params(
        x, tau0, e_fo, face_reconstruction=1))
    check("All cells frozen out" in log, "reference run did not freeze out")
    f_ref = surface_fluxes(read_surface(ref))

    bm = r.workdir("beastmode3")
    x = _write_2d_fireball(bm, 45, 8.36, tau0)
    check(woods_saxon_2d(x)[0].max() < 0.05 and woods_saxon_2d(x)[:, 0].max()
          < 0.05, "the initial state does not fit in the start grid")
    log, _, t_bm = r.run(bm, _fireball_2d_params(
        x, tau0, e_fo, beastMode=3, gridPadding=2.0))
    check("All cells frozen out" in log, "beastMode 3 run did not freeze out")
    nenlarge = log.count("BeastMode: enlarged the grid")
    ncoarsen = log.count(", dx = ")
    print("    beastMode 3: enlarged the grid {} times, doubled dtau and "
          "coarsened {} times".format(nenlarge, ncoarsen))
    check(nenlarge > 0, "beastMode 3 never enlarged the grid")
    check(ncoarsen > 0, "beastMode 3 never coarsened the grid")
    f_bm = surface_fluxes(read_surface(bm))
    print("    reference: {} cells, {:.1f} s; beastMode 3: {} cells, {:.1f} s"
          .format(f_ref["n"], t_ref, f_bm["n"], t_bm))
    report("|S_beast / S_ref - 1| (entropy flux)", rel(f_bm["S"], f_ref["S"]),
           5e-3)
    report("|E_beast / E_ref - 1| (dE/deta flux)", rel(f_bm["E"], f_ref["E"]),
           5e-3)
    report("|<u_perp>_beast / <u_perp>_ref - 1|",
           rel(f_bm["u_perp"], f_ref["u_perp"]), 5e-3)


def test_rerun(r):
    # a 27 x 27 grid (8 fm) is too small; each re-run adds 10 cells per side
    n, size, pad = 27, 8.0, 20
    dx = size / (n - 1)
    common = dict(Total_evolution_time_tau=30.,
                  **freeze_out_params(0.3, every=5))

    wd = r.workdir("rerun")
    write_TATB(wd, n, size, sigma=1.5, NA=30., NB=25., b=0.5)
    with open(os.path.join(wd, "freeze_list.dat"), "w") as f:
        f.write("# e_FO\n")
    with open(os.path.join(wd, "old_output.dat"), "w") as f:
        f.write("0\n")  # an output of an earlier run: must be removed
    p = tatb_params(n, size, 17, 8.0, 30., **common)
    p["freeze_list_filename"] = "freeze_list.dat"
    log, _, t_rr = r.run(wd, p)
    nstart = log.count("initArena")
    print("    hydro started {} times ({:.1f} s)".format(nstart, t_rr))
    check(nstart >= 2, "the run that reached the grid edge was not re-run")
    check("All cells frozen out" in log, "the re-run did not freeze out")
    for name in ("TA.dat", "TB.dat", "freeze_list.dat"):
        check(os.path.exists(os.path.join(wd, name)),
              name + " was removed by the output clean-up")
    check(not os.path.exists(os.path.join(wd, "old_output.dat")),
          "an old output file survived the clean-up")
    f_rr = surface_fluxes(read_surface(wd))

    # the same TA, TB zero-padded on a grid that is large from the start
    ref = r.workdir("rerun_reference")
    for name in ("TA.dat", "TB.dat"):
        T = np.loadtxt(os.path.join(wd, name))
        np.savetxt(os.path.join(ref, name), np.pad(T, pad), fmt="%.10e")
    nref = n + 2 * pad
    p = tatb_params(nref, dx * (nref - 1), 17, 8.0, 30., **common)
    log, _, _ = r.run(ref, p)
    check(log.count("initArena") == 1, "the reference run was re-run")
    f_ref = surface_fluxes(read_surface(ref))
    report("|E_rerun / E_ref - 1| (surface energy)",
           rel(f_rr["E"], f_ref["E"]), 1e-3)
    report("|S_rerun / S_ref - 1| (entropy flux)", rel(f_rr["S"], f_ref["S"]),
           1e-3)

    # a profile whose re-run would use the same grid must stop after one try
    st = r.workdir("rerun_same_grid")
    tau0 = 0.6
    x = np.linspace(-3., 3., 21)
    eta = np.linspace(-2., 2., 11)
    e, _ = gaussian_fireball(x, x, eta, e0=20., R=1.5, eta0=0.)
    zero = np.zeros_like(e)
    write_ipglasma_file(os.path.join(st, "ic.dat"), x, x, e, zero, zero, zero,
                        tau0, eta=eta)
    params = base_params(
        Initial_profile=95, Initial_Distribution_input_filename="ic.dat",
        boost_invariant=0, Initial_time_tau_0=tau0,
        Total_evolution_time_tau=10., Delta_Tau=0.02,
        X_grid_size_in_fm=6., Y_grid_size_in_fm=6., Grid_size_in_x=21,
        Grid_size_in_y=21, Eta_grid_size=4., Grid_size_in_eta=11,
        Eta_plateau_size=100., **freeze_out_params(0.3))
    log, status, _ = r.run(st, params, expect_success=False)
    check(status != 0, "a run truncated at the grid edge exited with 0")
    check("The grid is too small" in log,
          "profile 95 on a 6 fm grid did not reach the grid edge")
    check(log.count("initArena") == 1, "profile 95 was re-run on the same grid")
    check("would use the same grid. Not re-running" in log,
          "missing the message about not re-running")


def test_rerun_2d(r):
    # the boost-invariant finder must also stop at the grid edge and re-run
    n, size, pad = 27, 8.0, 20
    dx = size / (n - 1)
    common = dict(boost_invariant=1, **freeze_out_params(0.3, every=5))
    wd = r.workdir("rerun_2d")
    write_TATB(wd, n, size, sigma=1.5, NA=30., NB=25., b=0.5)
    log, _, _ = r.run(wd, tatb_params(n, size, 1, 0., 30., **common))
    nstart = log.count("initArena")
    print("    hydro started {} times".format(nstart))
    check(nstart >= 2, "the 2D run that reached the grid edge was not re-run")
    f_rr = surface_fluxes(read_surface(wd))

    ref = r.workdir("rerun_2d_reference")
    for name in ("TA.dat", "TB.dat"):
        T = np.loadtxt(os.path.join(wd, name))
        np.savetxt(os.path.join(ref, name), np.pad(T, pad), fmt="%.10e")
    nref = n + 2 * pad
    log, _, _ = r.run(ref, tatb_params(nref, dx * (nref - 1), 1, 0., 30.,
                                       **common))
    check(log.count("initArena") == 1, "the reference run was re-run")
    f_ref = surface_fluxes(read_surface(ref))
    report("|E_rerun / E_ref - 1| (surface energy)",
           rel(f_rr["E"], f_ref["E"]), 1e-3)
    report("|S_rerun / S_ref - 1| (entropy flux)", rel(f_rr["S"], f_ref["S"]),
           1e-3)


def test_rerun_source(r):
    # the source keeps TA, TB on the original grid; the re-run only pads
    # the hydro grid, so it must match a run on the padded TA, TB.
    # 23 cells (6.8 fm): the 112 fireball reaches the edge
    dx = 8.0 / 26
    n, pad = 23, 20
    size = dx * (n - 1)
    common = dict(Initial_profile=112, **freeze_out_params(0.3, every=5))
    wd = r.workdir("rerun_source")
    write_TATB(wd, n, size, sigma=1.5, NA=30., NB=25., b=0.5)
    log, _, _ = r.run(wd, tatb_params(n, size, 17, 8.0, 30., **common))
    nstart = log.count("initArena")
    print("    hydro started {} times".format(nstart))
    check(nstart >= 2, "the 112 run that reached the grid edge was not "
          "re-run")
    check("All cells frozen out" in log, "the re-run did not freeze out")
    f_rr = surface_fluxes(read_surface(wd))

    ref = r.workdir("rerun_source_reference")
    for name in ("TA.dat", "TB.dat"):
        T = np.loadtxt(os.path.join(wd, name))
        np.savetxt(os.path.join(ref, name), np.pad(T, pad), fmt="%.10e")
    nref = n + 2 * pad
    log, _, _ = r.run(ref, tatb_params(nref, dx * (nref - 1), 17, 8.0, 30.,
                                       **common))
    check(log.count("initArena") == 1, "the reference run was re-run")
    f_ref = surface_fluxes(read_surface(ref))
    report("|E_rerun / E_ref - 1| (surface energy)",
           rel(f_rr["E"], f_ref["E"]), 1e-3)
    report("|S_rerun / S_ref - 1| (entropy flux)", rel(f_rr["S"], f_ref["S"]),
           1e-3)


def test_threads(r):
    nthreads = 4
    out = {}
    for n in (1, nthreads):
        wd = r.workdir("threads_{}".format(n))
        write_TATB(wd, 41, 12.0, NA=15., NB=10., b=1.0)
        params = tatb_params(41, 12.0, 25, 6.0, 2.6,
                             **freeze_out_params(0.3, every=5))
        params.update(Viscosity_Flag_Yes_1_No_0=1,
                      Include_Shear_Visc_Yes_1_No_0=1, Shear_to_S_ratio=0.12)
        log, _, _ = r.run(wd, params, threads=n)
        if n > 1 and "using 1 threads" in log:
            raise TestSkipped("MUSIChydro was built without OpenMP")
        surf = read_surface(wd)
        order = np.lexsort(surf.T[::-1])
        out[n] = (read_conservation(wd), surf[order])
    c1, s1 = out[1]
    cn, sn = out[nthreads]
    check(c1.shape == cn.shape, "different number of time steps")
    check(s1.shape == sn.shape, "different number of surface cells: {} vs {}"
          .format(len(s1), len(sn)))
    report("max |dE|/E0 between 1 and {} threads".format(nthreads),
           np.max(np.abs(c1[:, 1] - cn[:, 1])) / c1[0, 1], 1e-12)
    report("max |dN_B|/N_B0", np.max(np.abs(c1[:, 5] - cn[:, 5])) / c1[0, 5],
           1e-12)
    report("max |d(surface)| (float32, all columns)",
           np.max(np.abs(s1 - sn)), 0.)


def _ipglasma_state(n=41, size=12.0):
    x = np.linspace(-size / 2., size / 2., n)
    e, _ = gaussian_fireball(x, x, None, e0=30., R=1.5)
    X, Y = np.meshgrid(x, x, indexing="ij")
    ux = 0.1 * X / 3.
    uy = 0.1 * Y / 3.
    ueta = 0.15 * np.ones_like(X)  # 1/fm, a uniform longitudinal boost
    return x, e, ux, uy, ueta


def _read_evolution(workdir, nx, ny, neta):
    data = np.loadtxt(os.path.join(workdir, "evolution_xyeta.dat"))
    check(len(data) % (nx * ny * neta) == 0, "incomplete evolution file")
    # rows: ieta, iy, ix with ix fastest; columns T, mu_B, v_x, v_y, v_z
    return data.reshape(-1, neta, ny, nx, 5)


def test_ipglasma_paths(r):
    tau0, dtau = 0.6, 0.01
    x, e, ux, uy, ueta = _ipglasma_state()
    eta = np.linspace(-0.4, 0.4, 5)
    fields = {}
    for profile in (9, 17):
        wd = r.workdir("ipglasma_profile_{}".format(profile))
        write_ipglasma_file(os.path.join(wd, "ic.dat"), x, x, e, ux, uy, ueta,
                            tau0)
        params = base_params(
            Initial_profile=profile,
            Initial_Distribution_input_filename="ic.dat",
            boost_invariant=0, Initial_time_tau_0=tau0,
            Total_evolution_time_tau=0.15, Delta_Tau=dtau,
            X_grid_size_in_fm=x[-1] - x[0], Y_grid_size_in_fm=x[-1] - x[0],
            Grid_size_in_x=len(x), Grid_size_in_y=len(x),
            Eta_grid_size=eta[-1] - eta[0], Grid_size_in_eta=len(eta),
            Eta_plateau_size=100., Eta_fall_off=1.0,
            output_evolution_data=1, outputBinaryEvolution=0,
            output_evolution_every_N_timesteps=1,
            output_evolution_every_N_x=1, output_evolution_every_N_y=1,
            output_evolution_every_N_eta=1)
        r.run(wd, params)
        fields[profile] = _read_evolution(wd, len(x), len(x), len(eta))
    # profile 17 starts one step earlier and deposits the source over the
    # steps tau0 - dtau -> tau0 -> tau0 + dtau; compare from tau0 + 2 dtau
    n = min(len(fields[9]) - 2, len(fields[17]) - 3)
    check(n >= 5, "too few common time slices")
    f9 = fields[9][2:2 + n]
    f17 = fields[17][3:3 + n]
    hot = f9[..., 0] > 0.1 * f9[..., 0].max()
    names = ("T", "mu_B", "v_x", "v_y", "v_z")
    for col in (0, 2, 4):
        scale = np.max(np.abs(f9[..., col][hot]))
        diff = np.max(np.abs(f17[..., col] - f9[..., col])[hot]) / scale
        report("max |{0}_17 - {0}_9| / max|{0}_9| (hot cells)".format(
            names[col]), diff, 1e-3)


def test_boost_invariance(r):
    tau0 = 0.6
    x, e, ux, uy, _ = _ipglasma_state()
    zero = np.zeros_like(e)
    tables = {}
    for label, neta in (("2d", 1), ("3d", 5)):
        wd = r.workdir("boost_invariance_" + label)
        write_ipglasma_file(os.path.join(wd, "ic.dat"), x, x, e, ux, uy, zero,
                            tau0)
        params = base_params(
            Initial_profile=9, Initial_Distribution_input_filename="ic.dat",
            boost_invariant=1 if neta == 1 else 0, Initial_time_tau_0=tau0,
            Total_evolution_time_tau=3.0, Delta_Tau=0.01,
            X_grid_size_in_fm=x[-1] - x[0], Y_grid_size_in_fm=x[-1] - x[0],
            Grid_size_in_x=len(x), Grid_size_in_y=len(x),
            Eta_grid_size=0.05 * (neta - 1), Grid_size_in_eta=neta,
            Eta_plateau_size=100., Eta_fall_off=1.0)
        r.run(wd, params)
        tables[label] = read_table(wd, "momentum_anisotropy_eta_-0.5_0.5.dat")
    a, b = tables["2d"], tables["3d"]
    n = min(len(a), len(b))
    check(n > 100, "too few time steps")
    a, b = a[:n], b[:n]
    report("max |tau_2d - tau_3d|", np.max(np.abs(a[:, 0] - b[:, 0])), 1e-9)
    scale = np.max(np.abs(a[:, 1:]), axis=0)
    cols = scale > 1e-8
    report("max |anisotropy_3d - anisotropy_2d| / max",
           np.max(np.abs(b[:, 1:][:, cols] - a[:, 1:][:, cols])
                  / scale[cols]), 1e-3)


def test_bjorken_rk_order(r):
    # uniform, boost-invariant fluid with shear and bulk viscosity: theta =
    # u^tau / tau exactly and nothing depends on the spatial resolution, so
    # the change of T(tau) with dtau measures the time-integration error.
    # The second RK stage used tau_n for theta and the transport
    # coefficients of tau_n while its fields are at tau_{n+1}: first order,
    # T changed by 6e-4 between dtau = 0.02 and 0.01 (rounding to the 7
    # digits of the output after the fix)
    n, size, tau0, tau_end = 21, 6.0, 0.4, 1.4
    x = np.linspace(-size / 2., size / 2., n)
    e = 30. * np.ones((n, n))
    zero = np.zeros_like(e)
    T_center = {}
    for dtau in (0.02, 0.01):
        wd = r.workdir("bjorken_rk_order_{}".format(dtau))
        write_ipglasma_file(os.path.join(wd, "ic.dat"), x, x, e, zero, zero,
                            zero, tau0)
        params = base_params(
            Initial_profile=9, Initial_Distribution_input_filename="ic.dat",
            boost_invariant=1, Initial_time_tau_0=tau0,
            Total_evolution_time_tau=tau_end - tau0, Delta_Tau=dtau,
            X_grid_size_in_fm=size, Y_grid_size_in_fm=size,
            Grid_size_in_x=n, Grid_size_in_y=n,
            Eta_grid_size=0., Grid_size_in_eta=1, Eta_plateau_size=100.,
            Viscosity_Flag_Yes_1_No_0=1, Include_Shear_Visc_Yes_1_No_0=1,
            Shear_to_S_ratio=0.2, Include_Bulk_Visc_Yes_1_No_0=1,
            Include_second_order_terms=1,
            output_evolution_data=1, outputBinaryEvolution=0,
            output_evolution_every_N_timesteps=int(round(0.1 / dtau)),
            output_evolution_every_N_x=1, output_evolution_every_N_y=1,
            output_evolution_every_N_eta=1)
        r.run(wd, params)
        fields = _read_evolution(wd, n, n, 1)
        check(len(fields) == 11, "expected 11 time slices")
        # the edges are 3 fm away: the center is still uniform Bjorken flow
        T_center[dtau] = fields[-1, 0, n // 2, n // 2, 0]
    report("|T(dtau = 0.02) / T(dtau = 0.01) - 1| at tau = 1.4",
           abs(T_center[0.02] / T_center[0.01] - 1.), 2e-6)


def test_tatb_source(r):
    results = {}
    NA, NB = 15., 10.
    for profile in (111, 112):
        wd = r.workdir("tatb_profile_{}".format(profile))
        # a large impact parameter puts cells with TA >> TB in the grid
        write_TATB(wd, 41, 12.0, NA=NA, NB=NB, b=3.0)
        params = tatb_params(41, 12.0, 49, 12.0, 0.5, Initial_profile=profile)
        r.run(wd, params)
        results[profile] = read_conservation(wd)
    c1, c2 = results[111], results[112]
    check(np.all(np.isfinite(c2)), "non-finite numbers with profile 112")
    # after the deposit the source run must hold the same totals
    i1 = 1
    i2 = int(np.argmin(np.abs(c2[:, 0] - c1[i1, 0])))
    check(abs(c1[i1, 0] - c2[i2, 0]) < 1e-6, "no common time step")
    print("    tau = {:.3f}: E = {:.5g} / {:.5g} GeV, Pz = {:.5g} / {:.5g} GeV,"
          " N_B = {:.5g} / {:.5g} (111 / 112)".format(
              c1[i1, 0], c1[i1, 1], c2[i2, 1], c1[i1, 4], c2[i2, 4],
              c1[i1, 5], c2[i2, 5]))
    ecm = 19.6
    report("|E_111 / (sqrt(s)/2 (N_A + N_B)) - 1|",
           rel(c1[0, 1], ecm / 2. * (NA + NB)), 1e-3)
    report("|E_112 / E_111 - 1|", rel(c2[i2, 1], c1[i1, 1]), 1e-3)
    report("|N_B,112 / N_B,111 - 1|", rel(c2[i2, 5], c1[i1, 5]), 1e-3)
    report("|Pz_112 / Pz_111 - 1|", rel(c2[i2, 4], c1[i1, 4]), 1e-3)


def test_memory_history(r):
    if not r.memory_driver:
        raise TestSkipped("no --memory-driver given")
    wd = r.workdir("memory_history")
    write_TATB(wd, 31, 12.0, NA=15., NB=10., b=1.0)
    params = tatb_params(31, 12.0, 17, 8.0, 1.8)
    params.update(store_hydro_info_in_memory=1, output_evolution_data=1,
                  outputBinaryEvolution=0,
                  output_evolution_every_N_timesteps=5,
                  output_evolution_every_N_x=2, output_evolution_every_N_y=2,
                  output_evolution_every_N_eta=2)
    log, status, _ = r.run(wd, params, exe=r.memory_driver,
                           expect_success=False)
    lines = [l for l in log.splitlines() if l.startswith("compared ")]
    check(lines, "the driver printed no summary (exit {})".format(status))
    print("    " + lines[-1])
    check(status == 0, "the in-memory history differs from the evolution "
          "file")


def _write_gubser_ideal(path, tau, n=201, size=10.0, q=1.0, e0=1.0,
                        nb0=0.5):
    """the ideal Gubser solution with net baryon density, in the format of
    tests/Gubser_flow/Gubser_solution_ideal.py: x, y, e [1/fm^4],
    n_B [1/fm^3], u^x, u^y with y running fastest"""
    x = np.linspace(-size / 2., size / 2., n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    e, nb, ux, uy = _gubser_ideal(X, Y, tau, q, e0, nb0)
    np.savetxt(path, np.column_stack([
        X.ravel(), Y.ravel(), e.ravel(), nb.ravel(), ux.ravel(), uy.ravel()]),
        fmt="%.8e")


def _gubser_ideal(X, Y, tau, q=1.0, e0=1.0, nb0=0.5):
    r = np.sqrt(X**2 + Y**2)
    den = 1. + 2. * q**2 * (tau**2 + r**2) + q**4 * (tau**2 - r**2)**2
    e = e0 / tau**(4. / 3.) * (2. * q)**(8. / 3.) / den**(4. / 3.)
    nb = nb0 / tau**3 * 4. * (q * tau)**2 / den
    ur = np.sinh(np.arctanh(2. * q**2 * tau * r
                            / (1. + q**2 * tau**2 + q**2 * r**2)))
    rs = np.where(r > 0., r, 1.)
    return e, nb, ur * X / rs, ur * Y / rs


def test_gubser_face_reconstruction(r):
    # ideal Gubser flow with n_B (conformal EoS) against the analytic
    # solution, with the KT face states limited in T^{tau mu}
    # (face_reconstruction = 0) and in e, n, u^i (1). Relative L1 errors
    # on y = 0 at tau = 1.5 fm/c: 1.9e-3 (e), 1.4e-3 (n_B), 1.3e-3 (u^x)
    # for mode 0 and 5.2e-4, 2.6e-4, 6.5e-4 for mode 1.
    hbarc = 0.19733
    tau_end = 1.5
    errors = {}
    for mode in (0, 1):
        wd = r.workdir("gubser_face_reconstruction_{}".format(mode))
        gdir = os.path.join(wd, "tests", "Gubser_flow")
        os.makedirs(gdir)
        for tau in (0.98, 1.0):
            _write_gubser_ideal(
                os.path.join(gdir, "y=0_tau={:4.2f}_ideal.dat".format(tau)),
                tau)
        r.run(wd, base_params(
            Initial_profile=0, boost_invariant=1, Grid_size_in_eta=1,
            Eta_grid_size=10.0, Initial_time_tau_0=1.0,
            Total_evolution_time_tau=tau_end - 1.0, Delta_Tau=0.02,
            UseCFL_condition=0, X_grid_size_in_fm=10.0,
            Y_grid_size_in_fm=10.0, Grid_size_in_x=201, Grid_size_in_y=201,
            Include_Rhob_Yes_1_No_0=1, face_reconstruction=mode))
        d = np.loadtxt(os.path.join(
            wd, "Gubser_flow_check_tau_{:g}.dat".format(tau_end)))
        d = d[np.abs(d[:, 1]) < 1e-6]
        e, nb, ux, _ = _gubser_ideal(d[:, 0], d[:, 1], tau_end)
        l1 = lambda num, ana: np.sum(np.abs(num - ana)) / np.sum(np.abs(ana))
        errors[mode] = (l1(d[:, 2] / hbarc, e), l1(d[:, 3], nb),
                        l1(d[:, 5], ux))
    for i, name in enumerate(("e", "n_B", "u^x")):
        print("    L1 error of {:4s} mode 0: {:.2e}, mode 1: {:.2e}".format(
            name, errors[0][i], errors[1][i]))
        report("mode 1 L1 error of " + name, errors[1][i], 1e-3)
        check(errors[1][i] < errors[0][i],
              "face_reconstruction = 1 is less accurate than 0 in " + name)


def test_gpu_backend(r):
    # the device kernels on the CPU (gpu_backend = 1, double) and, in a
    # build with metal=ON, on the GPU (2, float) against the CPU code of
    # Advance on the lumpy 3D event with n_B: conformal EoS, hotQCD (9) and
    # NEOS (14) and the 4D EoS (20) when their tables are installed.
    # Backend 1: with clang identical bit for bit (0 for all three EoSs).
    # g++ contracts multiply-adds across statements, differently in the two
    # code paths: 1e-16 for EoS 0 and 9; with NEOS the 1e-19 noise grows to
    # 4e-8 by tau = 2.3 fm/c once a capped or reverted cell flips. The bit
    # for bit check runs in the unit test (unittest_bugfix).
    # Metal (float, stored tau T^{tau mu}): ~1e-6 of E.
    # Then the Gubser test (face_reconstruction 1) on the GPU must be as
    # close to the analytic solution as the CPU code.
    n, size = 41, 12.0
    eos_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "..", "EOS")
    limits = {1: 1e-6, 2: 1e-4}
    have = {1: True, 2: True}
    for eos, table in [(0, None), (9, "hotQCD"), (14, "neos_bqs"),
                       (20, "neos4D")]:
        if table and not os.path.isdir(os.path.join(eos_dir, table)):
            print("    EoS {}: no tables in {}, skipped".format(eos, eos_dir))
            continue
        res = {}
        for backend in (0, 1, 2):
            if not have.get(backend, True):
                continue
            wd = r.workdir("gpu_backend_eos{}_{}".format(eos, backend))
            _write_lumpy_TATB(wd, n, size)
            if table:
                os.symlink(eos_dir, os.path.join(wd, "EOS"))
            log, _, t = r.run(wd, tatb_params(
                n, size, 17, 10.0, 1.0, EOS_to_use=eos,
                gpu_backend=backend))
            if backend > 0 and "built without" in log:
                print("    gpu_backend = {}: not in this build, skipped"
                      .format(backend))
                have[backend] = False
                continue
            if backend > 0:
                check("hydro on the device" in log,
                      "gpu_backend = {} fell back to the CPU code (EoS {})"
                      .format(backend, eos))
            res[backend] = (read_conservation(wd), t)
        a = res[0][0]
        for backend in sorted(res):
            print("    EoS {} gpu_backend {}: {:.1f} s".format(
                eos, backend, res[backend][1]))
            if backend == 0:
                continue
            b = res[backend][0]
            diff = np.max(np.abs(a[:, 1:] - b[:, 1:]) / np.abs(a[:, 1:2]))
            report("EoS {} backend {} max |conservation diff| / E".format(
                eos, backend), diff, limits[backend])
    if not have[2]:
        return
    hbarc = 0.19733
    tau_end = 1.5
    errors = {}
    for backend in (0, 2):
        wd = r.workdir("gpu_backend_gubser_{}".format(backend))
        gdir = os.path.join(wd, "tests", "Gubser_flow")
        os.makedirs(gdir)
        for tau in (0.98, 1.0):
            _write_gubser_ideal(
                os.path.join(gdir, "y=0_tau={:4.2f}_ideal.dat".format(tau)),
                tau)
        r.run(wd, base_params(
            Initial_profile=0, boost_invariant=1, Grid_size_in_eta=1,
            Eta_grid_size=10.0, Initial_time_tau_0=1.0,
            Total_evolution_time_tau=tau_end - 1.0, Delta_Tau=0.02,
            UseCFL_condition=0, X_grid_size_in_fm=10.0,
            Y_grid_size_in_fm=10.0, Grid_size_in_x=201, Grid_size_in_y=201,
            Include_Rhob_Yes_1_No_0=1, face_reconstruction=1,
            gpu_backend=backend))
        d = np.loadtxt(os.path.join(
            wd, "Gubser_flow_check_tau_{:g}.dat".format(tau_end)))
        d = d[np.abs(d[:, 1]) < 1e-6]
        e, nb, ux, _ = _gubser_ideal(d[:, 0], d[:, 1], tau_end)
        l1 = lambda num, ana: np.sum(np.abs(num - ana)) / np.sum(np.abs(ana))
        errors[backend] = (l1(d[:, 2] / hbarc, e), l1(d[:, 3], nb),
                           l1(d[:, 5], ux))
    for i, name in enumerate(("e", "n_B", "u^x")):
        print("    Gubser L1 error of {:4s} CPU code: {:.3e}, Metal: {:.3e}"
              .format(name, errors[0][i], errors[2][i]))
        report("Gubser Metal / CPU L1 error of " + name + " - 1",
               abs(errors[2][i] / errors[0][i] - 1.), 1e-2)

def test_gpu_backend_viscous(r):
    # the viscous update on the device: the lumpy 3D event with n_B, shear
    # and bulk with the second-order terms (as lumpy_second_order), for the
    # conformal EoS and hotQCD (9). Backend 1 against the CPU code as in
    # gpu_backend (bit for bit with clang, g++ contracts differently);
    # Metal in float: 8e-7 of E, and W / (e + P) and Pi in the cells above
    # 0.1 GeV/fm^3 to 2e-6 (L1) at tau = 2 fm/c. Also the 4D EoS (20), whose
    # tables are float on the host too.
    n, size = 41, 12.0
    ncell = n * n * 17

    def fields(wd):
        e = np.loadtxt(os.path.join(wd, "evolution_xyeta.dat"))[-ncell:, 0]
        W = np.loadtxt(os.path.join(
            wd, "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat"))[-ncell:]
        Pi = np.loadtxt(os.path.join(
            wd, "evolution_bulk_pressure_xyeta.dat"))[-ncell:, 1]
        return e, W, Pi

    eos_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "..", "EOS")
    limits = {1: 1e-6, 2: 1e-4}
    have = {1: True, 2: True}
    for eos, table in [(0, None), (9, "hotQCD"), (20, "neos4D")]:
        if table and not os.path.isdir(os.path.join(eos_dir, table)):
            print("    EoS {}: no tables in {}, skipped".format(eos, eos_dir))
            continue
        res = {}
        for backend in (0, 1, 2):
            if not have.get(backend, True):
                continue
            wd = r.workdir("gpu_viscous_eos{}_{}".format(eos, backend))
            _write_lumpy_TATB(wd, n, size)
            if table:
                os.symlink(eos_dir, os.path.join(wd, "EOS"))
            params = tatb_params(n, size, 17, 10.0, 1.0, EOS_to_use=eos,
                                 gpu_backend=backend)
            params.update({
                "Viscosity_Flag_Yes_1_No_0": 1,
                "Include_Shear_Visc_Yes_1_No_0": 1,
                "Shear_to_S_ratio": 0.12,
                "Include_Bulk_Visc_Yes_1_No_0": 1,
                "Include_second_order_terms": 1,
                "output_evolution_data": 1,
                "outputBinaryEvolution": 0,
                "output_evolution_every_N_timesteps": 10,
            })
            log, _, t = r.run(wd, params)
            if backend > 0 and "built without" in log:
                print("    gpu_backend = {}: not in this build, skipped"
                      .format(backend))
                have[backend] = False
                continue
            if backend > 0:
                check("hydro on the device" in log and "viscous" in log,
                      "gpu_backend = {} fell back to the CPU code (EoS {})"
                      .format(backend, eos))
            res[backend] = (read_conservation(wd), t, fields(wd))
        a = res[0][0]
        e0, W0, Pi0 = res[0][2]
        hot = e0 > 0.1  # GeV/fm^3
        for backend in sorted(res):
            c = res[backend][0]
            print("    EoS {} gpu_backend {}: {:.1f} s, max |E/E0 - 1| "
                  "{:.2e}".format(eos, backend, res[backend][1],
                                  np.max(np.abs(c[:, 1] / c[0, 1] - 1))))
            if backend == 0:
                continue
            check(len(c) == len(a), "backend {} wrote {} rows, the CPU code "
                  "{}".format(backend, len(c), len(a)))
            diff = np.max(np.abs(a[:, 1:] - c[:, 1:]) / np.abs(a[:, 1:2]))
            report("EoS {} backend {} max |conservation diff| / E".format(
                eos, backend), diff, limits[backend])
            report("EoS {} backend {} max |E/E0 - 1|".format(eos, backend),
                   np.max(np.abs(c[:, 1] / c[0, 1] - 1)), 2e-5)
            _, W, Pi = res[backend][2]
            report("EoS {} backend {} L1 |W/(e+P) diff| (hot cells)".format(
                eos, backend), np.sum(np.abs(W[hot] - W0[hot]))
                / np.sum(np.abs(W0[hot])), {1: 1e-5, 2: 1e-4}[backend])
            report("EoS {} backend {} L1 |Pi diff| (hot cells)".format(
                eos, backend), np.sum(np.abs(Pi[hot] - Pi0[hot]))
                / np.sum(np.abs(Pi0[hot])), {1: 1e-5, 2: 1e-4}[backend])


def _write_strings_event(path, ncopy=40, seed=11):
    """a strings event (Initial_profile 13) from the three strings of
    tests/test_source_terms: copies moved in the transverse plane and
    boosted in eta (positions and rapidities together, so each stays a
    valid string), with collision times tau_0 in [0, 0.4] fm/c so that the
    string ends deposit over the whole run. The first string carries its
    baryon number at eta_s = -0.4 (fraction 1 at the left end) so that n_B
    is deposited at mid-rapidity during the run."""
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                       "test_source_terms", "strings_event_0.dat")
    with open(src) as f:
        head = [f.readline(), f.readline()]
    rows = np.loadtxt(src, ndmin=2)
    rows[0, 23] = 1.  # baryon_fraction_l
    rng = np.random.default_rng(seed)
    # columns: 5-10 transverse positions, eta_s/rapidity columns below
    xcols, ycols = (5, 7, 9), (6, 8, 10)
    etacols = (4, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 27, 28)
    out = []
    for _ in range(ncopy):
        r = rows[rng.integers(len(rows))].copy()
        dx, dy = rng.normal(0., 1.2, size=2)
        r[list(xcols)] += dx
        r[list(ycols)] += dy
        # a baryon or charge eta of 0 marks no baryon or charge at that end
        keep = [c for c in etacols if c in (4, 11, 12, 13, 14, 17, 18)
                or r[c] != 0.]
        r[keep] += rng.uniform(-0.8, 0.8)
        r[3] = rng.uniform(0., 0.4)
        out.append(r)
    with open(path, "w") as f:
        f.writelines(head)
        for r in out:
            f.write("  ".join("{:.8e}".format(v) for v in r) + "\n")


def test_gpu_backend_sources(r):
    # the source terms with the update on the device: the host evaluates
    # them for each RK stage from the u of the device slice (the predictor
    # in the second stage), the device adds them to qi. The lumpy TATB
    # source (112, ideal and viscous) and a strings event (13, viscous with
    # n_B; its n_B source depends on u). Backend 1 against the CPU code:
    # bit for bit with clang (unit test), g++ contracts differently. Metal
    # in float. The difference of E, P and N_B is taken relative to the
    # largest E of the run: the strings run starts from an empty grid.
    n, size = 41, 12.0
    have = {1: True, 2: True}
    cases = []
    visc = {"Viscosity_Flag_Yes_1_No_0": 1, "Include_Shear_Visc_Yes_1_No_0": 1,
            "Shear_to_S_ratio": 0.12, "Include_Bulk_Visc_Yes_1_No_0": 1,
            "Include_second_order_terms": 1}
    for viscous in (0, 1):
        def tatb(wd, backend, viscous=viscous):
            _write_lumpy_TATB(wd, n, size)
            p = tatb_params(n, size, 17, 10.0, 1.0, Initial_profile=112,
                            gpu_backend=backend)
            if viscous:
                p.update(visc)
            return p
        cases.append(("TATB source" + (" viscous" if viscous else ""), tatb))

    def strings(wd, backend):
        _write_strings_event(os.path.join(wd, "strings.dat"))
        p = base_params(
            Initial_profile=13,
            Initial_Distribution_input_filename="strings.dat",
            string_source_sigma_x=0.5, string_source_sigma_eta=0.5,
            stringPreEqFlowFactor=1.0, Include_Rhob_Yes_1_No_0=1,
            boost_invariant=0, Initial_time_tau_0=0.4,
            Total_evolution_time_tau=2.0, Delta_Tau=0.02,
            Eta_grid_size=12.0, Grid_size_in_eta=49,
            X_grid_size_in_fm=size, Y_grid_size_in_fm=size,
            Grid_size_in_x=n, Grid_size_in_y=n, gpu_backend=backend)
        p.update(visc)
        p["Include_Bulk_Visc_Yes_1_No_0"] = 0
        return p
    cases.append(("strings viscous", strings))
    limits = {1: 1e-6, 2: 1e-4}
    for name, make in cases:
        res = {}
        for backend in (0, 1, 2):
            if backend and not have[backend]:
                continue
            wd = r.workdir("gpu_sources_{}_{}".format(
                name.replace(" ", "_"), backend))
            log, _, t = r.run(wd, make(wd, backend))
            if backend > 0 and "built without" in log:
                print("    gpu_backend = {}: not in this build, skipped"
                      .format(backend))
                have[backend] = False
                continue
            if backend > 0:
                check("source terms from the host" in log,
                      "gpu_backend = {} fell back to the CPU code ({})"
                      .format(backend, name))
            res[backend] = (read_conservation(wd), t)
        a = res[0][0]
        for backend in sorted(res):
            c = res[backend][0]
            print("    {} gpu_backend {}: {:.1f} s, E {:.4g} -> {:.4g} GeV, "
                  "N_B {:.4g}".format(name, backend, res[backend][1],
                                      c[0, 1], c[-1, 1], c[-1, 5]))
            if backend == 0:
                continue
            check(len(c) == len(a), "backend {} wrote {} rows, the CPU code "
                  "{}".format(backend, len(c), len(a)))
            diff = np.max(np.abs(a[:, 1:] - c[:, 1:])) / np.max(np.abs(a[:, 1]))
            report("{} backend {} max |conservation diff| / E_max".format(
                name, backend), diff, limits[backend])
            if np.max(a[:, 5]) > 0.:
                report("{} backend {} max |N_B diff| / N_B".format(
                    name, backend), np.max(np.abs(a[:, 5] - c[:, 5]))
                    / np.max(a[:, 5]), limits[backend])


TESTS = [
    ("conservation_ideal", test_conservation_ideal),
    ("conservation_viscous", test_conservation_viscous),
    ("conservation_outflow", test_conservation_outflow),
    ("lumpy_ideal", test_lumpy_ideal),
    ("lumpy_second_order", test_lumpy_second_order),
    ("lumpy_viscous_Du", test_lumpy_viscous_Du),
    ("lumpy_diffusion", test_lumpy_diffusion),
    ("freezeout_balance", test_freezeout_balance),
    ("beastmode3", test_beastmode3),
    ("rerun", test_rerun),
    ("rerun_2d", test_rerun_2d),
    ("rerun_source", test_rerun_source),
    ("threads", test_threads),
    ("ipglasma_paths", test_ipglasma_paths),
    ("boost_invariance", test_boost_invariance),
    ("bjorken_rk_order", test_bjorken_rk_order),
    ("gubser_face_reconstruction", test_gubser_face_reconstruction),
    ("tatb_source", test_tatb_source),
    ("memory_history", test_memory_history),
    ("gpu_backend", test_gpu_backend),
    ("gpu_backend_viscous", test_gpu_backend_viscous),
    ("gpu_backend_sources", test_gpu_backend_sources),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("tests", nargs="*", help="tests to run (default all)")
    parser.add_argument("--music", default="./MUSIChydro")
    parser.add_argument("--memory-driver", default=None)
    parser.add_argument("--workdir", default=None,
                        help="keep the run directories here")
    parser.add_argument("--threads", type=int, default=None,
                        help="OMP_NUM_THREADS for the runs")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_intermixed_args()

    if args.list:
        for name, _ in TESTS:
            print(name)
        return 0
    names = args.tests or [name for name, _ in TESTS]
    unknown = set(names) - {name for name, _ in TESTS}
    if unknown:
        parser.error("unknown tests: " + ", ".join(sorted(unknown)))
    check_exe = os.path.isfile(args.music) and os.access(args.music, os.X_OK)
    if not check_exe:
        parser.error("MUSIC executable not found: " + args.music)

    workroot = args.workdir or tempfile.mkdtemp(prefix="music_itest_")
    os.makedirs(workroot, exist_ok=True)
    runner = Runner(args.music, args.memory_driver, workroot,
                    args.workdir is not None, args.threads)
    print("work directory: " + workroot)

    failed, skipped = [], []
    for name, func in TESTS:
        if name not in names:
            continue
        print("[ RUN  ] " + name)
        t0 = time.time()
        try:
            func(runner)
            print("[  OK  ] {} ({:.1f} s)".format(name, time.time() - t0))
        except TestSkipped as err:
            skipped.append(name)
            print("[ SKIP ] {}: {}".format(name, err))
        except (TestFailure, subprocess.TimeoutExpired) as err:
            failed.append(name)
            print("[ FAIL ] {}: {}".format(name, err))

    nrun = len(names) - len(skipped)
    print("{} passed, {} failed, {} skipped".format(
        nrun - len(failed), len(failed), len(skipped)))
    if failed:
        print("failed: " + ", ".join(failed))
        print("run directories kept in " + workroot)
    elif args.workdir is None:
        shutil.rmtree(workroot, ignore_errors=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
