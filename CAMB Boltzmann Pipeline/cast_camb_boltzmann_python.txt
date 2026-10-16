import numpy as np
import camb
from camb import model

def run_cast_camb_production():
    print("=== CAST PRODUCTION BOLTZMANN PIPELINE (CAMB 2.0.4) ===")
    print("Executing High-Precision (Double Accuracy) Boltzmann Hierarchy...")

    # 1. Initialize High-Precision Lambda-CDM Baseline
    print("1. Initializing High-Precision Lambda-CDM Baseline...")
    l_max = 5000
    pars_lcdm = camb.CAMBparams()
    pars_lcdm.set_cosmology(
        H0=67.4, ombh2=0.02237, omch2=0.1200,
        mnu=0.06, omk=0, tau=0.0544
    )
    pars_lcdm.InitPower.set_params(As=2.1e-9, ns=0.9649, r=0)
    pars_lcdm.set_for_lmax(l_max, lens_potential_accuracy=1)
    pars_lcdm.set_accuracy(
        AccuracyBoost=2.0, lSampleBoost=2.0, lAccuracyBoost=2.0
    )

    results_lcdm = camb.get_results(pars_lcdm)
    powers_lcdm = results_lcdm.get_cmb_power_spectra(
        pars_lcdm, CMB_unit='muK'
    )
    totCL_lcdm = powers_lcdm['total']
    cl_tt_lcdm = totCL_lcdm[2:l_max+1, 0]
    ell_arr = np.arange(2, l_max+1)

    # 2. Construct Dense Primordial P_R(k) Grid (N_k = 5000)
    print("2. Constructing Dense Primordial P_R(k) Grid (N_k = 5000)...")
    k_min = 1e-4
    k_max = 5.0
    N_k = 5000
    k_grid = np.geomspace(k_min, k_max, N_k)
    k_c = 0.25           # Sub-Planckian UV cutoff scale in Mpc^-1
    xi_coupling = 0.055  # Boundary impedance coupling constant

    u = k_grid / k_c
    u_clipped = np.clip(u, 1e-6, 1.0)
    K_lattice = (np.pi**2 / 40.0) * (u_clipped**2)
    T_CA_sq = (1.0 - K_lattice)**2

    # Note on Riccati Phase Modulation:
    # Under full line-of-sight Boltzmann transfer, the microscopic
    # oscillatory mode term (+2|beta_k|cos(phi_k)) from Eq. (156) is
    # heavily damped by Last Scattering Surface thickness projection
    # smearing (<cos(phi_k)> ~ 0, Section 7.3, Item 2), leaving the
    # envelope impedance damping (1 - xi*(1 - T_CA^2)) as net transfer.
    P_R_lcdm = 2.1e-9 * (k_grid / 0.05)**(0.9649 - 1.0)
    P_R_cast = P_R_lcdm * (1.0 - xi_coupling * (1.0 - T_CA_sq))

    # 3. Inject into CAMB Fortran Engine via SplinedInitialPower
    print("3. Injecting into CAMB Fortran Engine via SplinedInitialPower...")
    pars_cast = camb.CAMBparams()
    pars_cast.set_cosmology(
        H0=67.4, ombh2=0.02237, omch2=0.1200,
        mnu=0.06, omk=0, tau=0.0544
    )
    pars_cast.set_initial_power_table(k_grid, P_R_cast)
    pars_cast.InitPower.effective_ns_for_nonlinear = 0.9649
    pars_cast.set_for_lmax(l_max, lens_potential_accuracy=1)
    pars_cast.set_accuracy(
        AccuracyBoost=2.0, lSampleBoost=2.0, lAccuracyBoost=2.0
    )

    # 4. Solving High-Precision Boltzmann Hierarchy & Transfer
    print("4. Solving High-Precision Boltzmann Integrals & Transfer...")
    results_cast = camb.get_results(pars_cast)
    powers_cast = results_cast.get_cmb_power_spectra(
        pars_cast, CMB_unit='muK'
    )
    totCL_cast = powers_cast['total']
    cl_tt_cast = totCL_cast[2:l_max+1, 0]

    # 5. Evaluating Relative Multipole Response
    print("5. Evaluating Observational Multipole Response (Delta C_l / C_l)...")
    delta_cl_rel = (cl_tt_cast - cl_tt_lcdm) / cl_tt_lcdm * 100.0
    mask_core = (ell_arr >= 3000) & (ell_arr <= 4000)
    mean_core = np.mean(delta_cl_rel[mask_core])
    mask_full = (ell_arr >= 2500) & (ell_arr <= 4500)
    mean_full = np.mean(delta_cl_rel[mask_full])

    print("\n================== BENCHMARK REPORT ==================")
    print(f"Core Resonance Window (l in [3000, 4000]): {mean_core:.2f}% "
          f"(Target: ~ -2.4%)")
    print(f"Full CMB-S4 Window (l in [2500, 4500]): {mean_full:.2f}% "
          f"(Target Band: -2.2% to -3.4%)")
    print("======================================================\n")

if __name__ == '__main__':
    run_cast_camb_production()