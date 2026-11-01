"""
================================================================================
Causal Automaton Spacetime (CAST) - Production CAMB Boltzmann Pipeline (v2.1)
================================================================================
Author: Dev S. Patel (Premise House / Independent Researcher)
Paper: Causal Automaton Spacetime II: Cosmological Perturbation Theory, 
       Planck Likelihood Constraints, and Discrete Lattice Regularization Signatures
Submission: Physical Review D (2026)

Description:
  Evaluates the cellular boundary impedance imprint on the high-multipole 
  CMB TT angular power spectrum (ell in [2000, 5000]) using full line-of-sight 
  numerical Boltzmann integration via CAMB.
  
  Updated Parameters (Planck 2018 plik_lite + low-ell MCMC Posterior Calibrated):
    - k_c = 0.370 Mpc^-1 (MCMC posterior mean, replacing legacy 0.25 Mpc^-1)
    - xi_coupling = 0.055 (Theoretical boundary impedance coupling constant)
    - ell_c ~ 5180 (Angular multipole projection of lattice cutoff)
================================================================================
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import camb
    from camb import model
except ImportError:
    print("[ERROR] CAMB library is not installed in this Python environment.")
    print("Please install CAMB using: pip install camb")
    sys.exit(1)


def run_cast_camb_production():
    print("=" * 70)
    print("=== CAST PRODUCTION BOLTZMANN PIPELINE (CAMB 2.0+) ===")
    print("=== MCMC CALIBRATED PARAMETERS: k_c = 0.370 Mpc^-1, xi = 0.055 ===")
    print("=" * 70)

    # -------------------------------------------------------------------------
    # 1. Initialize High-Precision Baseline Lambda-CDM Cosmology
    # -------------------------------------------------------------------------
    print("\n[1/5] Initializing High-Precision Planck 2018 Lambda-CDM Baseline...")
    l_max = 5000
    
    pars_lcdm = camb.CAMBparams()
    pars_lcdm.set_cosmology(
        H0=67.4,
        ombh2=0.02237,
        omch2=0.1200,
        mnu=0.06,
        omk=0.0,
        tau=0.0544
    )
    pars_lcdm.InitPower.set_params(
        As=2.1e-9,
        ns=0.9649,
        r=0.0
    )
    pars_lcdm.set_for_lmax(l_max, lens_potential_accuracy=1)
    pars_lcdm.set_accuracy(
        AccuracyBoost=2.0,
        lSampleBoost=2.0,
        lAccuracyBoost=2.0
    )

    print("      Solving linear Boltzmann hierarchy for baseline Lambda-CDM...")
    results_lcdm = camb.get_results(pars_lcdm)
    powers_lcdm = results_lcdm.get_cmb_power_spectra(pars_lcdm, CMB_unit='muK')
    totCL_lcdm = powers_lcdm['total']
    cl_tt_lcdm = totCL_lcdm[2:l_max + 1, 0]
    ell_arr = np.arange(2, l_max + 1)

    # -------------------------------------------------------------------------
    # 2. Construct Dense Primordial Curvature Spectrum Grid P_R(k)
    # -------------------------------------------------------------------------
    print("\n[2/5] Constructing Dense CAST Primordial Perturbation Spectrum (N_k = 5000)...")
    k_min = 1e-4
    k_max = 5.0
    N_k = 5000
    k_grid = np.geomspace(k_min, k_max, N_k)

    # Updated MCMC calibrated parameters
    k_c = 0.370          # Cellular UV cutoff scale in Mpc^-1 (MCMC Mean)
    xi_coupling = 0.055  # Boundary impedance coupling (Theoretical saturation)

    print(f"      Lattice cutoff scale k_c  = {k_c:.3f} Mpc^-1")
    print(f"      Boundary coupling xi     = {xi_coupling:.3f}")

    # Dimensionless normalized wavenumber
    u = k_grid / k_c
    u_clipped = np.clip(u, 1e-6, 1.0)

    # 3D isotropic Brillouin zone dispersion kernel: T_CA(u) = 1 - (pi^2 / 40) * u^2
    K_lattice = (np.pi**2 / 40.0) * (u_clipped**2)
    T_CA_sq = (1.0 - K_lattice)**2

    # Primordial spectra
    P_R_lcdm = 2.1e-9 * (k_grid / 0.05)**(0.9649 - 1.0)
    P_R_cast = P_R_lcdm * (1.0 - xi_coupling * (1.0 - T_CA_sq))

    # -------------------------------------------------------------------------
    # 3. Inject Modified Spectrum into CAMB via SplinedInitialPower
    # -------------------------------------------------------------------------
    print("\n[3/5] Injecting Modified SplinedInitialPower Table into CAMB...")
    pars_cast = camb.CAMBparams()
    pars_cast.set_cosmology(
        H0=67.4,
        ombh2=0.02237,
        omch2=0.1200,
        mnu=0.06,
        omk=0.0,
        tau=0.0544
    )
    pars_cast.set_initial_power_table(k_grid, P_R_cast)
    pars_cast.InitPower.effective_ns_for_nonlinear = 0.9649
    pars_cast.set_for_lmax(l_max, lens_potential_accuracy=1)
    pars_cast.set_accuracy(
        AccuracyBoost=2.0,
        lSampleBoost=2.0,
        lAccuracyBoost=2.0
    )

    # -------------------------------------------------------------------------
    # 4. Solve Full Line-of-Sight Integrals for CAST
    # -------------------------------------------------------------------------
    print("\n[4/5] Solving Line-of-Sight Boltzmann Hierarchy with LSS Visibility Window...")
    results_cast = camb.get_results(pars_cast)
    powers_cast = results_cast.get_cmb_power_spectra(pars_cast, CMB_unit='muK')
    totCL_cast = powers_cast['total']
    cl_tt_cast = totCL_cast[2:l_max + 1, 0]

    # Relative power deviation in percentage
    delta_cl_rel = (cl_tt_cast - cl_tt_lcdm) / cl_tt_lcdm * 100.0

    # Window metrics
    mask_core = (ell_arr >= 3000) & (ell_arr <= 4000)
    mean_core = np.mean(delta_cl_rel[mask_core])

    mask_full = (ell_arr >= 2500) & (ell_arr <= 4500)
    mean_full = np.mean(delta_cl_rel[mask_full])

    mask_deep = (ell_arr >= 4000) & (ell_arr <= 5000)
    mean_deep = np.mean(delta_cl_rel[mask_deep])

    # -------------------------------------------------------------------------
    # 5. Output Benchmark Report & Publication-Grade Plot
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("=== CAST OBSERVATIONAL BENCHMARK REPORT (PRD PART II) ===")
    print("=" * 70)
    print(f"Core Resonance Window  (ell in [3000, 4000]): {mean_core:+.2f}%  [Target: -2.44%]")
    print(f"Full CMB-S4 Window     (ell in [2500, 4500]): {mean_full:+.2f}%  [Target: -2.40%]")
    print(f"Deep High-ell Tail     (ell in [4000, 5000]): {mean_deep:+.2f}%  [Projected Damping]")
    print(f"Lattice Regularizer Projected Scale: ell_c ~ {int(k_c * 14000):d}")
    print("=" * 70)

    # Generate Publication Figure
    print("\n[5/5] Generating Publication-Grade Figure (cast_camb_deficit.png/pdf)...")
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=300)

    # Plot range: ell in [2000, 5000]
    mask_plot = (ell_arr >= 2000) & (ell_arr <= 5000)
    ell_plot = ell_arr[mask_plot]
    delta_plot = delta_cl_rel[mask_plot]

    # Target window shading (CMB-S4 Core)
    ax.axvspan(3000, 4000, color='gold', alpha=0.18, label=r'CMB-S4 Core Target Window ($\ell \in [3000, 4000]$)')

    # Theoretical 1-sigma target band
    ax.axhspan(-2.44 - 0.6, -2.44 + 0.6, color='gray', alpha=0.15, label=r'Theoretical $1\sigma$ Band ($-2.44\% \pm 0.6\%$)')

    # Central target line
    ax.axhline(-2.44, color='dimgray', linestyle='--', linewidth=1.4, alpha=0.85, label=r'CAST Central Deficit Target ($-2.44\%$)')

    # Null baseline
    ax.axhline(0.0, color='black', linestyle='-', linewidth=0.9, alpha=0.7)

    # CAST Full Line-of-Sight Curve
    ax.plot(
        ell_plot, delta_plot,
        color='#1b7837', linewidth=2.4,
        label=r'CAST Full Boltzmann Hierarchy (CAMB $k_c=0.370\,\mathrm{Mpc}^{-1}$)'
    )

    ax.set_xlim(2000, 5000)
    ax.set_ylim(-4.0, 0.8)
    ax.set_xlabel(r'Multipole Moment $\ell$', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'$\Delta C_\ell^{TT} / C_\ell^{\Lambda\mathrm{CDM}}\ [\%]$', fontsize=13, fontweight='bold')
    ax.set_title(r'Machine-Integrated CAST Imprint on CMB Temperature Anisotropy ($\ell \in [2000, 5000]$)',
                 fontsize=13, pad=12, fontweight='bold')

    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='lower left', frameon=True, framealpha=0.92, facecolor='white', fontsize=10.5)

    plt.tight_layout()
    
    png_filename = "cast_camb_deficit.png"
    pdf_filename = "cast_camb_deficit.pdf"
    
    fig.savefig(png_filename, dpi=300)
    fig.savefig(pdf_filename)
    print(f"      [OK] Saved high-res raster plot: {png_filename}")
    print(f"      [OK] Saved vector PDF plot:     {pdf_filename}")
    plt.close(fig)

    print("\n=== Pipeline Execution Completed Successfully ===")


if __name__ == '__main__':
    run_cast_camb_production()