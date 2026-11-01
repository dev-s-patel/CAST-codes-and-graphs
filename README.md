# Causal Automaton Spacetime (CAST) — Numerical Pipeline & Verification Suite

Official computational reproduction suite and Bayesian Markov Chain Monte Carlo (MCMC) pipeline for the companion research series:

* **Part I (Foundations):** *"Causal Automaton Spacetime I: Emergent Spacetime Continuum, Quantum Correlations, and Gravitational Dynamics from Reversible Cellular Networks"*  
  **DOI:** [10.5281/zenodo.23743449](https://doi.org/10.5281/zenodo.23743449) | *Target Journal: Foundations of Physics*

* **Part II (Cosmology):** *"Causal Automaton Spacetime II: Cosmological Perturbation Theory, Planck Likelihood Constraints, and Discrete Lattice Regularization Signatures"*  
  **DOI:** [10.5281/zenodo.23864786](https://doi.org/10.5281/zenodo.23864786) | *Target Journal: JCAP*

---

## 🌌 Overview & Key Physical Results

This repository contains the numerical verification algorithms, Boltzmann transfer components, and Bayesian likelihood sampling pipelines establishing the observational realization of Causal Automaton Spacetime (CAST):

* **Fundamental Cellular Regularizer Scale:** $l_{\rm cell} = 2\sqrt{\ln 2}\, l_P \approx 2.6914 \times 10^{-35}\text{ m}$.
* **Empirical Planck 2018 Fit:** Minimum-variance best-fit log-likelihood $-\ln L_{\min} = 463.18758$ ($\chi^2_{\rm total} = 1002.855$), yielding an absolute global improvement of $\Delta\chi^2 \approx -13.02$ over baseline $\Lambda$CDM ($\chi^2_{\Lambda\mathrm{CDM}} = 1015.88$).
* **Ab-Initio Concordance:** Posterior coupling $\xi_{\rm cast} = 0.05498 \pm 0.03366$ confirms the theoretical boundary prediction ($\xi_{\rm theory} \approx 0.050$) within a $0.15\sigma$ concordance core without post-hoc tuning.
* **Bayesian Model Selection:** Physical single-parameter extension ($N_{\rm par} = 7$, spatial cutoff $k_c = 0.370 \pm 0.136\text{ Mpc}^{-1}$) achieves decisive evidence support of $\Delta\mathrm{AIC} \approx -10.80$ and $\Delta\mathrm{BIC} \approx -5.0$ relative to continuous $\Lambda$CDM.
* **Acoustic Tail Damping:** Finite recombination visibility ($\Delta\chi_{\rm rec} \sim 150\text{--}250\text{ Mpc}$) smears out microscopic mode ringing, isolating a monotonic power suppression ($-1.48\%$ core window, $-2.34\%$ deep tail) projecting to $\ell_c \sim 5180$ directly accessible to CMB-S4, Simons Observatory, and SPT-3G.

---

## ⚡ Terminal Execution Commands

### 0. Clone and Install Dependencies

```bash
git clone https://github.com/dev-s-patel/CAST-codes-and-graphs.git
cd CAST-codes-and-graphs
pip install -r requirements.txt
```

---

### 1. Appendix A.1: Fourth-Order Runge-Kutta Riccati-Bogoliubov Solver

```bash
python run_cast_rk4_riccati_solver.py
```
* **Script:** `run_cast_rk4_riccati_solver.py`
* **Function:** Numerically integrates the coupled first-order Mukhanov-Sasaki Riccati-Bogoliubov system across 50,000 conformal de Sitter time slices with sub-cycle resolution ($\ge 30$ steps/oscillation).
* **Verification:** Validates cycle-averaged adiabatic projection onto the theoretical asymptotic envelope ($|\beta_k^{\rm asymp}| = 2.28 \times 10^{-4}$) while monitoring symplectic unitarity violation ($||\alpha_k|^2 - |\beta_k|^2 - 1|$) at machine precision ceiling $\mathcal{O}(10^{-14})$.

---

### 2. Appendix A.2: 3D Brillouin Zone Density of States Simpson Integrator

```bash
python run_cast_simpson_integrator.py
```
* **Script:** `run_cast_simpson_integrator.py`
* **Function:** Computes the orientation-averaged 3D lattice acoustic dispersion relation over $\mathbb{S}^2$ via 2D composite Simpson quadrature across the first Brillouin zone.
* **Verification:** Validates strict infrared continuum recovery ($\mathcal{T}_{\rm CA}(u \to 0) = 1.0000000026$, absolute error $< 2.6 \times 10^{-9}$), leading-order Taylor series consistency ($1 - \frac{\pi^2}{40}u^2$) with residuals $\mathcal{O}(10^{-6})$, and boundary impedance cutoff damping ($\mathcal{T}_{\rm CA}(1) \approx 0.7774$).

---

### 3. Appendix A.3: Standalone Production CAMB Boltzmann Pipeline

```bash
python run_cast_camb_production.py
```
* **Script:** `run_cast_camb_production.py`
* **Function:** Injects the 3D Brillouin cellular regularizer directly into the CAMB linear Einstein-Boltzmann differential hierarchy via `SplinedInitialPower`, solving full line-of-sight photon transfer functions up to $\ell_{\max} = 5000$.
* **Outputs:** 
  * `cast_camb_deficit.pdf` (Publication-grade vector plot)
  * `cast_camb_deficit.png` (High-resolution 300 DPI preview)

---

### 4. Appendix A.4: Vectorized Cobaya Theory Wrapper

```bash
# Verify component import inside Cobaya ecosystem
python -c "import cast_theory; print('[OK] CAST Cobaya theory module loaded successfully!')"
```
* **Script:** `cast_theory.py`
* **Function:** Subclasses `cobaya.theories.camb.CAMB` to intercept the parameter stream ($k_c, \xi_{\rm cast}$), constructs the regularized primordial spectrum $\mathcal{P}_{\mathcal{R}}^{\rm CAST}(k)$ across a dense 3000-node wavenumber grid ($k \in [10^{-4}, 5.0]\text{ Mpc}^{-1}$), and overrides CAMB's internal power table while locking non-linear matter corrections to background $n_s$.

---

### 5. Appendix A.5: Production Cobaya MCMC Sampling & Posterior Contours

```bash
# Execute production MCMC sampling (MPI multi-core recommended)
cobaya-run cast_mcmc_run.yaml

# Generate 1D posterior distributions and 2D credible contours (Figure 4)
python plot_contours.py
```
* **Configuration:** `cast_mcmc_run.yaml`
* **Dataset:** Co-added, foreground-marginalized Planck 2018 high-$\ell$ likelihood (`plik_lite` $TT,TE,EE$, $\ell \in [30, 2508]$) + low-$\ell$ Commander temperature ($\ell \in [2, 29]$) + low-$\ell$ SimAll $EE$ polarization ($\ell \in [2, 29]$).
* **Sampling Statistics:** Fast-slow Metropolis-Hastings sampling across 4 independent chains reaching 37,091 accepted steps ($N_{\rm eff} = 25,347$, Gelman-Rubin convergence $R - 1 = 0.0291 < 0.03$).
* **Outputs:**
  * `cast_planck_contours.pdf` (Vector contour triangle plot for LaTeX Figure 4)
  * `cast_planck_contours.png` (High-resolution raster plot)

---

## 📦 Requirements

* **Python:** $\ge 3.8$
* `numpy >= 1.24.0`
* `scipy >= 1.10.0`
* `matplotlib >= 3.7.0`
* `camb >= 1.4.0`
* `cobaya >= 3.3.0`
* `getdist >= 1.4.0`
* `pyyaml >= 6.0`

Install all dependencies via pip:
```bash
pip install numpy scipy matplotlib camb cobaya getdist pyyaml
```

---

## 📜 Citations

If you use or verify this computational pipeline, please cite the companion manuscripts:

```bibtex

@article{patel2026cast1,
  title={Causal Automaton Spacetime I: Emergent Spacetime Continuum, Quantum Correlations, and Gravitational Dynamics from Reversible Cellular Networks},
  author={Patel, Dev S.},
  journal={Submitted to Foundations of Physics},
  year={2026},
  doi={10.5281/zenodo.23743449}
}

@article{patel2026cast2,
  title={Causal Automaton Spacetime II: Cosmological Perturbation Theory, Planck Likelihood Constraints, and Discrete Lattice Regularization Signatures},
  author={Patel, Dev S.},
  journal={Submitted to Journal of Cosmology and Astroparticle Physics (JCAP)},
  year={2026},
  doi={10.5281/zenodo.23864786}
}
```

---

## ⚖️ License
This project is open-source and distributed under the [MIT License](LICENSE).