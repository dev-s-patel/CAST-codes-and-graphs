\# Causal Automaton Spacetime (CAST) — Numerical Pipeline \& Verification Suite



Official computational reproduction suite for the preprint:  

\*\*"Causal Automaton Spacetime: An Emergent Cellular Framework for Spacetime Dynamics, Quantum Correlations, and Cosmological Falsification"\*\*  

\*Dev S. Patel (2026)\* | \*\*DOI:\*\* \[10.5281/zenodo.22916869](https://doi.org/10.5281/zenodo.22916869)



\---



\## ⚡ Terminal Execution Commands



\### 0. Clone and Install Dependencies

```bash

git clone \[https://github.com/dev-s-patel/CAST-codes-and-graphs.git](https://github.com/dev-s-patel/CAST-codes-and-graphs.git)

cd CAST-codes-and-graphs

pip install -r requirements.txt

```



\---



\### 1. Appendix A.1: Riccati-Bogoliubov Solver

```bash

cd "Riccati-Bogoliubov Solver"

python run\_cast\_rk4\_riccati\_solver.py

cd ..

```

\* \*\*Function:\*\* Evaluates sub-Planckian mode emergence and monitors symplectic unitarity violation ($\\Delta \\mathcal{U} < 10^{-15}$).

\* \*\*Output:\*\* Generates `cast\_riccati\_rk4\_plot.png`.



\---



\### 2. Appendix A.2: Simpson Integrator

```bash

cd "Simpson Integrator"

python run\_cast\_simpson\_integrator.py

cd ..

```

\* \*\*Function:\*\* Computes 3D Brillouin zone orientation-averaged density of states and IR recovery ($T\_{\\mathrm{CA}}(u \\to 0) = 1.000$).

\* \*\*Output:\*\* Generates `cast\_lattice\_transfer\_plot.png`.



\---



\### 3. Appendix A.3: Production CAMB Boltzmann Pipeline

```bash

cd "CAMB Boltzmann Pipeline"

python run\_cast\_camb\_production.py

cd ..

```

\* \*\*Function:\*\* Integrates sub-Planckian boundary impedance into CAMB Boltzmann engine to compute high-$\\ell$ CMB power spectrum angular deficit ($\\Delta C\_\\ell / C\_\\ell \\approx -2.44\\%$).

\* \*\*Output:\*\* Generates `cast\_camb\_deficit.png`.



\---



\## 📦 Requirements

\* Python >= 3.8

\* numpy >= 1.24.0

\* scipy >= 1.10.0

\* matplotlib >= 3.7.0

\* camb >= 1.4.0

