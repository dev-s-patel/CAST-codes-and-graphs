import numpy as np

def run_cast_rk4_riccati_solver():
    print("=== CAST PRODUCTION RK4 RICCATI-BOGOLIUBOV SOLVER ===")
    print("Resolving Sub-Planckian Mode Emergence with Sub-Cycle Stepping...")

    # Scaled Benchmark Parameters for Exact Numerical Resolution
    H_star = 1.0e-3         
    M_subP = 1.0            
    u_val = 0.50            # Normalized wavenumber u = k / k_c
    
    # Boundary crossing at a_0 = 1.0; eta_0 = -1 / H_star
    eta_in = - 1.2 / H_star   # Prior to emergence threshold
    eta_out = - 0.2 / H_star  # Asymptotic post-emergence regime
    
    # Explicit sub-cycle RK4 stepping (>= 30 steps per oscillation)
    N_steps = 50000
    d_eta = (eta_out - eta_in) / N_steps
    
    def dispersion(eta):
        arg = - 0.5 * H_star * u_val * eta
        arg = np.clip(arg, 0.01, np.pi/2 - 0.01)
        omega = 2.0 * M_subP * np.sin(arg)
        d_omega = - H_star * u_val * M_subP * np.cos(arg)
        return omega, d_omega

    Y = np.array([1.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j], dtype=complex)
    
    def derivatives(eta, state):
        alpha, beta, theta = state[0], state[1], state[2].real
        omega, d_omega = dispersion(eta)
        g = d_omega / (2.0 * omega)
        
        d_alpha = g * beta * np.exp(2.0j * theta)
        d_beta  = g * alpha * np.exp(-2.0j * theta)
        d_theta = omega + 0.0j
        return np.array([d_alpha, d_beta, d_theta], dtype=complex)

    print(f"Executing RK4 across {N_steps} slices (Step size d_eta = {d_eta:.4f})...")
    eta = eta_in
    for step in range(N_steps):
        k1 = derivatives(eta, Y)
        k2 = derivatives(eta + 0.5 * d_eta, Y + 0.5 * d_eta * k1)
        k3 = derivatives(eta + 0.5 * d_eta, Y + 0.5 * d_eta * k2)
        k4 = derivatives(eta + d_eta, Y + d_eta * k3)
        Y += (d_eta / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        eta += d_eta

    alpha_final = Y[0]
    beta_final = Y[1]
    
    norm_check = np.abs(alpha_final)**2 - np.abs(beta_final)**2
    norm_error = np.abs(norm_check - 1.0)
    
    sinc_u = np.sinc(u_val)
    analytic_beta = (H_star / M_subP) * np.abs(1.0 - sinc_u) / (2.0 * np.sqrt(sinc_u))
    
    print("\n================ RK4 CONVERGENCE REPORT ================")
    print(f"Final |alpha_k|             : {np.abs(alpha_final):.8f}")
    print(f"Final |beta_k|  (Numerical) : {np.abs(beta_final):.6e}")
    print(f"Asymptotic Modulus (Theory) : {analytic_beta:.6e}")
    print(f"Symplectic Invariance Norm  : {norm_check.real:.12f} (Target: 1.000000000000)")
    print(f"Unitary Violation Error     : {norm_error:.2e}")
    print("========================================================\n")

if __name__ == '__main__':
    run_cast_rk4_riccati_solver()