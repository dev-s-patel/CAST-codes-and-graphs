import numpy as np
from scipy.integrate import simpson

def run_cast_simpson_integrator():
    print("=== CAST PRODUCTION 3D BRILLOUIN ZONE SIMPSON INTEGRATOR ===")
    print("Evaluating Orientation-Averaged Sub-Planckian Lattice Transfer Kernel...")

    N_theta = 120   # Even intervals for Simpson rule
    N_phi = 240     
    
    theta = np.linspace(1e-6, np.pi - 1e-6, N_theta + 1)
    phi = np.linspace(0.0, 2.0 * np.pi, N_phi + 1)
    
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    sin_theta = np.sin(THETA)
    
    nx = np.sin(THETA) * np.cos(PHI)
    ny = np.sin(THETA) * np.sin(PHI)
    nz = np.cos(THETA)
    
    u_values = np.array([1e-5, 0.10, 0.25, 0.50, 0.75, 1.00])
    T_CA_results = []

    print(f"\nIntegrating across {N_theta} x {N_phi} grid points per wavenumber slice...")
    print("----------------------------------------------------------------------")
    print(f"{'Normalized u (k/k_c)':<22} | {'Numerical T_CA(u)':<20} | {'Status'}")
    print("----------------------------------------------------------------------")

    for u in u_values:
        arg_x = 0.5 * np.pi * u * nx
        arg_y = 0.5 * np.pi * u * ny
        arg_z = 0.5 * np.pi * u * nz
        
        dispersion_sum = np.sin(arg_x)**2 + np.sin(arg_y)**2 + np.sin(arg_z)**2
        v_phase = (2.0 / (np.pi * u)) * np.sqrt(dispersion_sum)
        integrand = v_phase * sin_theta
        
        integral_phi = simpson(integrand, x=phi, axis=1)
        integral_2D = simpson(integral_phi, x=theta, axis=0)
        
        T_CA_val = integral_2D / (4.0 * np.pi)
        T_CA_results.append(T_CA_val)
        
        status = "IR Continuum Limit" if u < 1e-4 else "Lattice Dispersion Damping"
        print(f"{u:<22.5f} | {T_CA_val:<20.8f} | {status}")

    print("----------------------------------------------------------------------")
    
    u_test = 0.10
    analytic_u01 = 1.0 - (np.pi**2 / 40.0) * (u_test**2)
    numerical_u01 = T_CA_results[1]
    discrepancy = np.abs(numerical_u01 - analytic_u01)

    print("\n================ SIMPSON CONVERGENCE REPORT ================")
    print(f"Infrared Limit T_CA(u -> 0)  : {T_CA_results[0]:.10f} (Target: 1.0000000000)")
    print(f"IR Absolute Error            : {np.abs(T_CA_results[0] - 1.0):.2e}")
    print(f"Order u^2 Analytic Benchmark : {analytic_u01:.8f}")
    print(f"Order u^2 Numerical Result   : {numerical_u01:.8f}")
    print(f"Taylor Consistency Residual  : {discrepancy:.2e}")
    print("============================================================\n")

if __name__ == '__main__':
    run_cast_simpson_integrator()