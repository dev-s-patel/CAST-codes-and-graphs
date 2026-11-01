import numpy as np
from cobaya.theories.camb.camb import CAMB, CAMBOutputs

class cast_theory(CAMB):
    r"""
    Sub-Planckian Cellular Automata Spacetime (CAST) Theory Wrapper.
    Injects 3D Brillouin lattice dispersion transfer function into CAMB initial power table.
    """
    def get_can_support_params(self):
        return list(super().get_can_support_params()) + ["k_c", "xi_cast"]

    def calculate(self, state, want_derived=True, **params_values_dict):
        # 1. CAST parameters pop karo taaki standard CAMB power args mein leak na hon
        kc = params_values_dict.pop("k_c", 0.25)
        xi = params_values_dict.pop("xi_cast", 0.055)

        # 2. Transfer calculations fetch karo
        params, results = self.provider.get_CAMB_transfers()

        # 3. Standard baseline primordial parameters collect karo
        args = {
            self.translate_param(p): v
            for p, v in params_values_dict.items()
            if p in self.power_params
        }
        args.update(self.initial_power_args)

        # Primordial values seedha params dictionary se read karo (SplinedInitialPower attribute error se bachne ke liye)
        ns = params_values_dict.get("ns", args.get("ns", 0.9649))
        if "As" in params_values_dict:
            As = params_values_dict["As"]
        elif "logA" in params_values_dict:
            As = np.exp(params_values_dict["logA"]) * 1e-10
        else:
            As = args.get("As", 2.1e-9)
        pivot = args.get("pivot_scalar", 0.05)

        # 4. Appendix A.3: Sub-Planckian Brillouin dispersion modulation construct karo
        ks = np.geomspace(1e-4, 5.0, 3000)
        P_0 = As * (ks / pivot) ** (ns - 1.0)

        # Brillouin cutoff impedance
        u = np.clip(ks / kc, 1e-6, 1.0)
        K_lat = (np.pi**2 / 40.0) * (u**2)
        T_CA_sq = (1.0 - K_lat)**2
        P_mod = P_0 * (1.0 - xi * (1.0 - T_CA_sq))

        # Initial power table override via CAMBparams
        results.Params.set_initial_power_table(ks, P_mod)
        results.Params.InitPower.effective_ns_for_nonlinear = ns

        # 5. Non-linear settings
        if self.non_linear_sources or self.non_linear_pk:
            if self.use_non_linear_ratio:
                non_linear_ratio = self.provider.get_non_linear_ratio(results)
                results.Params.NonLinearModel.set_ratio(
                    non_linear_ratio["k_h"],
                    non_linear_ratio["z"],
                    non_linear_ratio["ratio"]
                )
            else:
                nl_args = {
                    self.translate_param(p): v
                    for p, v in params_values_dict.items()
                    if p in self.nonlin_params
                }
                nl_args.update(self.nonlin_args)
                results.Params.NonLinearModel.set_params(**nl_args)

        # 6. Final power spectra solve karo
        try:
            results.power_spectra_from_transfer()

            if "sigma8" in params_values_dict and hasattr(results.Params.InitPower, "As"):
                sigma8 = results.get_sigma8_0()
                results.Params.InitPower.As *= (params_values_dict["sigma8"]**2 / sigma8**2)
                results.power_spectra_from_transfer()

            # 7. Likelihood collectors ko Cl pass karo
            for product, collector in self.collectors.items():
                if collector:
                    state[product] = collector.method(results, *collector.args, **collector.kwargs)
                    if collector.post:
                        state[product] = collector.post(*state[product])
                else:
                    state[product] = results.copy()
        except self.camb.baseconfig.CAMBError as e:
            if self.stop_at_error:
                self.log.error(
                    "Computation error!\nParameters: %r, %r",
                    dict(state["params"]),
                    dict(self.extra_args),
                )
                raise
            else:
                return False

        # 8. Derived aur extra parameters populate karo
        intermediates = CAMBOutputs(
            params, results, results.get_derived_params() if results else None
        )
        if want_derived:
            state["derived"] = self._get_derived_output(intermediates)
        state["derived_extra"] = {
            p: self._get_derived(p, intermediates) for p in self.derived_extra
        }
