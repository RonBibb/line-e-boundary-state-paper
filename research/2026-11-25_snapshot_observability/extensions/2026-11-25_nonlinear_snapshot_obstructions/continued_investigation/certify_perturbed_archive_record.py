"""Extend the interval certificate to perturbations in all stored fields."""
import json
import platform

from certify_archive_tangent import (
    INPUTS,DURATION,I,bootstrap,comparison_flow,comparison_matrix,downward,
    tangent_bound,upward,validate_initial_data,iv,
)

P_NOISE, Q_SCALED_NOISE, MOMENTUM_SCALED_NOISE = I("0.002"),I("0.003"),I("0.00015")
RADIUS = I("0.06")
ARCHIVE_CENTER_MISMATCH = I("0.002610")


def certified_case(n):
    archive = INPUTS/"bb-cbe3b"/f"cbe3b_slices_N{n}.npz"
    initial = validate_initial_data(archive,n,RADIUS,
                                    (P_NOISE,Q_SCALED_NOISE,MOMENTUM_SCALED_NOISE))
    enclosed = bootstrap(initial)
    derivative = tangent_bound(enclosed)
    # Along a line in measured variables (P,Q,M=e^{2P}Q_tau), the
    # scaled derivative of Q_tau is e^{-P} dM - 2B dP.
    initial_y = iv.exp(P_NOISE)*Q_SCALED_NOISE
    initial_z = iv.exp(P_NOISE)*MOMENTUM_SCALED_NOISE+2*initial["B"]*P_NOISE
    terminal,*_ = comparison_flow(comparison_matrix(enclosed),
                                   [P_NOISE,I(0),initial_y,initial_z])
    data_change = (terminal-P_NOISE)/DURATION
    defect = I(derivative["observation_defect_upper_bound"])
    safe_additional_error = (1-defect)*RADIUS-ARCHIVE_CENTER_MISMATCH-data_change
    assert safe_additional_error.a > 0
    return {"bounds_for_data_change":upward(data_change),
            "still_available_normalized_secant_error":downward(safe_additional_error),
            "inverse_norm_upper_bound":derivative["inverse_norm_upper_bound"],
            "derivative_defect_upper_bound":derivative["observation_defect_upper_bound"]}


if __name__=="__main__":
    result = {"platform":platform.platform(),"velocity_radius":"0.06",
              "P_sup_radius":"0.002","weighted_Q_sup_radius":"0.003",
              "weighted_momentum_sup_radius":"0.00015",
              "results":{str(n):certified_case(n) for n in (1024,2048,4096)}}
    print(json.dumps(result,indent=2,allow_nan=False))
