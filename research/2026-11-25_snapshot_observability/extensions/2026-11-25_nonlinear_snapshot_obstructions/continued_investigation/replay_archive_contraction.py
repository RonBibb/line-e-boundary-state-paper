"""Reconstruct the archived hidden velocity by the certified secant iteration."""
import json
from pathlib import Path
import platform
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from check_archive_tangent import INPUTS, field_rhs


def replay(grid, noisy, perturb_measured=False):
    with np.load(INPUTS/"bb-cbe3b"/f"cbe3b_slices_N{grid}.npz") as archive:
        initial = np.stack([archive[f"tau18.0_{i}"] for i in range(4)])
        observed = archive["tau16.0_0"].copy()
    x = 2*np.pi*np.arange(grid)/grid
    dx = 2*np.pi/grid
    observed_secant = (initial[0]-observed)/2
    direction = np.cos(5*x)-0.65*np.cos(31*x)+0.37*np.sin(83*x)
    direction /= np.max(np.abs(direction))
    injected = noisy*direction
    known_velocity = initial[1].copy()
    if perturb_measured:
        P0,Q0,R0 = initial[0].copy(),initial[2].copy(),initial[3].copy()
        P_shift = 0.002*np.sin(11*x)
        Q_shift = 0.003*np.exp(-P0)*np.cos(8*x)
        momentum = np.exp(2*P0)*R0+0.00015*np.exp(P0)*np.sin(6*x)
        initial[0] = P0+P_shift
        initial[2] = Q0+Q_shift
        initial[3] = np.exp(-2*initial[0])*momentum
        injected = injected+P_shift/2
    target = observed_secant+injected
    guess = observed_secant.copy() if perturb_measured else target.copy()

    def forward(velocity):
        initial[1] = velocity
        solved = solve_ivp(lambda t,state:field_rhs(t,state.reshape(4,grid),dx).ravel(),
                           (18.,16.),initial.ravel(),method="DOP853",rtol=1e-12,
                           atol=3e-15,max_step=0.04)
        assert solved.success
        return (initial[0]-solved.y[:grid,-1])/2

    rows = []
    for iteration in range(12):
        residual = target-forward(guess)
        rows.append({"iteration": iteration,
                     "normalized_record_residual": float(np.max(np.abs(residual))),
                     "distance_to_archived_velocity": float(np.max(np.abs(guess-known_velocity)))})
        if rows[-1]["normalized_record_residual"] < 1e-12:
            break
        guess += residual
    assert rows[-1]["normalized_record_residual"] < 1e-10
    if not noisy:
        assert rows[-1]["distance_to_archived_velocity"] < 1e-10
    return rows, guess, float(np.max(np.abs(injected)))


if __name__=="__main__":
    result = {}
    for grid in (1024,2048,4096):
        noiseless_history, noiseless,_ = replay(grid,0.)
        noisy_history, noisy,_ = replay(grid,0.008)
        measured_history,measured,record_noise = replay(grid,0.005,perturb_measured=True)
        measured_difference = float(np.max(np.abs(measured-noiseless)))
        assert measured_difference < 1.122*(record_noise+0.000347)
        result[str(grid)] = {"noiseless":noiseless_history, "noisy":noisy_history,
                             "normalized_noisy_velocity_difference":float(
                                 np.max(np.abs(noisy-noiseless))/0.008),
                             "all_record_fields_perturbed":{
                                 "history":measured_history,
                                 "normalized_record_noise_sup":record_noise,
                                 "velocity_difference_sup":measured_difference}}
    print(json.dumps({"platform":platform.platform(),"observation_noise_amplitude":0.008,
                      "results":result},indent=2,allow_nan=False))
