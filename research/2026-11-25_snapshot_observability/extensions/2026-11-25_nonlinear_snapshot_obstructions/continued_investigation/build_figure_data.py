"""Export exact and measured scientific series for the manuscript PGFPlots."""
import json
from pathlib import Path
import numpy as np

from check_null_transport import tuned_collision

ROOT = Path(__file__).resolve().parent
DESTINATION = ROOT/"figures"


def read(name):
    return json.loads((ROOT/name).read_text())


def table(filename, header, rows):
    np.savetxt(DESTINATION/filename,np.asarray(rows),delimiter=",",header=header,
               comments="",fmt="%.16g")


def main():
    model, b_delay = tuned_collision()
    a, b = model.a,model.a+b_delay
    sample = np.linspace(a,b,321)
    rows = []
    for t in sample:
        shift = model.P0+np.log(t)
        plus, minus = (model.exact_fields(t,np.array([0.]),s)[0,0]+shift for s in (1.,-1.))
        rows.append(((t-a)/b_delay,plus,minus))
    table("collision_path.csv","fraction,branch_a,branch_b",rows)

    record = read("null_transport_checks.json")
    table("null_gains.csv","mode,nodal,placed,principal",[
        (int(k),m[1]["full_gain_for_sine_primitive"],m[0]["full_gain_for_sine_primitive"],
         m[0]["principal_gain_for_sine_primitive"])
        for k,m in record["high_mode_results"].items()])

    record = read("full_metric_holonomy_checks.json")["explicit_short_delay"]
    table("full_metric_gains.csv","mode,orbit,metric",[
        (row["mode"],row["orbit_gain"],row["full_metric_gain"])
        for row in record["high_modes"]])

    record = read("archive_contraction_replay.json")
    curve = record["results"]["4096"]["noiseless"]
    initial_error = curve[0]["normalized_record_residual"]
    table("archive_iteration.csv","iteration,residual,certified_upper",[
        (row["iteration"],row["normalized_record_residual"],initial_error*0.10725**row["iteration"])
        for row in curve])

    record = read("production_modulus_checks.json")
    truth = float(record["independent_continuum_reference"][-1]["modulus"])
    table("modulus_convergence.csv","grid,zero_tail_depth8,bessel_tail_depth8,bessel_reference",[
        (int(n),
         abs(row["estimates"]["8"]["earlier_anchor_constant"]-truth),
         abs(row["estimates"]["8"]["tail_corrected_matched_stencil"]["modulus"]-truth),
         abs(row["discrete_exact_reference"]["modulus"]-truth))
        for n,row in record["grids"].items()])


if __name__ == "__main__":
    main()
