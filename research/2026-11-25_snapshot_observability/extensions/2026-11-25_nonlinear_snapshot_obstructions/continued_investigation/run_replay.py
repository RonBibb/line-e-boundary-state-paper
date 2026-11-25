import concurrent.futures
import json
from pathlib import Path
import subprocess
import sys
import time

root=Path(__file__).resolve().parent
here=root/"build"/"replay"
here.mkdir(parents=True,exist_ok=True)
runs=[
    ("check_nonlinear_sampling.py","nonlinear_sampling_checks.json"),
    ("check_null_transport.py","null_transport_checks.json"),
    ("check_open_two_snapshot.py","open_two_snapshot_checks.json"),
    ("reconstruct_fixed_snapshots.py","fixed_snapshot_reconstruction.json"),
    ("check_dispersion_correction.py","dispersion_correction_checks.json"),
    ("check_both_velocity_holonomy.py","both_velocity_holonomy_checks.json"),
    ("check_full_metric_holonomy.py","full_metric_holonomy_checks.json"),
    ("check_low_velocity_tail.py","low_velocity_checks.json"),
    ("check_ernst_both_velocities.py","ernst_both_velocity_checks.json"),
    ("check_production_modulus.py","production_modulus_checks.json"),
    ("check_circular_loop.py","circular_loop_checks.json"),
    ("certify_asymmetry.py","asymmetry_certificate.json"),
    ("check_asymptotic_geometry.py","asymptotic_geometry_checks.json"),
    ("check_extended_predictions.py","extended_predictions.json"),
    ("certify_archive_tangent.py","archive_tangent_certificate.json"),
    ("certify_perturbed_archive_record.py","perturbed_archive_certificate.json"),
    ("replay_archive_contraction.py","archive_contraction_replay.json"),
    ("../check_obstruction.py","parent_obstruction_checks.json"),
    ("../check_exact_collision.py","parent_exact_collision_checks.json"),
    ("../check_rank_two.py","parent_rank_two_checks.json"),
    ("../check_cubic.py","parent_cubic_checks.json"),
    ("../check_archive_tangent.py","parent_archive_tangent_checks.json"),
    ("../../../check_geometry.py","parent_geometry_checks.json"),
]


def execute(task):
    driver,result=task
    start=time.monotonic()
    temporary=here/(result+".tmp")
    with temporary.open("w") as output,(here/(result+".err")).open("w") as error:
        run=subprocess.run([sys.executable,"-B",str(root/driver)],cwd=root,
                           stdout=output,stderr=error,check=False)
    valid=False
    if run.returncode==0:
        try:
            saved=json.loads(temporary.read_text())
            valid=saved["platform"].startswith("Linux")
        except Exception:
            pass
    if valid:
        temporary.rename(here/result)
    return {"driver":driver,"evidence":result,"seconds":round(time.monotonic()-start,2),
            "linux_json":valid,"exit_code":run.returncode}


with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results=[]
    tasks=[pool.submit(execute,task) for task in runs]
    for future in concurrent.futures.as_completed(tasks):
        result=future.result()
        results.append(result)
        print(json.dumps(result),flush=True)
        (here/"SUMMARY.json").write_text(json.dumps(results,indent=2)+"\n")
if any(not entry["linux_json"] for entry in results):
    sys.exit(1)
