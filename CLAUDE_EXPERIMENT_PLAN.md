# Claude experiment plan for two metric slices in vacuum Gowdy spacetimes

**Exploration goal.** Determine when two measured spatial geometries recover both missing initial velocities in a vacuum Gowdy spacetime. The working manuscript gives exact pairs that agree at two observation times and a stable inverse for a circular-loop family at other times. A valuable next result would explain which genuinely inhomogeneous backgrounds admit the same recovery, or exhibit a background where a velocity direction remains invisible even though the high frequencies are controlled. A prior theorem that already settles the question, or a counterexample that forces a correction, is a successful outcome.

Claude, use the priorities below as a starting point. Add tests, change their order, or abandon a route when a finding gives you a better one. Record what changed your decision. The prior-art and manuscript-annotation audit is part of the assignment whichever route you take. Run **all experiments on `linux.vm`**, including numerical runs, symbolic checks, interval calculations, and reproductions.

Work from the November continuation. Its TeX master includes sections from the parent directory. The August E1 manuscript and its citation matrix document the earlier claim about a finite-window reconstruction. Use the correction ledger when comparing the two versions. [November manuscript][pdf] · [TeX master][master] · [Research readout][readout] · [Correction ledger][corrections]

## Experiments in suggested priority order

### 1. Audit prior art against the actual manuscript

List the claims in the abstract, introduction, and theorem statements, then follow the master's TeX includes. For each claim, ask what was already known for the **same observations, unknowns, admissible data, and stability norm**. Search across mathematical relativity, inverse wave problems, and control theory. Useful query pairs include “Gowdy” with “two time inverse problem” or “Gowdy-to-Ernst”, and “finite-instant observability” or “two-time snapshots” with “coupled waves”, “polarization transport”, or “Gramian”. Include restricted singular-limit or scalar results wherever the current draft uses them.

Start from the November review and historical citation matrix. Search arXiv and a scholarly citation index, follow references backward and citations forward, and inspect the closest primary theorems. A first check is Wang's finite-instant observability survey, which appears in the November review and has no bibliography entry in the November master. Trace the most relevant underlying theorem. Also check whether the introduction's citation to Cui, Laurent and Wang tells a reader of the metric theorem exactly which part is inherited. [November prior-art notes][prior] · [Historical citation matrix][matrix]

Return a dated claim-to-source table. For each relevant prior, give its identifier, theorem or page, overlap and difference, and the exact TeX file and passage where it is annotated. Add missing attribution at the point of use and supply the bibliography entry. If a source subsumes a result, narrow the novelty language and identify what still survives. Rebuild the manuscript on the VM, check that citations resolve, and read the changed passages. Log queries, databases, exclusions that affect the conclusion, and inaccessible sources. State coverage as of the search date.

### 2. Try to break the two arguments the paper depends on

- **Stable recovery.** Derive the extra measurement independently from the momentum constraint, `λ_θ = −2t⟨U_t,U_θ⟩`, where `U=(P,Q)`. Check the transported projection, derivative counts, rational cycle closure, the irrational weak-mode construction, and the short-delay estimate that controls every mode. Use the circular loop with `a=1.6`, winding 4, initial radius 0.32, and zero radial velocity. At delay 0.02 the saved example includes a nearly invisible orbit-only mode at frequency 1575. Compare orbit fields alone and the full metric on that same mode. Use the longer placement with accumulated frame rotation π as the adverse control. Identify precisely which inequality would fail if a proof step is wrong. [Orbit observations][orbit] · [Full metric proof][metric] · [Replay drivers][replay]
- **Exact nonuniqueness.** Reconstruct the two-mode Ernst example with frequencies 6 and 13. Verify periodicity, agreement of both induced metrics and the known velocity component on both slices, and the differing initial velocity. Then allow the quotient-circle shifts, reflections and changes of Killing basis admitted by the theorem. Check the mean-curvature argument with the corrected induced-metric determinant. Keep the time orientation and boundary hypotheses explicit. The circular-loop pair exchanged by reflection is a useful control for how marking changes the conclusion. [Rank-two construction][rank] · [Circular-loop comparison][circle]

### 3. Resolve the finite-dimensional ambiguity on an inhomogeneous background

This is the first new theorem to pursue if the preceding checks hold. The transported metric projections can control high frequencies while leaving a finite-dimensional kernel. The paper removes that kernel in the circular case and already covers nearby data. Try the exact inhomogeneous Ernst family at its common initial node, with the second slice moved to a short, nonresonant delay. The short-delay irrational result gives a reason to start there. Withhold **both** initial velocities and use exactly two full induced metrics.

1. Specify the background, marking, time gap, norms and periodic momentum constraint. Check the high-frequency Gramian along the whole circle. Retain both metric profiles, including their means with a common initial anchor, when testing the full observation.
2. Compute the smallest gains and candidate invisible directions for the exact linearized equations in successively larger Fourier truncations. Include the constant modes and coupling between frequencies. Refine the time integration independently. Track the direction as well as its gain to see whether it converges to a low mode or escapes toward the cutoff.
3. If a positive lower bound persists, seek a proof coupling the finite-mode estimate to a controlled high-frequency remainder, or adapt the circular backward estimate. State which geometric assumption provides the bound and how it extends to the nonlinear map. If a gain collapses, try to construct an admissible continuum kernel and test it against both absolute metric profiles. An exact nonlinear collision requires two distinct histories with identical data.

Extra early orbit coefficients are an existing repair. If you need them, specify the additional measurements and which ambiguity each resolves. A stable finite matrix, a weak sequence at increasing frequency, a linearized kernel, and an exact nonlinear collision support four different conclusions. [Metric criterion and circular estimate][metric] · [Earlier recovery with a measured velocity component][ernst]

### 4. Determine whether the gain survives finite measurement accuracy

If injectivity survives, test how well the velocities can actually be recovered. Perturb the measured metric potential before differentiating it. Compare orbit-only and full-metric reconstructions at matched noise, resolution, and delay, including placements slightly off the long-delay resonance. Report the input noise norm, bandwidth, differentiation or regularization rule, and velocity error. Check whether the gain persists under refinement. A bandwidth-dependent error law would clarify what a finite-resolution measurement buys even if the general injectivity question stays open.

### 5. Follow the result to an observable at the singular limit

If a new stable background satisfies the manuscript's low-velocity smallness condition, determine which asymptotic fields its two slices control, in which norms, and with which metric anchor. If an ambiguity survives instead, test whether the invisible direction changes the limiting speed or the limiting metric potential. The draft already has pairs that separate these targets. Use them as controls and articulate what the new example adds. [Low-velocity theorem][tail] · [Equal-speed example][modulus]

## Linux baseline and experiment records

Connect with `ssh linux.vm`. If its configured address times out on the research Mac, the following route uses the running Parallels guest's current address and checks the saved SSH key for `linux.vm`.

```sh
guest_ip=$(prlctl exec 'Ubuntu 24 ARM64' hostname -I | awk '{print $1}')
ssh -o HostName="$guest_ip" -o HostKeyAlias=linux.vm -o StrictHostKeyChecking=yes linux.vm
```

In the guest, use a clean checkout of `RonBibb/line-e-boundary-state-paper`. Run the commands below from the checkout root to record the source and environment and reproduce the controls. The research environment on the guest matches the November dependency list. [Pinned dependencies][deps]

```sh
hostname
uname -sm
git rev-parse HEAD
. /home/parallels/codex-p344-audit-venv/bin/activate
python -V
python -m pip freeze | grep -Ei '^(numpy|scipy|mpmath|sympy)=='
(cd reproduction && sha256sum -c SHA256SUMS)
cd research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation
OPENBLAS_NUM_THREADS=1 python -B run_replay.py
OPENBLAS_NUM_THREADS=1 python -B check_long_delay_metric.py > build/replay/long_delay_metric_checks.json
```

The baseline is 23 replay records plus the separate long-delay control. Compare the scientific payloads under `build/replay/` with the committed JSON. For the six `parent_` outputs, the reference is beside the listed driver with the prefix removed. The runner checks Linux output and successful execution. Inspect the values, bounds and residuals themselves. Keep the frozen inputs intact and use a scratch checkout for historical producers that write beside their source. [Reproduction instructions][replay] · [Frozen inputs][inputs]

## What to return

Keep new scripts, commands, parameters, machine-readable evidence and a short research report together in a dated folder under `research/`. Include the citation table and manuscript edits. Lead the report with the strongest claim Ron can defend and what it changes in the paper. Mark each finding as proved, numerically observed, conditional, refuted, or unresolved. For an unresolved route, give the smallest remaining proof obligation or discriminating test. Close with one recommendation for the paper's direction and the evidence that would change it.

[pdf]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/output/pdf/Gowdy_Finite_Time_Observability.pdf
[master]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/finite_time_observability.tex
[readout]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/READOUT.md
[corrections]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/CORRECTION_LEDGER.md
[prior]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/PRIOR_ART_AND_REVIEW.md
[matrix]: reviews/2026-11-25_novelty_validity_strategy/CITATION_MATRIX.md
[orbit]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/both_velocities_holonomy.tex
[metric]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/full_metric_holonomy.tex
[replay]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/REPRODUCE.md
[rank]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/rank_two.tex
[circle]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/circular_loop_collision.tex
[ernst]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/exact_rank_two_recovery.tex
[tail]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/low_velocity_tail.tex
[modulus]: research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/identical_velocity_different_modulus.tex
[deps]: reproduction/requirements.txt
[inputs]: reproduction/README.md
