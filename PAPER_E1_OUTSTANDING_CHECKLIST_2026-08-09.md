# PAPER E1 — OUTSTANDING ITEMS CHECKLIST
**Created:** 2026-08-09
**Manuscript:** `PAPER_E1_MANUSCRIPT_DRAFT_01_2026-08-09.md` (Draft 0.2, post major-referee revision)
**Target:** *Classical and Quantum Gravity*, regular research article
**Science status:** frozen except item §0. Evidence chain 14/14.

**STATUS TAGS:** `[OPEN]` `[BLOCKED-ON-x]` `[DECIDED]` `[DONE]` `[KILLED]`

---

## §0. GATE ITEM — run before any text pass

The only remaining computation. Its outcome changes prose in four places, so
running it first avoids writing the same paragraphs twice.

- [ ] `[OPEN]` **PRODUCTION-SIDE MODULUS ESTIMATE.**
  Compute `c_hat_zero_remainder = a_f - q_hat(0) * tau_f` on the **production**
  holdout (not the analytic benchmark), using extracted `q_hat(0)` and the
  collapse-side anchor `a_f`, at depths 16 and 24, resolutions N = 256/512/1024.

  **Why:** the existing `E1-CLAMBDA-ESTIMATOR` run uses analytic inputs
  (`q0`, `anchor`, `remainder` are bit-identical across all three N in the JSON).
  It validates the estimator *formula*. It does not yet show recovery from
  numerically extracted quantities.

  **No known truth available** — score on depth-consistency instead:
  - GATE A: depth-16 and depth-24 estimates agree.
  - GATE B: the gap decreases at second order under refinement.

  **Pre-registered error budget** (freeze before running):
  - neglected remainder ~1e-13 at depth 16, ~1e-20 at depth 24 → negligible;
  - `delta_q_hat * tau_f` with holdout velocity error ~1e-10 → ~1e-9;
  - `delta_a_f` from second-order collapse-side quadrature → dominant term.
  Therefore the depth gap should be set by `delta_a_f` and converge at 2nd order.

  **OUTCOME FORK:**
  - PASS → the production campaign *recovers* the modulus. Delete the
    "imports a finite-depth representative rather than recovering its singular
    limit" qualification from Abstract, §6.2, §7, §8. This closes the last
    circularity objection (that the datum proven necessary is the datum
    supplied by hand).
  - FAIL → change nothing. Qualification stands as currently written; report
    the attempt and the measured obstruction in §6.2.

  Preserve the machine output either way.

---

## §1. TEXT AND CONSISTENCY (no analysis; ~1 hr total)

- [ ] `[OPEN]` **Remainder decay rate sentence, §6.3.** The estimator errors
  6.1277e-8 / 3.3291e-11 / 1.4544e-14 give ratios 1841 and 2289 across
  Δτ_f = 4, against 1987 and 2236 predicted by a τ·e^(−2τ) remainder. State
  this — it is an independent confirmation of the assumed two-term asymptotic
  form, currently unreported.

- [ ] `[OPEN]` **Reframe the modulus-estimator claim.** Present it as an
  analytic-input validation of the estimator formula with the neglected
  remainder falling below roundoff by depth 16 — *not* as a refinement result.
  Say explicitly that the inputs are analytic, so the archived resolution
  columns in `E1_CLAMBDA_ESTIMATOR_RESULT.json` do not invite a referee to ask
  what the refinement was testing.

- [ ] `[OPEN]` **Route-relabel sweep.** `λ_rec = q_hat·τ_f + λ̂_∞` is now the
  primary protocol. Grep for residual "substitution" / "variant" / "bounded
  route" language describing it. Confirm the finite-slice quadrature is
  labelled consistently as the error-isolation arm everywhere (§2.2, §5, §6.1,
  §6.3, §8).

- [ ] `[OPEN]` **Jacobian conditioning explanation, §6.2.** Shape is worse
  conditioned than full (121.8 / 564.6 vs 44.9 / 30.7) — counterintuitive as
  printed. If the mechanism is that removing `c_λ` and the scale direction
  leaves the retained directions more nearly degenerate, say so: it is
  independent corroboration of the §6.2 finding, currently reading as an
  unexplained oddity.

- [ ] `[OPEN]` **Verify the eight parameter-vector entries are listed** by name
  in §6.2. Rank-6 / rank-7 claims are unauditable without them.

- [ ] `[OPEN]` **Confirm the coarse-grid momentum note survives** the revision:
  1.1162e-3 exceeds the 1e-3 gate but the gate applies to the fine grid. If
  that clause was lost, a referee reads a failed gate reported as a pass.

- [ ] `[OPEN]` **Internal-vocabulary strip.** No test IDs in manuscript-facing
  text (`E1-CLAMBDA-ESTIMATOR`, `E1-FG1-...`, gate codes). They belong in the
  archive only. NOTE: the draft ledger line "Add the exceptional known-answer
  control **with its precise test ID**" would reintroduce exactly this — that
  item is now scoped out to future work, so delete the line rather than
  carrying it forward.

---

## §2. APPENDICES

- [ ] `[OPEN]` **Appendix A — homogeneous controls.** State variables, clock
  interpolation, mismatch vector, analytic derivations, and the
  Schwarzschild/KS-vs-stiff-FLRW adverse-control termination criterion.
  Collection, not computation; §3 already points at it.

- [ ] `[OPEN]` **Appendix B — reduced equations from the normalized action.**
  Derive the §4 system so the conventions (8πG = 1, symmetric scalar
  normalization) are auditable.

---

## §3. FIGURES (5)

- [ ] `[OPEN]` **Fig. 1 — protocol schematic.** The
  `C₋ → E → B_adm → R → C₊ → (O,T)` chain with the discarded-history gap and
  the parity map marked.
- [ ] `[OPEN]` **Fig. 2 — production convergence.** Promoted route at both
  depths (3.7478e-5 / 9.3715e-6 / 2.3432e-6), with the error-isolation arm
  overlaid to show the metric-route contribution and the depth-24 growth it
  explains.
- [ ] `[OPEN]` **Fig. 3 — ablation.** Shape layer flat under `c_λ` removal;
  absolute layer at 29 / 106 and non-decreasing under refinement. This is the
  paper's headline finding and currently has no figure.
- [ ] `[OPEN]` **Fig. 4 — exact benchmark.** Field, λ, and momentum errors
  through N = 256/512/1024 with observed orders 1.9988–2.0000 annotated.
- [ ] `[OPEN]` **Fig. 5 — modulus estimator.** Error vs depth on log axes with
  an e^(−2τ) reference slope. Add the production-side points if §0 passes.

---

## §4. BIBLIOGRAPHY

- [ ] `[OPEN]` Verify all 14 entries against published versions (volume, page,
  year, DOI).
- [ ] `[OPEN]` **Ref [6] specifically** (arXiv:1011.0270) — confirm it is
  distinct from [5] and supply full journal details. It carried a
  `[DETAILS TO VERIFY]` placeholder through Draft 0.1.

---

## §5. RELEASE AND ARCHIVE  *(shared operation — see §7)*

- [ ] `[OPEN]` Curated repository: source, frozen specifications, machine
  outputs, adjudication record, independent 14/14 result-chain verifier.
- [ ] `[OPEN]` Environment lock file, checksum manifest, one-command
  clean-environment reproduction entry point.
- [ ] `[OPEN]` **Preserve adverse and superseded evidence explicitly:**
  - the original malformed strict-decrease gate output and its adjudication;
  - the endpoint-copy baseline estimator;
  - the finite-slice metric quadrature arm;
  - the KS/FLRW invariant-mismatch termination.
- [ ] `[OPEN]` Mint Zenodo DOI; replace the `[REPOSITORY AND DOI TO BE
  INSERTED]` placeholder in *Data and code availability*.
  **Do not submit with a live placeholder** — it is a return-without-review
  trigger on its own.

---

## §6. FRONT MATTER AND SUBMISSION

- [ ] `[OPEN]` Author block: name, no institutional affiliation, corresponding
  author, ORCID 0009-0004-1153-2464.
- [ ] `[OPEN]` Acknowledgments.
- [ ] `[OPEN]` Convert to CQG class (iopart); check figure formats and length.
- [ ] `[OPEN]` Data-availability statement per current CQG policy (distinct
  from the free-text §11 block).
- [ ] `[OPEN]` Cover letter. Lead with the *method* contribution and the
  shape-vs-marked-geometry separation; do not lead with "completeness test,"
  which reads as a numerical-methods note to a screening editor.
- [ ] `[OPEN]` Suggested referees + oppositions. Natural pool is Beyer,
  LeFloch, Isenberg, Ames, Rendall, Ringström — note these are also the authors
  of refs [1]–[9], so the suggestion list and the likely-referee list coincide.
- [ ] `[OPEN]` **arXiv endorsement, gr-qc.** Required as an unaffiliated
  author. Cross-program dependency: the halo-spheres ledger carries the same
  ask for astro-ph.GA (natural endorsers Sellwood / Spekkens). Arrange both.

---

## §7. CROSS-PROGRAM NOTES

**BATCHING.** §5 is the same release-engineering operation now pending on four
manuscripts: E1, the Palatini methods paper, the halo-spheres companion
(data-availability block is that paper's single hard blocker), and the junction
submission package. Build the checksum/lock/entry-point template once and
instantiate it four times.

**QUEUE POSITION.** The framework checkpoint (§15.1, CF-2026-08-08-PUB1) still
lists Line E as "scientifically promising but not manuscript-ready; requires
its own fidelity closure and a coherent bounded claim." Both conditions are now
met. E1 currently outranks the Palatini paper on manuscript maturity and
acceptance odds; the checkpoint should be updated so the queue reflects that.

**COMPETING CLOCK.** Only the halo-spheres companion has an external clock
running (AAS78771 at ApJL, desk-passed 2026-07-06, referee report plausibly
imminent). Nothing in this checklist is time-critical; that one is.

---

## LOG

2026-08-09  Created from the Draft 0.2 review. Science frozen except §0.
            Draft-ledger items carried forward and expanded; the
            "exceptional known-answer control" line retired to future work
            per the §8 limitation rewrite.
