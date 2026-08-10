# Paper E1 Novelty and Prior-Art Audit

**Date:** 2026-08-09  
**Evidence cutoff:** searches and sources reviewed through 2026-08-09  
**Status:** initial admission audit; bibliography-level verification still required  
**Scope:** Paper E1 as frozen in `PAPER_E1_ADMISSION_AND_SCOPE_FREEZE_2026-08-09.md`

## 1. Executive ruling

Paper E1 has a plausible publishable contribution, but its novelty is narrower
than the existence of collapse/bang asymptotic equivalence.

The established literature already contains:

1. Einstein–scalar quiescent singularities described by freely specifiable
   asymptotic data;
2. Fuchsian singular initial-value formulations for Gowdy and related
   spacetimes;
3. numerical algorithms that construct Gowdy solutions from data prescribed
   at the singularity;
4. geometric singular initial-data formulations with uniqueness of
   development; and
5. formal contracting-to-expanding singularity scattering maps.

Consequently, Paper E1 must not claim novelty for time reversal, AVTD
asymptotics, reconstruction from singular data, or the abstract idea of
relating data on two sides of a singular hypersurface.

The potentially novel contribution is the following combined validation
protocol:

> Starting from an ordinary finite-time contracting Cauchy evolution, extract
> a deliberately reduced candidate boundary representation, discard the bulk
> history, independently initialize an expanding branch from only the retained
> data, and test representational completeness through normalized geometric
> observables evaluated on multiple invariant clocks, with positive, negative,
> convergence, constraint, and information-loss controls.

This is best presented as an **a posteriori extraction–reconstruction and
boundary-data-completeness test**, rather than as a new singular initial-value
theory or a physical bounce law.

## 2. Nearest prior-art families

### 2.1 Quiescent Einstein–scalar singularities

Andersson and Rendall proved the existence of broad analytic families of
Einstein–scalar or stiff-fluid spacetimes with quiescent asymptotics and the
expected number of free functions. Near the singularity, the evolution at
different spatial points decouples in the stated Gaussian description.

**Implication for E1:** neither quiescent Einstein–scalar behavior nor a
functional asymptotic parametrization is new.

Primary source:

- Lars Andersson and Alan D. Rendall, *Quiescent cosmological singularities*,
  [arXiv:gr-qc/0001047](https://arxiv.org/abs/gr-qc/0001047).

### 2.2 Gowdy AVTD asymptotics and singular data

Kichenassamy and Rendall constructed Gowdy solutions using Fuchsian reduction,
with precise Kasner-like asymptotics and freely specifiable functions. Later
work extended AVTD existence and gauge control in smooth Gowdy families.

**Implication for E1:** the four-function-type asymptotic representation and
the idea that such functions determine smooth Gowdy solutions are established
foundations. E1 must identify exactly how its Einstein–scalar polarized data
specialization and scored observables relate to these formulations.

Primary sources:

- Satyanad Kichenassamy and Alan D. Rendall, *Analytic description of
  singularities in Gowdy spacetimes*,
  [arXiv:1709.09710](https://arxiv.org/abs/1709.09710).
- Ellery Ames, Florian Beyer, James Isenberg, and Philippe G. LeFloch, *A
  class of solutions to the Einstein equations with AVTD behavior in
  generalized wave gauges*,
  [arXiv:1602.03018](https://arxiv.org/abs/1602.03018).

### 2.3 Numerical Fuchsian reconstruction

Beyer and LeFloch established a second-order hyperbolic Fuchsian framework and
introduced a numerical algorithm for approximating Gowdy solutions from
singular asymptotic data. This is the closest direct methods competitor.

**Implication for E1:** “we numerically construct a Gowdy solution from
asymptotic data” is not a novel methods claim. E1 must distinguish its workflow
as a test of data extracted from a separate finite-time Cauchy evolution, its
independent opposite-orientation reconstruction, its invariant comparison
observables and clocks, and its explicit test of data loss.

Primary sources:

- Florian Beyer and Philippe G. LeFloch, *Second-order hyperbolic Fuchsian
  systems. Gowdy spacetimes and the Fuchsian numerical algorithm*,
  [arXiv:1006.2525](https://arxiv.org/abs/1006.2525).
- Florian Beyer and Philippe G. LeFloch, *Second-order hyperbolic Fuchsian
  systems and applications*,
  [arXiv:1011.0270](https://arxiv.org/abs/1011.0270).

### 2.4 Cauchy-data/asymptotic-data correspondence

Ames, Beyer, and Isenberg established, for a linearized Einstein–scalar
subsystem on contracting backgrounds, a homeomorphism between ordinary Cauchy
data and asymptotic data.

**Implication for E1:** the general concept of extracting complete asymptotic
data from Cauchy data already has a rigorous precedent. E1's distinction must
be its nonlinear declared test systems, operational finite-resolution
protocol, opposite-branch reconstruction, and invariant scoring—not the bare
concept of a Cauchy/asymptotic correspondence.

Primary source:

- Ellery Ames, Florian Beyer, and James Isenberg, *Contracting asymptotics of
  the linearized lapse-scalar field sub-system of the Einstein-scalar field
  equations*,
  [arXiv:1904.02854](https://arxiv.org/abs/1904.02854).

### 2.5 Geometric singular initial data and unique development

Recent work proves a geometric singular initial-value formulation for the
Einstein–nonlinear-scalar system, without symmetry or analyticity under its
declared hypotheses, including uniqueness of development from singular data.

**Implication for E1:** E1 should not imply that it establishes the first
geometric or complete singular-data formulation. Its reduced state vector is
an operational diagnostic representation, and its relationship to geometric
singular data must be stated carefully.

Primary source:

- Andrés Franco-Grisales, *Developments of initial data on big bang
  singularities for the Einstein-nonlinear scalar field equations*,
  [arXiv:2409.17065](https://arxiv.org/abs/2409.17065).

### 2.6 Two-sided singularity scattering maps

Le Floch, LeFloch, and Veneziano define data spaces on the contracting and
expanding sides of quiescent singularities and classify covariant ultralocal
singularity scattering maps under their assumptions.

**Implication for E1:** the idea of a map relating contracting and expanding
singularity data is explicitly established. E1 must distinguish its identity
map/time-parity comparison as a diagnostic hypothesis being tested, not a
selected junction law. E1 also supplies no physical scattering prescription.

Primary sources:

- Bruno Le Floch, Philippe G. LeFloch, and Gabriele Veneziano, *Cyclic
  spacetimes through singularity scattering maps. The laws of quiescent
  bounces*, [arXiv:2005.11324](https://arxiv.org/abs/2005.11324).
- Philippe G. LeFloch, *Singularity scattering laws for bouncing
  cosmologies: a brief overview*,
  [arXiv:2106.07958](https://arxiv.org/abs/2106.07958).

### 2.7 Spike and scale-invariant state-space literature

Gowdy spikes and the distinction between true and false spikes are established
topics, including descriptions in Hubble-normalized scale-invariant variables.

**Implication for E1:** spike formation and normalized variables are not new.
The bounded E1 contribution is the known-answer demonstration that a declared
restricted representation can lose localized exceptional information. The
production spike's nature remains outside E1.

Primary source:

- Lars Andersson, Henk van Elst, and Claes Uggla, *Gowdy phenomenology in
  scale-invariant variables*,
  [arXiv:gr-qc/0310127](https://arxiv.org/abs/gr-qc/0310127).

## 3. Claim-by-claim novelty ruling

| Candidate claim | Ruling | Reason |
|---|---|---|
| Einstein–scalar contraction can have quiescent/AVTD asymptotics | Established | Rigorous asymptotic literature predates E1. |
| Smooth Gowdy solutions can be parametrized by singular asymptotic functions | Established | Central Fuchsian result. |
| Gowdy solutions can be numerically generated from singular data | Established | Beyer–LeFloch is a direct precedent. |
| Contracting and expanding data can be related across a singular boundary | Established as a formal framework | Singularity scattering maps already treat two-sided data. |
| Time reversal changes canonical parity while preserving intrinsic data | Established/basic | Consequence of reversible classical equations and canonical definitions. |
| Homogeneous stiff-scalar equivalence is stratified by normalized shear | Likely a useful synthesis, not yet shown novel | May be an elementary consequence of Jacobs/Kasner-like dynamics; requires targeted audit. |
| Candidate boundary data extracted from an ordinary contraction can be tested by independent opposite-branch reconstruction | Potentially novel combination | No identical protocol located in this initial audit. |
| Reconstruction can be scored through a declared invariant state vector on two independent invariant clocks | Potentially novel implementation | Requires proof that every scored quantity and clock is genuinely invariant and nonredundant. |
| Positive and negative controls can distinguish compatibility rather than mere reversibility | Potentially novel validation design | Strength depends on using identical machinery and on the negative control being appropriately posed. |
| Exceptional known-answer controls diagnose representational incompleteness | Potentially novel application | Spike information loss itself is known; its role as a controlled failure of the completeness test may be distinctive. |
| The tested smooth polarized family reconstructs with pointwise second-order convergence toward zero residual | New numerical result if independently verified | Must remain restricted to the declared data family and discretization. |

## 4. The novelty statement E1 should use

Recommended manuscript language:

> Existing Fuchsian and geometric singular-initial-value theories establish
> that suitable asymptotic data can parameterize quiescent Einstein and Gowdy
> solutions, and numerical Fuchsian algorithms construct solutions from such
> prescribed data. We address the complementary validation problem: given a
> finite-time contracting Cauchy evolution and a proposed reduced boundary
> representation, does that representation retain enough information to
> reconstruct the asymptotic invariant state of an independently initialized
> expanding branch? We formulate an extraction–reconstruction protocol,
> evaluate it on multiple invariant clocks, and subject it to compatible,
> incompatible, convergence, constraint, and exceptional-information
> controls.

This statement is safer and more differentiated than claiming a new
collapse–bang correspondence.

## 5. Required differentiation from the closest competitors

Before submission, the manuscript must include a comparison table covering at
least:

| Feature | Fuchsian existence theory | Numerical Fuchsian algorithm | Singularity scattering maps | Paper E1 |
|---|---:|---:|---:|---:|
| Data prescribed at singularity | Yes | Yes | Yes, on each side | Candidate data extracted from finite-time contraction |
| Ordinary Cauchy-to-asymptotic extraction tested | Not generally the numerical focus | Not the central published test | Assumed data space | Yes |
| Independent opposite-orientation branch | Not the central validation target | Not the central validation target | Formal map framework | Yes, for the declared reconstruction |
| Multiple invariant clocks | Not identified in this audit | Not identified in this audit | Geometric data rather than this scoring design | Yes, subject to audit |
| Aggregate and pointwise completeness residual | Not this form | Numerical solution error, not this declared two-branch residual | Not a numerical completeness statistic | Yes |
| Matched positive and negative controls | Not this protocol | Not identified in this audit | Classification rather than numerical validation | Yes |
| Exceptional-data failure control | Spike theory exists separately | Not identified as this protocol | Not the E1 numerical test | Yes, in bounded known-answer scope |
| Physical continuation selected | No | No | Prescribed/classified maps | No |

The phrases “not the central target” and “not identified in this audit” must
not be converted into absolute priority claims until full-text and citation
searches are complete.

## 6. Publication-risk assessment

### Risk A — method viewed as a repackaged Fuchsian numerical algorithm

**Severity:** high.  
**Mitigation:** make the extraction from a separate Cauchy evolution, data
discard, independent reconstruction, invariant scoring, and failure controls
mathematically explicit. Compare algorithms directly rather than merely citing
them.

### Risk B — result viewed as tautological time reversal

**Severity:** high.  
**Mitigation:** prove branch independence, run the same protocol on an
incompatible negative control, and show controlled failure when boundary data
are incomplete.

### Risk C — invariant language exceeds what the variables support

**Severity:** medium to high.  
**Mitigation:** separate genuine spacetime scalars and operational clocks from
gauge-fixed, frame, topology, and constraint diagnostics. Rename the method if
the full state vector is not invariant.

### Risk D — one tested polarized family is too narrow for a standalone paper

**Severity:** medium.  
**Mitigation:** emphasize validation methodology, exact/adverse controls, two
depths, three resolutions, component-level closure, and representational
failure controls. Do not broaden the sector-level claim.

### Risk E — physical interpretation overshadows the actual result

**Severity:** high if “collapse–bang” leads the title or abstract.  
**Mitigation:** retain the conservative title and state on page one that no
crossing, bounce, or continuation law is derived.

## 7. Journal-fit ruling

The provisional *Classical and Quantum Gravity* target remains reasonable only
if the method is formally distinguished from the numerical Fuchsian algorithm
and every headline result survives independent verification. If the paper
reduces to a numerical confirmation of known polarized AVTD reconstruction,
*General Relativity and Gravitation* or a substantially narrower technical
article would be the more realistic target.

## 8. Step-2 disposition

**Provisional result:** `NOVELTY-PASS-WITH-NARROWING`.

Paper E1 is admitted to continue, with novelty frozen provisionally as:

\[
\boxed{
\text{Cauchy extraction}
\rightarrow
\text{reduced boundary representation}
\rightarrow
\text{independent opposite-branch reconstruction}
\rightarrow
\text{invariant completeness audit}
}
\]

The following are explicitly retired as novelty claims:

- existence of quiescent Einstein–scalar singularities;
- AVTD behavior in smooth Gowdy models;
- construction of solutions from singular data;
- the general concept of contracting/expanding singularity maps; and
- the existence or phenomenology of Gowdy spikes.

The next required work is the equation-and-convention audit, followed by a
formal comparison of E1's actual reconstruction algorithm against the
Beyer–LeFloch Fuchsian numerical algorithm.
