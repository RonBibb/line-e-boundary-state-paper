# Research direction after the Gowdy audit

**Lead with finite-time observability.** The strongest new result is a classification of what can be recovered stably from the same observations. The nonlinear simulation and corrected Einstein–scalar benchmark give the theorem a concrete motivation and a boundary for its application.

## Mathematical result

For the polarized Gowdy field with full periodic spatial snapshots at areal times $a,b$ and $b-a=\pi p/q$ in lowest terms, the reconstruction is injective outside a locally finite set of choices of $a$. The target determines the stability cost.

| Target recovered from the same two fields | Sharp loss relative to observation Sobolev index |
|---|---:|
| Field on $a+(\pi/q)\mathbb Z$ at positive time | 0 derivatives |
| Field on any other fixed positive-time slice | 1 derivative |
| Unscaled Cauchy velocity at the reference slice | 2 derivatives, one of which is the usual wave-energy scaling |
| Limiting polarized Gowdy velocity | 1.5 derivatives |
| Linearized marked expansion asymmetry at a lattice slice | 2 derivatives |

The inverse-square term shifts the nominally invisible Fourier modes by order $1/k$. This revives uniqueness. At lattice targets the same small factor occurs in the numerator, which permits stable recovery of the field at those times. A time error of order $1/K^2$ is critical at spectral bandwidth $K$. For selected observation dates a true kernel exists. An exact zero-snapshot mode can satisfy the periodic vacuum constraint, and its two marked histories can have the same metric anchor.

A second result gives a constructive alternative. Divide Fourier modes into dyadic bands. Store each mode once at the reference time and once at its band's time. Uniform velocity recovery in the natural energy normalization uses two coefficients per mode and order $\log K$ distinct observation times. A simultaneous-approximation argument proves that order is necessary for the stated per-mode maximum-norm observation guarantee. This addresses a storage budget for displacement-only data. If full Cauchy velocity can be saved directly, conventional Cauchy logging remains the simplest way to store the full state.

The same finite-time proof works with a time-only potential under the stated nondegeneracy hypotheses. For the periodic Euler–Poisson–Darboux family with damping $\alpha/t$, $0<\alpha<2$, the limiting flux costs $2-\alpha/2$ derivatives from the fixed pair and $1-\alpha/2$ with the band schedule. This supplies an application beyond a single cosmological field equation.

## Closest literature and remaining priority question

- Szíjártó and Hegedűs, Theorem 1, and Vest, Theorems 4.1 and 5.1, already give sharp or near-sharp arithmetic estimates for fixed-time observations of vibrating systems. Constant loading already lifts rational-gap degeneracy with a one-derivative cost in the wave-energy normalization. The positive inverse-square potential, target-specific sharp bounds, and perturbation/timing statement are the parts to compare directly. [Klein–Gordon observation](https://doi.org/10.14232/ejqtde.2012.1.7) · [Vest](https://www.math.u-szeged.hu/ejqtde/p2411.pdf)
- Christensen et al., Theorems 2.6–2.7, already recover a lattice of target snapshots from two observations for the free Euclidean wave equation. In the current problem the rational lattice comes from the integer spatial spectrum and the cancellation of perturbed resonances. It can include a finer spacing than $b-a$. The unstable complement, the singular endpoint, and the effect of timing errors carry the mathematical distinction. [The snapshot problem for the wave equation](https://arxiv.org/pdf/2308.12208)
- Amorim and LeFloch (2008) establish sharp regularity for the periodic singular Euler–Poisson–Darboux Cauchy problem. Li gives the Kasner scattering isomorphism in Theorem 1.3. These establish important relationships between asymptotic coefficients and solutions. The extra loss here arises before that step, from inversion of two displaced-field samples. [Periodic EPD estimates](https://ems.press/content/serial-article-files/44605) · [Kasner scattering](https://arxiv.org/pdf/2401.08437)
- Budzinskiy, §1 Eqs. (1) and (3), describes the classical Bessel cross-products and annular root asymptotics. Aldroubi et al., §1, set up sampling a known evolution at chosen times from partial spatial data. Bessel root computations and the general principle of distributed sensing are established. The specific two-coefficient schedule and matching logarithmic lower bound are candidates for specialist confirmation. [Bessel cross-products](https://arxiv.org/pdf/1909.00293) · [Dynamical sampling](https://arxiv.org/pdf/1409.8333)
- Amorim, Bernardi and LeFloch, §§3–4.2.5, give the closest bidirectional Gowdy computation. Cacciatori and Ursino, Lemma 5.2 in v3, use the classical conserved Fourier–Bessel coefficients to describe polarized asymptotics in a different geometric problem. Niedermaier and Sedighi Jafari's 2026 work treats Dirac observables and velocity domination. None of the passages examined establishes the same two-snapshot target and timing estimates. This is a targeted novelty check pending expert review. [Bidirectional Gowdy](https://arxiv.org/pdf/0811.1266) · [Conserved polarized coefficients](https://arxiv.org/pdf/2308.00080v3) · [Dirac observables](https://doi.org/10.3842/SIGMA.2026.097)

## Evidence and suggested feedback

The first reader should check the arithmetic bounds, the data norm in the logarithmic lower bound, and whether the same schedule is known in dynamical sampling. An applied mathematics or PDE/control specialist is the appropriate reader. A mathematical relativist can independently verify the corrected Einstein–scalar reduction and the role of marking. A university discussion can center on one figure and the difference between a stable earlier slice and an ill-conditioned velocity from the same data.

A successor is worth developing if those readers confirm the finite-snapshot and storage results add something to the established free and loaded-wave results. The claims on the singular flux, metric signs and decoder comparison should remain explicit. The current scope is polarized analytic results and finite-window nonlinear evidence. A general nonlinear inverse estimate near a true Gowdy spike would be a separate research contribution.

[Working article source](snapshot_observability.tex) · [Corrections note](CORRECTIONS.md) · [Modal computation](snapshot_checks.json) · [Corrected geometry computation](geometry_checks.json)

## Nonlinear finding that changes the earlier paper

The original depth-18 production tape reproduces the mask score of 0.822 at the finest grid. Matching the saved window while changing the restart stencil to fourth order reduces that score to 0.00298 on identical observations. The corresponding score with the numerical endpoint velocity supplied is 0.00302. Three-grid improvement is visible in the fourth-order restart against the historical reference. The largest complementary-mask error drops from 3.034 to about 0.00295. This supersedes the interpretation of intrinsic information loss for that particular tested output. The inverse uses the specified PDE and finite window, with no exact Lim template, hidden location, or endpoint velocity.

The separate fourth-order forward tape on a different, constrained periodic IVP gives a score drop from 0.23614 to 0.002405 with the second-order inverse. The privileged comparison is 0.002479. The two tested noise directions give approximate comparison gains of 5.17 and 8.36 and leave a fine-grid momentum residual of $4.87\times10^{-4}$ on the underlying reference evolution. Results at depth 12, reflected components, and local synthetic fits should remain bounded to their observed cases.

[Primary nonlinear runs](nonlinear_checks.json) · [Distinct-IVP and noise test](nonlinear_robustness.json)
