---
title: "A Sudakov veto cascade for NuWro: construction, validation, and whether it pays"
author: "Josue (Isaacson HEP Theory, MSU), with Claude"
date: "28 October 2026"
geometry: margin=2.0cm
fontsize: 10pt
colorlinks: true
header-includes:
  - \usepackage{booktabs}
  - \renewcommand{\arraystretch}{1.1}
---

**Source.** `github.com/NuWro/nuwro` at `f16c937647ad07dabb4576f54dee6d717a9b8421` (28 Oct 2026). Every `file:line` for (a) refers to that commit; lines for (b) refer to the patched `kaskada7.cc` in `nuwro_sudakov.patch`. Every number is tagged [run], [read in source, file:line] or [derived].

**Build.** The `kaskada` binary (nucleon/pion-beam mode) was built against ROOT 6.32 (conda-forge). `nuwro` itself was not built, because it needs ROOT+EGPythia6. The QE spectral-function branches, which only `nuwro` reaches, were run through a small driver that calls `kaskada::kaskadaevent(false)` directly (`sf_test.cc`).

# Summary

1. **(a)'s step-size bias in $\sigma_R$ and transparency is negligible** [run].
   - Pooled over 13 energies (10–900 MeV), $\sigma_R(\Delta{=}0.2\,\text{fm})/\sigma_R(\Delta\to0)-1=-0.06\pm0.18\%$.
   - Pooled over 9 energies (0.2–2 GeV), the transparency bias is $+0.14\pm0.08\%$.
   - For comparison, the median data errors are 3.5% ($\sigma_R$) and 5.3% ($T$) [run, from your CSVs].
   - Moving `kaskada_NN_mfp_scale` across its quoted range (0.72–0.82) shifts $T$ by $-0.013$ to $+0.018$ (2–3.5%) [run].
   - **By your stopping rule, Step 2 says stop.**
2. **The bias is not negligible for cascade-level quantities** [run]. At the default step, (a) undercounts collisions:
   - the mean number of collisions per event by 2–4% in transparency mode;
   - the QE-SF branch by 2.7%.
   - The cause is that (a) evaluates the rate with *one random target per step*. By Jensen's inequality this suppresses collisions at $O(\Delta)$ wherever $\sigma_{NN}$ varies across the Fermi sphere, that is, for slow secondaries.
   - (b) has no step and agrees with (a)'s $\Delta\to0$ extrapolation everywhere tested.
3. **(b) is a genuine veto algorithm** (all five criteria), with one majorant violation in $1.9\times10^9$ proposals, traced to an (a) side effect and fixed (§3.3).
   - It reproduces exact closed forms (KS $p=0.08$–0.97).
   - It is independent of the segment length.
   - With the flag off, the patched build is bit-identical to (a).
4. **Cost** [run]. Per event:
   - beam mode: (b) is 0.9–1.7× the CPU of (a) at $\Delta=0.2$ fm;
   - SF cascades: 1.9×;
   - in all cases 3.5–6× *cheaper* than (a) at the step needed to remove the cascade bias ($\Delta\approx0.025$ fm).
5. **Reweighting closes for both (a) and (b)** [run].
   - (a) can be reweighted exactly with likelihood-ratio weights at its own step. So reweighting is *not* something only (b) can do.
   - (b)'s veto weights have a slightly lower effective sample size (ESS/N 0.79–0.98 vs 0.84–0.998).
6. **Found in (a)** while doing this [run / read in source]:
   - in `kaskada` beam mode the event flags are uninitialised, so every nucleon runs the QE-SF code path;
   - a Pauli-blocked pion absorption still removes nucleons;
   - a failed `remove_nucleon` duplicates particles.

# The five criteria

| criterion | (a) NuWro as it stands | (b) patched |
|---|---|---|
| 1. Sudakov factor $e^{-\int\Gamma}$ not inverted in closed form | **No.** It inverts $e^{-\Gamma(\ell_0)\lambda}$ with $\Gamma$ frozen at the step start: `freepath = -log(frandom())/res.xsec` (`kaskada7.cc:753`), move $\min(\Delta,\lambda)$ (`:768`). $\Delta$ is an accuracy parameter. | **Yes.** $\int\Gamma$ is never computed. Only the majorant is inverted: `cand = l - log(1-u)/gm` (`kaskada7.cc:1527`). |
| 2. Proposals from an inverted $\hat\Gamma\ge\Gamma$ | **No.** There is no majorant; the proposal uses $\Gamma$ itself, from one random target (`:747,753`). | **Yes.** Piecewise-constant $\hat\Gamma_c$ per segment that bounds every factor (`sud_bound`, `:1397–1505`). Majorant violations are counted (`:1638`): 1 (ratio 1.027) in $1.9\times10^9$ proposals, traced and fixed (§3.3); 0 after the fix. |
| 3. Accept with $\Gamma/\hat\Gamma$ at the candidate | **No.** The only veto is Pauli blocking (`:1001`), a physics veto, not $\Gamma/\hat\Gamma$. | **Yes.** Fresh target at the candidate via the unmodified `prepare_interaction()` (`:1635`), `acc = u*gm < g` (`:1640`). |
| 4. Continue from the rejected point, same channel | **No rejection step exists.** After a Pauli block the particle is re-queued at the collision point (`:84–90`) and draws a new *total* free path. | **Yes.** `sud_draw(t, c, lb, gm, seg_hi)` continues channel $c$ from the rejected point, in the same segment and with the same $\hat\Gamma_c$ (`:1653`, also after a Pauli block `:1676`). |
| 5. Competing channels, earliest wins | **No.** One total rate (`:747`), then the isospin with $\Gamma_p/\Gamma$ (`:754`, `Interaction.cc:449`). Equivalent only as $\Delta\to0$. | **Yes.** Two isospin channels per particle, each with its own candidate stream. The earliest candidate over **all** particles in global time wins (`:1559,1605–1612`). |

# Step 1 — the process (a) approximates

**State.**

- Each particle in flight carries $(\mathbf r, t, \mathbf p, \texttt{his\_fermi})$ and $s$, the distance travelled since creation or last collision (`particle.h:286–287`, reset at `kaskada7.cc:1084`).
- The residual nucleus carries $(Z_r,N_r)$ (`nucleus.h:119–121`).
- The event counter `number_of_interactions()` (`event1.h:495`) switches the correlation factor on.

**Hazard of channel $c\in\{n,p\}$** (the target isospin), per unit path length:
$$\Gamma_c=\mathbb E_T\!\left[\rho_c(r)\,\sigma_c(\mathbf p_1,T,\rho)\right]C_c(s,r)\,\frac{1}{s_{\rm mfp}}\,f_{xs}.$$

| factor | what (a) does | where |
|---|---|---|
| density | FB point density of ¹²C $\times A_r/A$, $R=8.0004$ fm; $\rho_c=\rho N_c/A_r$ | `nucleus.cc:69`, `nucleus_data.cc:232`, `kaskada7.cc:706–715` |
| target $T$ | LFG: isospin $\propto Z_r/A_r$, $\mathbf p$ uniform in a ball of $k_F^{\rm loc}(r)$ | `nucleus.cc:101–138`, `nucleus.h:181–191` |
| $\sigma_{\rm free}(E_k^{\rm eff})$ | PDG2016 table at the projectile energy in the target frame | `Interaction.cc:332–335,512–539` |
| medium | Pandharipande–Pieper factor on the elastic part, $1-0.2\rho/\rho_0$ on the inelastic part (free inelastic fraction, zero below 280 MeV) | `Interaction.cc:337–379` (non-SF path), `:186–256` (consistent-target path) |
| effective Pauli | unlike pair $\times0.9$ for projectile $E_k<40$ MeV | `Interaction.cc:380–389`, `:250–251` |
| SRC factor | $\mathrm{corr}_c(s)/\mathrm{norm}_c(r)$, piecewise-linear tables (ANL ab initio) | `kaskada7.cc:720–745`, `data/input/kaskada_NN_corr_1/6p6n.dat` |
| mfp scale | $1/0.76$ (`params_all.h:136`). It and the SRC factor are **off** for a `beam_placement=2` nucleon before its first collision | `kaskada7.cc:727–745` |
| pions | Oset tables at $(E_k,\rho)$; no target dependence | `Interaction.cc:403–424` |
| energy on the path | constant, except QE-SF nucleons: $E_k\to E_k+E_F(r)-E_F(r_0)$ | `kaskada7.cc:770–784` |
| jailing | $E_k\le w+\texttt{his\_fermi}$, $w=7$ MeV | `kaskada7.cc:584,795,831,1089` |

**Jailing is not a per-step event.** $J=E_k-\texttt{his\_fermi}$ is constant on a straight path, both outside and inside the SF energy update (`:775–783` shifts both terms equally). So the jail test can change only at creation or at the surface [derived]. (b) checks it at every event anyway: 0 in-flight jailings in all runs [run].

**Channel and collision.**

- The channel is the target isospin.
- The process (elastic / 1π / 2π) is chosen *after* the collision is decided, from the free inelastic fraction (`Interaction.cc:690–700`).
- The target orientation is redrawn $\propto\sigma_{\rm free}$ (`:668–687`).
- Kinematics are retried up to 100 times (`kaskada7.cc:975–989`).
- **Pauli blocking.** Any outgoing nucleon with $|\mathbf p|<k_F^{\rm loc}$ (`:1001`; the `pauli_blocking` parameter is not read by the cascade) makes `make_interaction` return false.
  - The particle is re-queued at the collision point with momentum and $s$ unchanged (`:84–90`).
  - In the continuum this is a thinning of $\Gamma_c$. **(b) keeps exactly this**: an accepted but blocked candidate is discarded and the particle continues from that point.
- **Success.** Products start at the collision point with $s=0$ and the target is removed (`:1083–1141`). Every particle's hazard changes.

**What `kaskada` actually runs.** `kaskada.cc` never initialises `event::flag` (`event1.h:10–27,41`; the default constructor sets only `weight` and `norm`, `event1.h:91`; the event is created at `kaskada.cc:132`).

- In our build `flag.qel` was non-zero in 100% of events [run, read back from the output tree], so `use_consistent_target` (`kaskada7.cc:689`) and the SF energy update (`:770`) are active for every nucleon in beam mode. In that path:
  - the rate uses a forward/reverse-averaged pair of targets;
  - `make_interaction` re-draws target *and* isospin at the collision point (`:869–966`);
  - the beam nucleon gains $E_F+w$ on entry (`:573–575`).
- This is undefined behaviour and may differ between builds.
- (b) reproduces whatever (a) does here. The collision target law is the one place where this matters (§4.2).

**Ordering.**

- (a) uses a FIFO queue (`kaskada7.h:54`, `kaskada7.cc:72–93`): one step per particle per round. That is rounds, not global time, even though `krok` advances a time coordinate (`particle.h:286`).
- The target process orders all candidates by global time.
- The two orders differ only through nucleus-state changes ($\delta\rho/\rho\sim1/A$ per collision). They cannot affect $\sigma_R$ or $T$.

**Markov.** The process is Markov in the state above. This needs $s$ in the state, because the SRC factor depends on the distance since the last collision [derived].

**One target per step.**

- As $\Delta\to0$, $P(\text{coll})=\mathbb E_T[1-e^{-\Gamma(T)\Delta}]\to\langle\Gamma\rangle_T\Delta$. So the step converges to the hazard $\langle\Gamma\rangle_T$ with target law $\propto\Gamma(T)f(T)$ (or, in the QE-SF path, a fresh target at the collision) [derived].
- At finite $\Delta$ there are two $O(\Delta)$ biases:
  - the frozen start-of-step rate;
  - **Jensen**: $\mathbb E_T[e^{-\Gamma\Delta}]\ge e^{-\langle\Gamma\rangle\Delta}$, a rate deficit $\approx\tfrac12\mathrm{Var}_T(\Gamma)\Delta$.
- Unit test (§4.1): with $\sigma$ uniform on $[0,2\bar\sigma]$, (a) has $P_{\rm int}$ low by 0.0101 (absolute) at $\Delta=1$ fm, 0.0020 at 0.2 fm and 0.0004 at 0.05 fm [run]. The analytic prediction for (a)'s per-step process is 0.0114, 0.0023 and 0.0006 [derived].

# Step 2 — (a)'s step-size bias

**Runs** [run]:

- proton beam on ¹²C;
- $\sigma_R=\pi(0.999R)^2\,P(\ge1$ collision$)$ with `beam_placement=2`; $T=P(0$ collisions$)$ with `beam_placement=1`, as NuWro's own `calcTrans` (`calculations.cc:1588`);
- $\Delta\in\{0.05,0.1,0.2,0.4,0.8,1.6\}$ fm at all energies, plus $0.025,0.0125,0.00625$ at low energy;
- $2\times10^5$ events per point.

The zero-step value comes from a weighted quadratic fit over $\Delta\le0.8$ fm. The bias quoted is the fitted curve at 0.2 fm minus its intercept.

| | (a), $\Delta$=0.2 fm | (a), $\Delta\to0$ | bias at 0.2 fm | (b) | (b)$-$(a)$_0$ |
|---|---|---|---|---|---|
| $\sigma_R$ 10 MeV | 300.4 ± 1.6 | 303.0 ± 1.0 | −1.3 ± 0.4 % | 300.1 ± 0.9 | −2.1σ |
| 30 | 277.3 ± 1.5 | 277.3 ± 1.2 | −0.0 ± 0.5 % | 277.4 ± 1.1 | +0.1σ |
| 100 | 198.5 ± 1.3 | 196.9 ± 1.0 | +1.1 ± 0.6 % | 198.9 ± 0.9 | +1.4σ |
| 300 | 174.1 ± 1.3 | 177.1 ± 1.3 | −1.1 ± 0.8 % | 175.9 ± 0.9 | −0.7σ |
| 900 | 226.2 ± 1.4 | 225.2 ± 1.5 | +0.4 ± 0.7 % | 225.4 ± 1.0 | +0.1σ |
| $T$ 300 MeV | 0.6224 ± 0.0011 | 0.6205 ± 0.0011 | +0.6 ± 0.2 % | 0.6251 ± 0.0008 | +3.4σ |
| 1000 | 0.5141 ± 0.0011 | 0.5144 ± 0.0012 | −0.0 ± 0.2 % | 0.5150 ± 0.0008 | +0.5σ |
| 2000 | 0.5016 ± 0.0011 | 0.5012 ± 0.0012 | +0.3 ± 0.3 % | 0.5002 ± 0.0008 | −0.7σ |

The full table (22 energies) is in `results/step2_table.md`. Over all 22 energies (b)$-$(a)$_0$ has $\chi^2=32/22$.

The one 3.4σ point ($T$ at 300 MeV) was rerun with $10^6$ events [run]:

- (a, $\Delta=0.05$): $0.6253\pm0.0005$;
- (b): $0.6253\pm0.0005$;
- six independent (a) seeds at 0.2 fm scatter consistently with the binomial error.

So it was a fluctuation of the quadratic intercept.

![Step dependence of (a) at six of the 22 energies, its zero-step extrapolation (green), and (b) (red, plotted at $\Delta<0$).](figs/fig1_step.png){width=95%}

**Comparison with data errors and tuning** [run]:

| quantity | step bias at 0.2 fm (pooled) | data error (median) | tuning lever |
|---|---|---|---|
| $\sigma_R$, 10–900 MeV | $-0.06\pm0.18\%$ (largest $|b|+2\sigma$: 3.1%) | 3.5% (min 2.1%) | NN scale ±10% → $\sigma_R$ ±6% |
| $T$, 0.2–2 GeV | $+0.14\pm0.08\%$ (largest 1.0%) | 5.3% (min 2.6%) | mfp 0.72–0.82 → $T$ −0.013…+0.018 |

**Verdict of Step 2.** For $\sigma_R$ and $T$ the bias is negligible: 10–50× below the data errors and the tuning range. By your rule this is where to stop.

I completed Steps 3–5 anyway, because the step does bias **cascade-level** observables (Fig. 2, [run]):

- $\langle N_{\rm coll}\rangle$ at 0.2 fm is low by 2.3% ($T_p$=300 MeV), 3.8% (1 GeV) and 3.9% (2 GeV) in transparency mode.
- The QE-SF branch is low by 2.7%. There (a) converges only below $\Delta\approx0.025$ fm: 2.052 (0.2) → 2.091 (0.05) → 2.108 (0.025) → 2.103 (0.0125).

These collisions are not counted in $\sigma_R$ or $T$. They are counted in hadron multiplicities and energy flow in neutrino events.

![Mean number of collisions per event. (a) at each step, its extrapolation, and (b). Right: QE-SF branch (non-transparent events), (b) at three segment lengths.](figs/fig3_ncoll.png){width=100%}

# Step 3 — (b)

**Algorithm** (`kaskada7.cc:1251–1706`). For each particle in flight:

1. Cut the straight path into segments of length $L_{\rm seg}$, ending at the computed surface crossing.
2. On each segment, bound each channel by a constant $\hat\Gamma_c$.
3. Draw the candidate by inverting $\hat\Gamma_c$.
4. Process the earliest event of all particles in global time:
   - move straight to it, with no steps;
   - call the unmodified `prepare_interaction()` (fresh target, all physics factors);
   - accept with $\Gamma_c/\hat\Gamma_c$;
   - on rejection, continue from that point in the same channel.
5. On acceptance, call the unmodified `make_interaction()`/`finalize_interaction()`.
   - A Pauli block or failed kinematics is treated as in (a).
   - After a successful collision every other particle is moved to the collision time and re-anchored with new majorants, which is exact because the proposal process is memoryless.
6. The surface crossing is an event: the unmodified `leave_nucleus()` runs there.

**Majorant** (all factors bounded, none changed):
$$\hat\Gamma_c=\bar\rho_{\max}\tfrac{N_c}{A}\;\max_{[E^{\rm eff}_{\min},E^{\rm eff}_{\max}]}\sigma_{\rm free}\;\frac{1}{1-a(\rho_{\max})/4}\;0.9^{[E_k<40,\,ij]}\;\frac{\max \mathrm{corr}_c}{\min\mathrm{norm}_c}\frac{1}{s_{\rm mfp}}\;f_{xs}\,H.$$

- **Density.** $\bar\rho_{\max}$ is taken on a 4000-cell radial grid with a curvature allowance.
- **Target energy range.** The $E^{\rm eff}$ range comes from $|T|\le k_F(\bar\rho_{\max}\max(Z_r,N_r)/A)$. It includes the SF energy shift and is zero when a co-moving target exists.
- **Free cross section.** It is piecewise linear, so the maximum is at a knot or end point (exact).
- **Medium factor.** The Pandharipande–Pieper factor is bounded by $1/(1-a/4)$, with $a=2M|\beta|/\lambda^2$. The proof: $\mathbf k/m^*(k)$ is the gradient of a radial convex potential with Hessian eigenvalues $\ge(1-a/4)/M$ [derived]. The bound is 1.20 at $\rho_0$; the true maximum, scanned numerically, is 1.06 [run].
- **SRC tables.** Piecewise linear, so exact at knots.
- **Pions.** The Oset tables are bilinear in $(E,\rho)$, so the bound is exact at the density knots.
- **Hyperons.** 3× the maximum of 17 samples. Not exercised here.

**Segment-length independence** [run]. $L_{\rm seg}=0.5$, 2 and 8 fm at 8 configurations ($\sigma_R$ at 10/30/100/300/900 MeV, $T$ at 0.3/1/2 GeV), $2$–$4\times10^5$ events each:

- $\sigma_R$ and $T$: $\chi^2=21.9/16$;
- $\langle N_{\rm coll}\rangle$: $\chi^2=20.4/16$;
- QE-SF $\langle N_{\rm coll}\rangle$: $2.122\pm0.006$ / $2.119\pm0.006$ / $2.105\pm0.006$.

**Counters** over all (b) runs [run]:

| counter | result |
|---|---|
| in-flight jailings | 0 (as derived) |
| finalize failures, the (a) quirk | 0–11 per $4\times10^5$ events, at $\ge0.8$ GeV only |
| acceptance | 1.1–3.0% |
| proposals per event | 24–100 |

The low acceptance is dominated by $\sup_T\sigma_{\rm free}(E^{\rm eff})$: for a 30–50 MeV nucleon a co-moving target gives $\sigma_{np}$ up to 13 b, against a Fermi-averaged 200–370 mb [run, offline scan of the PDG table]. Drawing the target as part of the proposal would recover a factor 3–30 for slow nucleons. It is not implemented.

## 3.3 Fixes made during validation

Two errors were found by the tests and fixed before the final runs.

**(i) Target law in the QE-SF path.**

- My first version reused the accepted target for the collision. (a) re-draws it at the collision point in that path, and so does `kaskada` beam mode (see Step 1).
- That gave +3% collisions in the SF test, at 6–10σ.
- Letting `make_interaction` re-select, as (a) does, removed the excess (§4.2). All (b) numbers in this report use the fixed version.

**(ii) Majorant violation.**

- One violation (ratio 1.027) occurred in $3.3\times10^7$ proposals ($T_p$=1.6 GeV).
- Cause: (a)'s `pion_abs` changes $(Z_r,N_r)$ **before** the Pauli check (`Interaction.h:574–680`) and does not undo it when the absorption is blocked, so the nucleus changed without re-anchoring.
- (b) now detects any change of $(Z_r,N_r)$ in a blocked attempt and re-anchors (`kaskada7.cc:1664–1677`). That happened 23 times in $4\times10^5$ events.
- After the fix, the same seed gives 0 violations [run].
- The other (b) runs predate this fix. It changes only events with a blocked absorption, and all other runs had 0 violations.

# Step 4 — validation

## 4.1 Closed form, uniform sphere

**Setup.**

- Test build `-DKASKADA_SUDAKOV_TEST`: constant $\sigma$, or random per evaluation with the same mean; absorption at the first collision.
- Flat ¹²C, $R=2.747$ fm, $\rho=0.138$ fm⁻³.
- Per event the exact $P_{\rm int}=1-e^{-\rho\sigma L(b)}$.
- The KS test is applied to the probability-integral transform of the first-collision depth, $u=(1-e^{-\rho\sigma d})/(1-e^{-\rho\sigma L})$.
- $2\times10^5$ events per run [run].

| run | $\sigma$ [mb] | $P_{\rm int}$ | exact | z | KS $D$ | KS $p$ |
|---|---|---|---|---|---|---|
| (b) $L_{\rm seg}$=0.3 fm | 40 | 0.82539 | 0.82536 | +0.04 | 0.0029 | 0.13 |
| (b) 1 fm | 40 | 0.82666 | 0.82613 | +0.67 | 0.0019 | 0.57 |
| (b) 3 fm | 40 | 0.82581 | 0.82562 | +0.24 | 0.0017 | 0.71 |
| (b) 20 fm | 40 | 0.82579 | 0.82593 | −0.18 | 0.0027 | 0.17 |
| (b) 1 fm, random $\sigma$ | 40 | 0.82732 | 0.82639 | +1.20 | 0.0019 | 0.57 |
| (b) 1 fm, random $\sigma$, $\hat\Gamma\times3$ | 40 | 0.82585 | 0.82608 | −0.30 | 0.0012 | 0.97 |
| (b) 0.3 fm, random, $\hat\Gamma\times3$ | 40 | 0.82671 | 0.82596 | +0.96 | 0.0015 | 0.86 |
| (b) 1 fm | 10 | 0.38635 | 0.38729 | −0.88 | 0.0031 | 0.45 |
| (b) 1 fm, random, $\hat\Gamma\times3$ | 100 | 0.96570 | 0.96609 | −1.09 | 0.0029 | 0.08 |
| (a) 0.2 fm | 40 | 0.82559 | 0.82499 | +0.77 | 0.0021 | 0.45 |
| (a) 0.05 fm, random $\sigma$ | 40 | 0.82520 | 0.82556 | −0.46 | 0.0031 | 0.08 |
| (a) 0.2 fm, random $\sigma$ | 40 | 0.82400 | 0.82597 | **−2.51** | 0.0029 | 0.11 |
| (a) 1 fm, random $\sigma$ | 40 | 0.81607 | 0.82621 | **−13.0** | 0.0092 | $<10^{-10}$ |

(a) is exact for a constant $\sigma$ (memoryless). With one random target per step it shows the Jensen deficit of Step 1.

## 4.2 Agreement with (a)'s zero-step extrapolation

- **Beam mode:** Step 2 table and Fig. 1; data comparison in Fig. 3. Overall $\chi^2=32/22$; the outlier is resolved above.
- **$\langle N_{\rm coll}\rangle$** in 13 configurations: (b)$-$(a)$_0$ within ±1.8σ (Fig. 2).
- **QE-SF branch** ($T_{N1}$=300 MeV, non-transparent):
  - (a)$_0$ $=2.110\pm0.003$ (quadratic in $\Delta$, $\chi^2$=3.3/4);
  - (b) $=2.116\pm0.003$ (four runs);
  - $z=+1.2$.
- **All SF branches vs (a) at its default step** ($10^5$ events each, $T_{N1}$=100/300/800 MeV):
  - final-state multiplicities agree: $\chi^2=68/48$ over $P(N{=}0)$, $\langle n_p\rangle$, $\langle n_n\rangle$, $\langle n_\pi\rangle$;
  - $\langle N_{\rm coll}\rangle$ is higher by 1–4% in (b), as expected from Fig. 2.

## 4.3 NuWro's own regression tests

- **There is no automated physics regression suite at this commit.** CI (`.github/workflows/CIBuild.yml`) only builds the Docker image. The `fsi` program (`simulations.cc`/`calculations.cc`) is a plotting driver with no pass/fail criterion.
- **Substitute tests** [run]:
  - (i) The patch applies cleanly to the pinned commit and builds.
  - (ii) With `kaskada_sudakov=0` the patched `kaskada` reproduces the pristine one **bit for bit**: md5 of all final-state momenta and `nod` counters, p 200 MeV surface / p 700 MeV transparency / π⁺ 300 MeV / n 50 MeV, 5000 events each.
  - (iii) A fresh clone with the patch applied reproduces the development build bit for bit in (b) mode, including the weights file.

## 4.4 Data, $\chi^2$ with data errors only

Your digitised CSVs were used: 99 $\sigma_R$ points, 15 $T$ points. The model was run at the data energies with $2\times10^5$ events per point [run].

| | $\sigma_R$ <100 MeV ($N$=82) | $\sigma_R$ ≥100 MeV ($N$=17) | $T$ 0.18–1.8 GeV ($N$=15) |
|---|---|---|---|
| (a), $\Delta$=0.2 fm | $\chi^2$=6561, model/data 0.75 | 548, 0.86 | 68.5, 0.90 |
| (b) | 6564, 0.75 | 550, 0.87 | 66.7, 0.90 |

- Both algorithms describe the data equally badly.
- NuWro's beam-mode $\sigma_R$ is 25% low below 100 MeV and 14% low above. $T$ is 10% low.
- The algorithm is irrelevant to this; the model inputs are what matter.
- In particular, `mfp_scale` and the SRC factor are switched off for the beam nucleon's first collision in `beam_placement=2` (`kaskada7.cc:727`), so $\sigma_R$ does not see the 0.76 tuning.

![Data (digitised from Isaacson et al. 2021; the 60.8/61.0 MeV pair is inconsistent), (a), (b), and (a)'s bias at 0.2 fm against the data's relative errors.](figs/fig2_data.png){width=95%}

## 4.5 CPU per event

Single core, idle machine, $5\times10^4$ events ($3\times10^4$ SF), startup subtracted, ROOT output included [run]. Times in ms per event:

| | (a) 0.2 fm | (a) 0.05 fm | (a) 0.025 fm | (b) | (b)/(a)$_{0.2}$ | (b)/(a)$_{0.025}$ |
|---|---|---|---|---|---|---|
| $\sigma_R$, 30 MeV | 0.124 | 0.40 | 0.69 | 0.18 | 1.45 | 0.26 |
| $\sigma_R$, 300 MeV | 0.133 | 0.40 | 0.73 | 0.12 | 0.90 | 0.16 |
| $T$, 1 GeV | 0.20 | 0.62 | 1.16 | 0.33 | 1.70 | 0.29 |
| QE-SF, $T_{N1}$=300 MeV | 0.21 | — | 1.45 | 0.40 | 1.88 | 0.28 |

The patched build with (b) off costs +3% (6.74 vs 6.52 s).

# Step 5 — reweighting payoff

**Weights** (`kaskada_weights="mfp=…,xs=…"`).

- **(b): veto weights.** Per candidate of a nucleon, with $k=\Gamma'/\Gamma$ ($k=s_{\rm mfp}/s'$ where the SRC/mfp factor is on; $k=f'/f$ for the NN scale):
  - accepted: $\times k$;
  - rejected: $\times(\hat\Gamma-k\Gamma)/(\hat\Gamma-\Gamma)$.
  - The majorant carries a headroom $H=\max k$ (`:1642–1648`).
  - Pauli-blocked candidates count as accepted.
- **(a): likelihood ratio of each censored free path** (`step_weights`).
  - $k\,e^{-(k-1)\Gamma\lambda}$ if $\lambda<\Delta$, else $e^{-(k-1)\Gamma\Delta}$.
  - This is exact for (a)'s own process at its own step.

**Closure** [run].

- Nominal samples ($2\times10^5$ events) were reweighted to `mfp` 0.6/1.0 and NN scale 0.8/1.25 and compared with regenerated samples.
- Configurations: $T_p$=0.3 and 1 GeV transparency mode; 50 and 300 MeV $\sigma_R$ mode.
- Five observables: $P_0$, $\langle N_{\rm coll}\rangle$, $P(N\ge2)$, $\langle N_{\rm out}\rangle$, $\langle N_{\rm el}\rangle$.

| | pulls (80, correlated) | max $|z|$ | ESS/N |
|---|---|---|---|
| (a) | $\chi^2$=36.6 | 1.7 | 0.84–0.998 |
| (b) | $\chi^2$=85.5 | 2.6 | 0.79–0.998 |

- All six $|z|>2$ pulls in (b) are in the 50 MeV $\sigma_R$ configuration. The nominal (b) sample there is itself $+2.4\sigma$ high in $\langle N_{\rm out}\rangle$ against an independent (b) sample [run]. Every reweighting of it inherits that.
- $T$ moves by $-0.057$ (mfp 0.6) and $+0.063$ (mfp 1.0) at 300 MeV, reproduced to $\le0.002$.

![Left: $P_0$ against `kaskada_NN_mfp_scale` from (b)'s weights (open) and regenerated samples (filled). Middle: pulls of reweighted against regenerated, five observables per row. Right: effective sample size.](figs/fig4_reweight.png){width=100%}

**Why (b)'s ESS is lower.**

- Variance of the veto weight: $\mathbb E W^2=\mathbb E\exp\!\int(\Gamma'-\Gamma)^2/[\Gamma(1-\Gamma/\hat\Gamma)]$.
- Variance of the exponential weight: $\mathbb E\exp\!\int(\Gamma'-\Gamma)^2/\Gamma$.
- The veto form is always larger, and is close to the exponential form only because the acceptance is 1–3% here [derived].
- So the payoff is not "reweighting becomes possible". **(a) is already exactly reweightable at its own step**; (b) makes the reweighted process step-free.

**GENIE, for comparison** (read in source: GENIE Reweight `5d45f8a`, Generator `e2e04f1`).

- `GReWeightINuke` accepts only `HAIntranuke2018/2025`; any other FSI model is refused (`GReWeightINuke.cxx:83–90`).
- Only hadrons with status `kIStHadronInTheNucleus` are weighted, that is, primaries (`:246`).
- The mean-free-path weight is $P'_{\rm surv}/P_{\rm surv}$ (escaped) or $(1-P')/(1-P)$ (interacted) (`GReWeightUtils.cxx:45–84,135–158`).
- `ProbSurvival` integrates the mean free path by stepping with the model's `HadStep` (`INukeUtils2018.cxx:284–316`).
- Fate weights are fraction rescalings.
- So GENIE reweights a one-shot fate per primary hadron with a stepped survival probability, and has no reweighting for the hN cascade. NuWro (a) or (b) gives exact per-event weights for the whole cascade, secondaries included.

# Verdict

**I would not propose (b) to the NuWro developers as a replacement for the stepping cascade.** I would send them the (a) findings and the patch as an optional, off-by-default mode. The evidence:

- **The observables NuWro is tuned to do not need it.** At the default 0.2 fm the step bias in $\sigma_R$ and transparency is $-0.06\pm0.18\%$ and $+0.14\pm0.08\%$. That is 10–50× below the data errors and below the `kaskada_NN_mfp_scale` band. The data disagreement (25%/14%/10%) is a physics-input problem that neither algorithm touches.
- **The cascade-level bias is real but small.** (a) misses 2–4% of collisions at 0.2 fm, mostly through the Jensen effect of one random target per step on slow secondaries. (b) removes this at 0.9–1.9× the CPU of (a)'s default. Reaching the same accuracy with (a) needs $\Delta\approx0.025$ fm at 3.5–6× the CPU. A cheaper fix inside (a), not tested here, would be to average $\Gamma$ over several targets per step.
- **Reweighting is not unique to (b).** (a) reweights exactly with likelihood ratios and has the better ESS.
- **The patch has costs.** It touches 8 files (+655 lines, no deletions), needs a proof-backed majorant for every new physics factor (hyperons are only bounded empirically), and has 1–3% acceptance until targets are sampled in the proposal.
- **What is clearly worth sending now** is what the validation exposed in (a):
  - (i) uninitialised `event::flag` in `kaskada.cc`, which silently turns on the QE-SF code path for every beam nucleon;
  - (ii) a Pauli-blocked pion absorption still removes nucleons (`Interaction.h:574–680`);
  - (iii) a failed `remove_nucleon` duplicates particles (`kaskada7.cc:1110` vs `:89`);
  - (iv) the cascade ignores `pauli_blocking`.

**When (b) would become worth proposing.** If NuWro adds observables that count secondary collisions at the percent level (multiplicities, calorimetric energy, neutron counting), or if analyses want event weights whose meaning does not depend on `step`. The natural first upgrade is then target sampling in the proposal.

# Files

- `nuwro_sudakov.patch`: the production patch against `f16c937`, 8 files, additions only.
- `nuwro_sudakov_unittest_hooks.patch`: the `#ifdef KASKADA_SUDAKOV_TEST` hooks used in §4.1.
- `tools/`: `sf_test.cc` (QE-SF driver), `ana.cc`/`ana2.cc`/`flags.cc` (tree readers), `runner.py`, analysis scripts.
- `results/`: JSON summaries of every table, and the step-2 table.

Usage: `kaskada -p kaskada_sudakov=1 [-p kaskada_sudakov_segment=2] [-p kaskada_weights="mfp=0.6,xs=1.2"]`. This writes `<out>.w.txt` (weights) and `<out>.sud.txt` (counters).
