# Log

Chronological journal. Failed attempts are recorded with the same prominence as successful ones,
per `CONTRIBUTING.md` §5. Each entry states the date, the issue worked on, what was done, what
came out, and the commit.

---

## 2026-08-11 — repository initialised

Scaffolding created: `PROBLEM.md` (formal problem statement, method, acceptance criteria),
`CONTRIBUTING.md` (working protocol), `CONVENTIONS.md` (stub, to be filled by issue #3),
`tex/main.tex` (skeleton with epistemic status tags), `src/`, `lit/notes/`, `results/`.

Issue tree opened under epic #1. No mathematical work done yet.

## 2026-08-11 — step 0: issue hierarchy and `gap` label

Attached issues #4–#12 to epic #1 as sub-issues (GitHub sub-issue hierarchy; #2 and #3 were
already attached). The GitHub MCP tooling does not expose issue database ids directly; they were
recovered by decoding the msgpack pagination cursors of the issue list and verified against the
ids returned for the existing sub-issues #2/#3. Created the `gap` label (no label-creation API
available: set it transiently on #12 via the issues API, which auto-creates missing labels,
verified with the label-read endpoint, then removed it from #12).

## 2026-08-11 — issue #2: literature status

Established that the problem is **still open** as of 2026-08-11. arXiv:2512.10101v2 (27 Apr 2026)
argues, does not prove, the strange-series AdS/dS interpolation (their own outlook: "in some sense
interpolates"; "It is currently unknown what the correct physics interpretation of such an
amplitude would be"). All ten INSPIRE-listed citing papers full-text grepped: none engages the
strange series (two "strange" hits were false positives). No newer Schouten/Isachenkov preprint;
the announced de Groot–Isachenkov–Posthuma noncommutative-geometry work has not appeared (only the
classical precursor arXiv:2601.22171 by de Groot). Isachenkov's SCGP talk (15 May 2026,
"Non-commutative geometry of the DSSYK model") is on video (SCGP event 496, video 7653) with no
associated preprint found. arXiv/Semantic Scholar searches for strange-series work in math.QA/OA:
nothing new. Deliverable: `lit/README.md`. Failed/limited checks recorded there (no full-text
check possible for one record without arXiv id; math-journal coverage weaker).

## 2026-08-11 — issue #3: conventions

Filled `CONVENTIONS.md` and the mirrored Conventions section of `tex/main.tex` (Conventions
2.1–2.8), following arXiv:2512.10101v2 verbatim, with the GKK (arXiv:0905.2830v2), BINN
(arXiv:2212.13668) and MMNNSU dictionaries worked out explicitly. Every dictionary row that could
be checked numerically was checked: `src/check_conventions.py`, 63/63 checks pass
(`results/check_conventions.out`); largest deviation among passing checks 5.2e-14, most below
1e-45 (60-digit arithmetic where exact cancellations of size q^{-2n} occur).

Findings worth recording:
- Generator map SI↔GKK is E_SI = E_GKK, F_SI = -F_GKK, K equal; symmetrised Casimir
  Ω̃_SI = -Ω_GKK. BINN's Casimir equals SI's Ω exactly (A=K, B=E, C=F, q̂=q); MMNNSU's C equals
  SI's Ω.
- In SI conventions the series are separated by the Ω̃ eigenvalue: principal [-1,1]; discrete and
  complementary > 1 (distinguished by one- vs two-sided K-spectrum); **strange < -1, and Ω̃ < -1
  happens only for the strange series** — "this representation is strange" is decidable from
  (Ω̃, spec K). Strange labels: a > 0 continuous in the abstract classification; a = ℓ - 1/2 with
  SI's ℓ; **in the Plancherel decomposition of the regular representation the strange series
  occurs only at discrete labels a ∈ ½Z_{>0}** (GKK Thms 4.6–4.9 analogues).
- GKK's stated principal-series Casimir eigenvalue μ(q^{2ib}) disagrees by a sign with their own
  printed matrix elements; settled numerically (check C1) in favour of the matrix elements:
  Ω_GKK(π_{b,ε}) = -cos(2b ln q). No effect at the level of the series family. The
  complementary-series eigenvalue (unstated in GKK) computed the same way: -μ(q^{1+2λ}).
- First run of check (C1) in float64 FAILED with deviations up to 1e46 at q=0.3 — catastrophic
  cancellation between q^{-2n}-sized terms, not a mathematics error; re-run in 60-digit mpmath
  passes. Recorded as required by CONTRIBUTING.md §5 (failed checks reported).
- Also fixed `.gitignore`: a bare `*.out` was silently excluding `results/*.out` from version
  control; LaTeX artefact patterns are now scoped to `tex/`.

Carried forward (flagged for issues #4 and #8, per coordinator note):
1. "dS region" is not agreed in the literature — upper spectral edge θ = |ln q²|k
   (arXiv:2411.16922, arXiv:2505.08116) vs middle of the band θ = π/2 (arXiv:2310.16994 line).
   If the claim's truth value depends on the choice, that is itself the finding for #8.
2. The strange series enters the Plancherel formula only at discrete labels a ∈ ½Z_{>0}, while
   "interpolates" is prima facie a continuous notion; #4 must confront whether a discretely
   indexed family can interpolate at all.

## 2026-08-11 — issue #4: formalising the claim

Wrote `lit/notes/2512.10101.md` (verbatim quotation ledger with locations, proved-vs-argued
split, conventions map, and a list of suspected typos in v2) and Section 1 of `tex/main.tex`.

Structure of the formalisation:
- Definitions 1.1–1.3 fix the edge regimes, the AdS region (lower edge; uncontroversial), and
  the two rival dS identifications (edge-dS of arXiv:2411.16922/2505.08116, used by the paper;
  centre-dS of arXiv:2310.16994).
- Lemma (q-Mehler form): both bilocal kernels are q-Mehler kernels; the strange one at
  t = -q^{2l}, i.e. strange matter = alternating chord weights (-1)^k q^{2kl}. DERIVED.
- Lemma (reflection identity): K^S(th1,th2) = K^O(th1, pi-th2) EXACTLY at finite q — the strange
  kernel is the discrete kernel with one spectral argument reflected through E -> -E. This
  sharpens the paper's two q->1 statements (4.19)/(4.20) into one finite-q identity. DERIVED.
- Proposition (limits): the paper's (4.19)/(4.20) re-derived as normalisation-free ratio
  statements; (4.20) upgraded to an exact identity K^S(L(k1),U(k2)) = K^O(L(k1),L(k2)).
  Also the converse (iv): discrete-series matter does NOT connect opposite edges. DERIVED;
  numerics in `src/check_claim_kernels.py`.
- Claim C-amp (the paper's literal claim): settled by the Proposition, conditional ONLY on the
  edge-dS identification. Claim C-geo (geometric interpolation via the radial Casimir problem,
  what "genuinely interpolate" ought to mean): CONJECTURED, open — the project's target.
- Remarks: dS-identification sensitivity (under centre-dS the amplitude argument does not go
  through as stated — quantitative version deferred to #8/#11); the discreteness tension (which
  continuous parameter is supposed to cross: abstract label a — discrete in the Plancherel
  decomposition; spectral parameter theta; radial coordinate); two degenerate readings excluded
  as contentless; open sub-question flagged for #6 (principal-series bilocals' mixed-edge
  kernels — not computed anywhere we know of, needed for the selectivity statement).

Suspected typos found in v2 while re-deriving (recorded in the lit note §7): the Gamma_{q^2}
shift in (4.19) should be i pi/(2 ln q), not i pi/ln q (the full period, which would not vanish);
mechanism unaffected.

Numerics record for #4 (failures included, per CONTRIBUTING.md §5):
- First run of `check_claim_kernels.py`: 6 of 36 K1 checks FAILED at q=0.9 (deviations up to
  2.5e-18). Diagnosis: the (a;q)_inf implementation capped its product at a dps-independent
  number of terms; the truncated tail contributed relative error ~1e-25, amplified by kernel
  values ~1e7. Fixed by scaling the term count with working precision; all K1 checks then pass
  (deviations < 6e-33).
- The reflection-identity check near q=1 was first written with an absolute tolerance while the
  unrenormalised kernels grow beyond 1e+10000; rewritten as a relative deviation (now < 4e-38 at
  all sampled q).
- The first version of the dS-centre remark (Remark on dS-identification dependence) argued the
  wrong direction: a sign error in the ln(u;q^2)_inf ~ -Li_2(u)/|ln q^2| asymptotics suggested
  the AdS<->centre strange amplitude diverges; direct computation shows it VANISHES like
  exp(-3pi^2/(4|ln q^2|)) (measured rates -4.97/-7.07/-7.22 at q=0.9/0.99/0.995 vs predicted
  -7.40). The remark in tex/main.tex states the corrected version with the measured numbers.
  Conclusion unchanged: under the centre-dS identification the amplitude argument fails.

## 2026-08-11 — issue #5: the quantum homogeneous space

Section 3 of `tex/main.tex`; CONVENTIONS.md §7 updated to the fixed choice.

The DSSYK space is the twisted-primitive ("(s,t)-Gauss") double quotient of SU_q(1,1)⋊Z_2 with
invariant algebra ⟨ρ_st, e⟩ (SI (3.47)-(3.48)), in the rescaled regime: two-sided version = the
chord sector at ε=-1 (SI (4.1)-(4.5)); one-sided version = SI's reduced quantum AdS_{2,q}, the
restricted q-lattice R²_{q²}(ξ) (4.22) with the q-exponential-weighted inner product (4.24),
Casimir action (4.25) (the radial-type q-difference operator), classical coordinates (4.30).
Excluded with sources: the quantum disk (Cartan quotient; SI §5 verbatim: "distinct from the
DSSYK coset considered in this paper") and the quantum hyperboloid (occurs in no source's DSSYK
reduction). Recorded SI's own caveat that the position-space picture rests on an *assumed*
q-Fourier transform (their §4.4 "is currently not available to us") — we stay on the momentum
lattice; a gap issue is required if any later step needs position space.

Sharpest content: where the strange series lives. PROVED (elementary, via the determinate
q²-Hermite moment problem / bounded Jacobi operator): the chord-sector state decomposition is
purely principal (a.c. spectrum [-1,1], multiplicity one) — the strange series does NOT occur as
states of the DSSYK sector. It occurs (GKK, quoted) as discrete summands of L²_q(SU(1,1)⋊Z_2)
in the sectors (+,-),(-,+),(+,+), and acts on the chord sector through the bilocals. Consequence
recorded for #6-#8: the radial problem is the same q-difference operator for every series;
series membership enters only through the spectral value Ω̃ (strange ⇔ Ω̃ < -1, below the
principal band). Any interpolation mechanism through strange *states* of the transfer matrix is
excluded.

## 2026-08-11 — issue #15: the normaliser origin of the reflection identity

New subsection of `tex/main.tex` §1 (`sec:normaliserorigin`), connecting Lemma 1.5 to the
Z₂-grading of GKK's dual von Neumann algebra. Sources re-read from the arXiv LaTeX (GKK
0905.2830v2 and SI 2512.10101v2), not from paraphrases.

Results (tags as in the tex):
- PROVED (lem:gradingchord): the extra factor in SI's strange bilocal (4.18) is
  Π = (2qρ̂_st)^{-iπ/(2 ln q)} = chord parity (-1)^k — the nontrivial character of the
  spectral lattice q^{-2Z} — a self-adjoint unitary with ΠΩ̃Π = -Ω̃ and Π|P^θ⟩ = |P^{π-θ}⟩.
- PROVED (prop:factorisation): S^ℓ(0) = O^ℓ(0)Π. The answer to "is Lemma 1.5 the statement
  that S^ℓ is O^ℓ conjugated by a normaliser element?" is NO: conjugation by Π fixes O^ℓ(0)
  (diagonal operators commute); the correct statement is the FACTORISATION through the
  grading unitary. The reflection identity in chord-sum form is now proved without q-Mehler.
- PROVED (in GKK; transcription + label computation, thm:fusion): for every ℓ ∈ ½Z_{≥1} the
  strange and discrete series are FUSED into a single irreducible corepresentation of the
  normaliser quantum group: W_{p,x}|_{U_q} ≅ π^S_{ℓ-½,ε(p)} ⊕ D⁻_ℓ ⊕ D⁺_ℓ — with exactly
  SI's label map a = ℓ - ½. Every strange occurrence in the Plancherel decomposition is of
  this form (GKK's decomposition theorem + principal corepresentations contain no strange
  constituents). GKK's own verbatim summary: "the discrete series are no longer split up
  into a positive discrete series and a negative discrete series."
- PROVED (thm:covobstruction): any AdS/dS criterion covariant under the normaliser grading
  assigns the same type to π^S_{ℓ-½} and D^±_ℓ. Radial-problem version
  (cor:radialparity): the parity map intertwines the chord-lattice radial problems at ω and
  -ω preserving all normalisability data. ARGUED (rem:obstructionreading): this is a
  candidate obstruction for Claim C-geo; the only grading-breaking datum identified so far
  is the orientation of the spectral axis, i.e. the external dS proposal itself. Issue #8
  reformulated accordingly (comment on #8; body edit pending permissions, see below).
- GAP (issue #16): whether Π is the literal r→∞ rescaled limit of a concrete element of
  M̂₋ (SI's rescaling is done on states and ρ̂_st, not on the dual algebra). The structural
  identification (anticommutation + essential uniqueness, lem:uniqueness via simplicity of
  the chord spectrum) does not depend on it.

Numerics: `src/check_normaliser_grading.py`, 72/72 checks pass
(`results/check_normaliser_grading.out`); worst deviations ~1.5e-59 on 60-digit arithmetic.
Failed first attempts, recorded per CONTRIBUTING.md §5:
- (G7) first run FAILED at (q,χ,j)=(0.6,0,2) with deviation 1.37e-48 against tol 1e-50: the
  coefficient arrays grow like q^{-2m} (size ~1e13), so an absolute tolerance was wrong;
  comparisons made relative, after which the same points pass at ≤4e-61.
- (G4) the "conjugation reading is false" check was first written as a MINIMUM of absolute
  differences over odd k, which decays like q^{2ℓk} and FAILED (2.7e-11, 6.6e-5 at the two
  sample points) — a check-design error, not a mathematics error: the correct invariant is
  the relative difference |ΠOΠ⁻¹ - S|/|S| = 2 on odd k, which passes exactly (deviation 0).

Environment facts recorded:
- LaTeX cannot be compiled locally: pdflatex fails at `\usepackage{mathtools}` (TeX Live too
  old; `tlmgr install` cross-release error). Document integrity verified instead by a
  label/ref/cite/environment consistency check (no duplicate labels, no unresolved refs, no
  cites missing from refs.bib, all environments balanced). A committed checker script is
  planned with the CI work.
- The working GitHub account (`olegnpavlov`) turned out to have read-only access:
  `git push` → 403, issue-body edits → permission denied, `gap` label on #16 silently
  dropped at creation. Issue #17 opened for the owner. Commits are made locally on the
  tracked branch and will be pushed when access is granted; scope changes are recorded as
  explicit issue comments (done for #8); gap issues carry "Gap:" in the title.

## 2026-08-11 — issue #18: the incomplete proof of Proposition 1.6(iv) finished

The proof of parts (ii)/(iv) of `prop:limits` contained a literally unfinished formula
(`|1+q^{2x+2n}| ≥ 1+cos(2k ln q⋯)q^{2ℓ+2n}` with an ellipsis) followed by a second,
gestural argument ("A cleaner route..."). Replaced by ONE complete elementary proof: since
the phase φ_q = 2(±k₁±k₂)ln q of q^{2x+2n} is independent of n, each per-factor modulus
ratio is (1-2ρcosφ+ρ²)/(1+2ρcosφ+ρ²); under the explicit hypothesis 2(k₁+k₂)|ln q| ≤ π/3
one has cosφ_q ≥ ½, giving per-factor bound g(ρ) = (1-ρ+ρ²)/(1+ρ+ρ²) (g strictly
decreasing, derivative numerator 2ρ²-2 < 0), hence ≤ g(½) = 3/7 on the ⌊ln2/(2|ln q|)-ℓ⌋+1
lattice sites with ρ_n ≥ ½ and ≤ 1 elsewhere; so |r_x(q)|² ≤ (3/7)^{⌊N_ℓ(q)⌋+1} → 0. The
q-gamma "cleaner route" gesture dropped (the factual note on SI's (4.19) misprint kept).
Tag bookkeeping: (ii)/(iv) vanishing now proved in full; prop:limits keeps DERIVED because
(i)/(iii) rest on Γ_{q²}→Γ and the q-Mehler lemma.

Numerics: new check (K5) in `src/check_claim_kernels.py` verifies the smallness hypothesis,
every per-factor inequality, and the final bound at q = 0.9/0.99/0.999: worst |r_x|² =
2.2e-7 vs bound 7.9e-2, 8.3e-103 vs 3.1e-13, 3.2e-1066 vs 4.8e-128 (bound valid, far from
tight, as expected of a counting bound). Output regenerated exactly as emitted; the K1-K4
portion is bit-identical to the previously committed file except that a trailing `EXIT=0`
line, which the script never printed (a hand edit), is gone. The systematic
committed-outputs policy is handled under the repository-integrity issue.

## 2026-08-11 — issue #19: PROBLEM.md attribution fix

PROBLEM.md §5 attributed arXiv:2511.03779 (JHEP 05(2026)080) to Jiuci Xu, contradicting
lit/README.md (issue #2), which had already established the author as Sergio E.
Aguilar-Gutierrez. Re-verified against the arXiv API today: sole author Aguilar-Gutierrez,
title "Cosmological Entanglement Entropy from the von Neumann Algebra of Double-Scaled SYK
& Its Connection with Krylov Complexity". Corrected, with the actual title spelled out.
Also verified the van der Heijden et al. reference at PROBLEM.md §5 against arXiv:2511.08743:
title "Quantum Symmetry and Geometry in Double-Scaled SYK" (van der Heijden, E. Verlinde,
J. Xu) — was already correct; arXiv id added so the reference is checkable. No other
PROBLEM.md changes.

## 2026-08-11 — issue #20: bit-for-bit reproducible committed outputs

`results/check_conventions.out` did not reproduce across machines: the (C2)-(C5) sections
used numpy float64, whose ~1e-17 residuals depend on BLAS summation order (committed
1.787e-18 for "K F = q^{-1} F K" reran as 3.490e-21 elsewhere). Ported run_C234 and run_C5
to pure mpmath (60 digits, software arithmetic, platform-independent), with the checks
semantically unchanged (same operators, truncation N=40, interior window, tolerances).
Numpy dependency dropped from the script. Deviations now ~1e-61. Verified: the (C1) section
of the regenerated output is bit-identical to the committed one, and an immediate rerun of
the regenerated script is bit-identical to its own output. Policy now holding for all four
check scripts (all pure mpmath): a committed output is exactly what the committed script
emits, reproducible bit-for-bit.

## 2026-08-11 — issue #21: uncited refs.bib keys

Five keys were defined in tex/refs.bib but cited nowhere in tex/main.tex, including
Schlosser:2024matsuki — the method paper PROBLEM.md §4 rests on. All five belong in the
document, so they are now cited where they belong rather than removed: the Section 5
(radial part) placeholder states the fixed method with citations to Schlosser:2024matsuki,
Isachenkov:2016superint, Isachenkov:2017integrability; the Section 7 (verification)
placeholder states the fixed cross-check targets with citations to
Blommaert:2024sinedilaton and Xu:2024vnalgebras. Neither sentence pre-empts issues #7 or
#9–#11. Zero uncited keys and zero unresolved cites, verified.

## 2026-08-11 — issue #6: matrix coefficients of the strange series

Section 4 of `tex/main.tex`. Sub-issues opened rather than solved inline: #13 (principal-series
bilocal kernels / selectivity, flagged in Remark 1.12) and #14 (`gap`: no Plancherel/orthogonality
theory for the strange family after the DSSYK rescaling — the project-level instantiation of the
Plancherel caveat PROBLEM.md quotes).

Results:
- The strange series transcribed into our conventions (PROVED; transcription verified by #3's
  checks C1-C2).
- Radial coefficients on the chord lattice at the strange value Ω̃ = -mu(q^{2a}): the regular
  solution is the continued q²-Hermite (-1)^k H_k(mu(q^{2a});q²), grows like q^{-2ak}
  (no strange states — consistent with #5); NEW closed form for the minimal (decaying) solution
  c^S_k(a) = (-1)^k q^{2ak} 2phi1(0,0; q^{4a+2}; q², q^{2k+2}), PROVED by an explicit
  coefficient recursion with a convergence argument; Wronskian W_k = (q²;q²)_k W_0 (PROVED) with
  W_0 = (s^{-1}-s)/((s²;q²)_inf(q²;q²)_inf) and the connection formula H_k(mu(s);q²) =
  c^min_k(s^{-1})/(s²;q²)_inf + c^min_k(s)/(s^{-2};q²)_inf (both DERIVED: asymptotic matching,
  verified to 1e-33..1e-40 at multiple (q,s) — `src/check_strange_coefficients.py`,
  `results/check_strange_coefficients.out`, all checks pass).
- Resonance observation: the connection formula degenerates exactly at s² ∈ q^{2Z} — which
  includes every strange label in the Plancherel decomposition (a ∈ ½Z_{>0}); the minimal
  solution itself stays regular there. Recorded in the theorem.
- The rho-shift (operator weight q^{2lk} vs wavefunction decay q^{(2l-1)k}) recorded as a
  remark with a numerical check (constant = 1.0), so it is not mistaken downstream for an error.
- Failed first attempt recorded: the S2 growth-rate check used a fixed absolute tolerance at
  k=60 and FAILED at (q,s)=(0.8,0.9) and (0.9,0.5); the deviations matched the predicted
  correction scale max(s²,q²)^k exactly (3e-6 at 0.81^60), i.e. a premature-asymptotics
  tolerance, not a formula error; fixed by comparing at k=200 with the predicted scale.
- No analytic continuation from compact quotients is used anywhere (the 2412.19681 caveat does
  not apply); stated in the tex.
