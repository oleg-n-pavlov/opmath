# Conventions

**Status: fixed by issue #3 (2026-08-11). Verified numerically by `src/check_conventions.py`
(output: `results/check_conventions.out`, 63/63 checks pass). Mirrored in `tex/main.tex` §2;
the two must agree exactly.**

Unless stated otherwise, we follow **arXiv:2512.10101v2** (Schouten–Isachenkov, "SI") verbatim.
Every deviation from any cited source is stated at the point where it occurs, and collected in the
translation table (§10). Equation references "SI (x.y)" are to the compiled v2 numbering
(LaTeX labels are also given, since v2's source is public); "GKK" = arXiv:0905.2830v2
(Groenevelt–Koelink–Kustermans); "BINN" = arXiv:2212.13668 (Berkooz–Isachenkov–Narayan–
Narovlansky); "BK" = Burban–Klimyk, J. Phys. A 26 (1993) 2139–2151; "MMNNSU" = Masuda–Mimachi–
Nakagami–Noumi–Saburi–Ueno, Lett. Math. Phys. 19 (1990) 187–194 and 195–204.

## 1. The deformation parameter and the DSSYK dictionary

- `q` is **real** with `0 < q < 1` [SI §3.1, first paragraph].
- The DSSYK chord parameter is `q_SYK := q^2` [SI §4.2: "This is, of course, exactly the
  partition function for the double-scaled SYK model (under the redefinition q^2 → q)"].
- `q = e^{-λ/2}` where λ is the DSSYK coupling, `λ = 2p^2/N`; equivalently `q_SYK = e^{-λ}`
  [SI §4.1: "where q = e^{-λ/2} and λ is the renormalization parameter"; PROBLEM.md §1].
- DSSYK energy: `E(θ) = 2 cos θ / √(1-q^2) = 2 cos θ / √(1-q_SYK)`, `θ ∈ [0, π]` [SI §4.2,
  below the partition function]. The **lower edge** of the spectrum is `θ = π` (E = -E_max),
  the **upper edge** is `θ = 0` (E = +E_max). The Schwarzian/AdS₂ regime is the lower edge:
  `θ = π - |ln q^2| k`, `k ≥ 0` fixed, `q → 1` [SI §4.2, "Schwarzian limit"]. The proposed dS
  regime of Blommaert–Levine–Mertens–Papadoulaki–Papalini (arXiv:2411.16922) and Okuyama
  (arXiv:2505.08116) is the upper edge: `θ = |ln q^2| k` [SI §4.3]. The alternative
  (Narovlansky–Verlinde, arXiv:2310.16994 etc.) puts dS at the middle of the band,
  `θ = π/2`, i.e. `E = 0` [SI §4.3, footnote]. Which (if either) notion of "dS region" gives the
  claim content is the business of issue #4, not of this file.

## 2. The quantum universal enveloping algebra U_q(su(1,1))

Following SI §3.1 (eq. (3.3) [label `eq:quantums2lrelations`] and the two displays after it, (3.4)), which
agrees with the standard convention of Klimyk–Schmüdgen:

Generators `E, F, K, K^{-1}`, relations

    K K^{-1} = K^{-1} K = 1,   K E = q E K,   K F = q^{-1} F K,
    E F - F E = (K^2 - K^{-2}) / (q - q^{-1}),

*-structure

    K* = K,   E* = -F,

Hopf structure

    Δ(K) = K ⊗ K,   Δ(E) = K ⊗ E + E ⊗ K^{-1},   Δ(F) = K ⊗ F + F ⊗ K^{-1},
    S(K) = K^{-1},  S(E) = -q^{-1} E,  S(F) = -q F,   ε(K) = 1,  ε(E) = ε(F) = 0.

Classically `K = q^{H/2}`, `q → 1` [SI §3.1]. The classical algebra [SI §2.1]:
`[H,E] = 2E`, `[H,F] = -2F`, `[E,F] = H`, `H* = H`, `E* = -F`, classical Casimir
`Ω_cl = ¼(H^2 - 2H + 4EF)` (SI eq. (2.6), label `eq:classicalCasimir`).

## 3. The Casimir and the DSSYK transfer matrix

- Casimir [SI eq. (3.5), label `eq:CasimirElement`]:

      Ω := (q^{-1} K^{-2} + q K^2 - 2)/(q - q^{-1})^2 + F E .

- Symmetrised Casimir [SI §3.3 and §3.5; the two defining displays agree]:

      Ω̃ := ½ ( (q - q^{-1})^2 Ω + 2 ) = ½ (q K^2 + q^{-1} K^{-2}) + ½ (q - q^{-1})^2 F E .

- DSSYK transfer matrix [SI §4.2]:  `½ √(1-q^2) T = Ω̃`, i.e. `T = 2 Ω̃ / √(1-q^2)`, and on
  chord states

      Ω̃ |k⟩ = ½ √(1-q^{2k+2}) |k+1⟩ + ½ √(1-q^{2k}) |k-1⟩ ,   k ∈ Z_{≥0} ,

  whose eigenvectors are `|P^θ⟩ = Σ_k H_k(cos θ; q^2)/√((q^2;q^2)_k) |k⟩` with eigenvalue
  `cos θ` (continuous q²-Hermite; verified numerically, check (C5)).
- `μ(y) := ½(y + y^{-1})` [SI App. A.2]; `[z]_q := (q^z - q^{-z})/(q - q^{-1})` [SI eq. (A.1), label
  `eq:qnumber`]; q-Pochhammer `(z;q)_n`, q-Gamma `Γ_q` [SI eqs. (A.2), (A.3)].

## 4. The coordinate algebra and the normaliser extension ⋊ Z₂

- `A_q(SU(1,1))` [SI eq. (3.1), label `eq:algebraicrelations`]: generators α, γ,

      αγ = qγα,  αγ* = qγ*α,  γγ* = γ*γ,  αα* - q^2 γ*γ = 1 = α*α - γ*γ,
      Δ(α) = α⊗α + qγ*⊗γ,  Δ(γ) = γ⊗α + α*⊗γ,  ε(α)=1, ε(γ)=0, S(α)=α*, S(γ)=-qγ.

- By Woronowicz's no-go theorem [SI §3.1, citing W1991], `SU_q(1,1)` does **not** exist as a
  locally compact quantum group: the tensor product of the irreducible *-representations π_θ of
  A_q is ill-defined (self-adjoint extensions of Δ(γ*γ) are not preserved by the coproduct
  action).
- **The normaliser extension.** Classically, the normaliser of SU(1,1) in SL(2,C) is
  `SU(1,1) ⋊ Z_2 = SU(1,1) ∪ SU(1,1)·v` [SI §2, eq. (2.25), label `eq:definitionnormalisersu11`], matrices
  `[[α, e γ̄],[γ, e ᾱ]]` with `|α|^2 - |γ|^2 = e ∈ {±1}`. Quantum version [SI §3.1, eq. (3.14),
  label `eq:e-relations`, and the display after it (3.15)]: add a central group-like self-adjoint
  generator
  `e`, `e^2 = 1`, commuting with α, γ, with the relations deformed to

      αα* - q^2 γ*γ = e = α*α - γ*γ,
      Δ(α) = α⊗α + q e γ*⊗γ,  Δ(γ) = γ⊗α + e α*⊗γ,  S(α) = eα*,  S(γ) = -qγ,

  and dual pairing `⟨X, e⟩ = ε(X)` for X ∈ U_q. **Why the extension is needed** [SI §3.1]: on the
  extended representation space the position operator is `e γ*γ` (spectrum in `±q^{2Z}`), the
  deficiency indices of Δ(e γ*γ) become (2,2), and there exist self-adjoint extensions whose
  domains — the "smoothness conditions" `f(∞) = λ f(-∞)`, `(D_{q^2}f)(∞) = λ (D_{q^2}f)(-∞)`,
  λ ∈ T [SI eq. (3.19), label `eq:smoothnessConditionsnormalizer`] — are preserved by the coproduct action,
  making SU_q(1,1)⋊Z₂ a genuine von Neumann algebraic quantum group (Koelink–Kustermans; SI §3.2).
  SI argue moreover that the reduction to DSSYK "works exclusively at the level of the
  normaliser": the chord subsector lies at `ε = -1` (the `SU(1,1)·v` component), see §7 below.

## 5. Representation series of U_q(su(1,1))

Admissible = *-representation with K-eigenvalues in `q^{½Z}` and finite-dimensional K-eigenspaces
[GKK §4.4]. GKK (following BK and MMNNSU) list the irreducible admissible *-representations; "each
is completely determined by the eigenvalue of the Casimir operator and the spectrum of K"
[GKK §4.4, before the list]. In **our (= SI) conventions**, using the dictionary of §10
(`Ω̃ = -Ω_GKK`, verified numerically, check (C2)):

| Series | Label range | K-spectrum | Ω̃ eigenvalue | classical q→1 fate |
|---|---|---|---|---|
| principal unitary `π_{b,ε}` | `b ∈ [0, -π/(2 ln q)]`, `ε ∈ {0,½}`, `(b,ε) ≠ (0,½)` | `q^{Z+ε}` (two-sided) | `-cos(2b ln q) = cos θ`, `θ := π + 2b ln q ∈ [0,π]`; sweeps `[-1, 1]` | survives (principal series of SL(2,R)) |
| positive/negative discrete `D_ℓ^±` | `ℓ ∈ ½ Z_{>0}` | `q^{±(ℓ+Z_{≥0})}` (one-sided) | `+μ(q^{2ℓ-1}) ∈ [1, ∞)` | survives (discrete series) |
| complementary | `λ ∈ (-½, 0)`, `ε = 0` | `q^Z` (two-sided) | `+μ(q^{1+2λ}) ∈ (1, μ(q))` | survives (complementary series) |
| **strange** `π^S_{a,ε}` | `a > 0`, `ε ∈ {0,½}` | `q^{Z+ε}` (two-sided) | `-μ(q^{2a}) ∈ (-∞, -1)` | **no classical analogue** [GKK §1; SI §4.3] |

Matrix elements (GKK conventions, §10 below; GKK labels `eq:posdiscrserrep`,
`eq:negdiscrserrep`, `eq:princunitaryserrep`, `eq:strangeserrep`): e.g. the **strange series**
acts on ℓ²(Z) by

    K e_n = q^{n+ε} e_n,
    (q^{-1}-q) E_GKK e_n = q^{-n-ε-½} √( (1+q^{2n+2ε+1+2a})(1+q^{2n+2ε+1-2a}) ) e_{n+1},
    F_GKK = (E_GKK)*   [GKK *-structure],

manifestly unitary for **all** `a > 0` because of the `+` signs (for the principal series the
same formula has `-` signs and unitarity restricts `b`). The substitution generating the strange
series from the discrete-series parameter in SI is `ℓ → ℓ + iπ/(2 ln q)` [SI §4.3, citing BK],
under which `Ω̃ = μ(q^{2ℓ-1}) ↦ μ(q^{2ℓ-1} e^{iπ}) = -μ(q^{2ℓ-1})`; the correspondence of labels
is `a = ℓ - ½`. SI take `ℓ ∈ ½Z_{>0}`; note `ℓ = ½` gives `a = 0`, which is **not** in the
strange family (it is the degenerate principal-series point `Ω̃ = -1`, the θ = π band edge).
Integer ℓ ↔ ε = 0; half-odd ℓ ↔ ε = ½.

**Decidability.** Among irreducible admissible unitary *-representations, membership is decided
by (Ω̃ eigenvalue, K-spectrum) alone [GKK §4.4]:
- `Ω̃ ∈ [-1, 1]`, two-sided K-spectrum → principal;
- `Ω̃ ∈ (1, μ(q))`, two-sided → complementary; `Ω̃ ≥ 1`, one-sided → discrete;
- **`Ω̃ < -1` (two-sided K-spectrum is then automatic) → strange.**

**Warning (sign slip in GKK).** GKK state `π_{b,ε}(Ω_GKK) = μ(q^{2ib}) = cos(2b ln q)`; direct
computation from their printed matrix elements gives `Ω_GKK = -cos(2b ln q)` (numerically
confirmed to 60 digits, check (C1); the discrete and strange eigenvalues as printed *are*
reproduced). Since `b ↦ -π/(2 ln q) - b` flips the sign within the same family, no statement about
the *family* is affected; pointwise dictionaries (e.g. `θ = π + 2b ln q` above) use the
matrix-element-derived sign. Similarly we record that the complementary-series eigenvalue
(not stated in GKK) is `Ω_GKK = -μ(q^{1+2λ})`, computed the same way.

## 6. The regular representation and the Plancherel decomposition

GKK Theorems (labels `thm:decompKp+-asUqmod`, `thm:decompKp-+asUqmod`, `thm:decompKp--asUqmod`,
`thm:decompKp++asUqmod`) decompose the GNS space K of the Haar weight of
`L^∞_q(SU(1,1)⋊Z_2)` as a U_q(su(1,1))-module, sector by sector (`p ∈ q^Z`, signs (ε,η)):
each sector contains the **full principal integral** `∫_0^{-π/(2 ln q)} π_{b,ε(p)} db` with
`ε(p) = ½ χ(p) mod 1`, `χ(p) = log_q p`; discrete series appear as discrete sums (e.g.
`D^+_{l+½χ(p)}`, `2l + χ(p) > 1` in the (+,-) sector); the **strange series appear as discrete
sums** with labels `π^S` at `q^{2a} = p q^{1+2l}` (i.e. `2a = χ(p) + 1 + 2l`, so `a ∈ ½Z_{>0}`)
in the (+,-), (-,+) and (+,+) sectors. The **complementary series does not appear** (measure
zero), exactly as classically. The strange series thus enters the Plancherel formula of the
*normaliser* quantum group with discrete-series-like (counting) measure. [GKK; summarised in SI
§3.3: "the left regular representation decomposes into a direct integral of the principal unitary
series, the discrete series and the strange series".]

Caveat inherited from the sources (record, do not resolve here): SI flag that for the associated
quantum *semigroup* split-offs, "the mathematical status of the corresponding Plancherel-type
theory is not clear to us" [SI §5]; any use of Plancherel theory beyond GKK's theorems must be
declared as an assumption (CONTRIBUTING.md §2, PROBLEM.md §2).

## 7. Quantum homogeneous space objects (choice fixed by issue #5)

**The DSSYK space is the twisted-primitive double quotient** `(Y_s)\G_q/(Y_t)` of
`G_q = SU_q(1,1)⋊Z_2` (invariant algebra `⟨ρ_st, e⟩`, SI (3.47)-(3.48)), in the rescaled
regime SI (4.2): two-sided version = the chord sector (SI (4.1), (4.5)); one-sided version =
SI's reduced quantum AdS_{2,q}, the restricted q-lattice `R²_{q²}(ξ)` (SI (4.22)) with inner
product (4.24), Casimir action (4.25), classical coordinates (4.30). NOT the quantum disk
(Cartan quotient — SI §5 call it "distinct from the DSSYK coset considered in this paper"),
NOT the quantum hyperboloid. The radial coordinate is `y` (equivalently `φ`, `e^φ = ½y^{-½}`).
Everything is defined on the momentum-space lattice; the position picture would need SI's
assumed q-Fourier transform (their §4.4) and is not used here. On the chord sector the state
decomposition is principal-only (tex/main.tex, Prop. on chord-sector states); the strange
series enters as bilocal operators and as summands of L²_q(G_q) (GKK sectors). See
tex/main.tex §3.

Objects defined by SI (any construction must reproduce these):

- Representation space of the normaliser coordinate algebra: `L²_{q^2}(I_{q^2})`,
  `I_q = (-∞,-1)_q ∪ [0,∞)_q`, position operator `e γ*γ` [SI §3.1, eq. (3.18), label `eq:repsNormalizer`].
- The `(s,t)`-spherical element `ρ_{st}` (twisted-primitive-invariant element), its normalised
  version `ρ̃_{st}`, rescaled `ρ̂_{st} = q^{2r} ρ̃_{st}`, `r → ∞` [SI eqs. (3.47), (3.62), (4.13); labels
  `eq:sphericalelement`, `eq:normalizedSphericalElement`, `eq:actionsphericalelement`], with

      ρ̂_{st} |k⟩ = ½ q^{-1} q^{-2k} |k⟩   (boundary conditions η = ξ = 1).

- **Chord states** [SI §4.1, eq. (4.5), label `eq:chordstates`]: the chord subsector is spanned by
  `|k⟩_{η,ξ} = lim_{r→∞} δ^{-,(st)}_{q^r η, q^{2r} y, q^r ξ}` with `y = q η ξ q^{2k}`,
  `k ∈ Z_{≥0}`; it sits at `ε = -1` (the non-identity component of the normaliser — this is SI's
  "the reduction works exclusively at the level of the normaliser").
- Coordinates ↔ classical variables [SI §4.1, eq. (4.4), label `eq:rescaledClassicalrelation`]:
  `p_L = ½ η^{-1}/(q^{-1}-q)`, `e^φ = ½ y^{-½}`, `p_R = ½ ξ^{-1}/(q^{-1}-q)`; `2φ = 2L + ln ε_r`
  (renormalised geodesic length); chord number k ↔ renormalised length via `y = qηξ q^{2k}`,
  so the DSSYK length variable is `ℓ_length = λ k` in the triple-scaling regime (`λ = -2 ln q`).
  Note: we write `ℓ` for the discrete/strange series label (SI's usage) and `ℓ_length` for the
  bulk geodesic length; SI use ℓ for both — we deviate to avoid a clash, stated here once.
- Bilocal operators [SI §4.2, §4.3]: discrete series `O^ℓ(τ) = e^{τT} (2q ρ̂_{st})^{-ℓ} e^{-τT}`;
  strange series `S^ℓ(τ) = e^{τT} (2q ρ̂_{st})^{-ℓ - iπ/(2 ln q)} e^{-τT}`, `ℓ ∈ ½Z_{>0}`.
  Two-point overlap (discrete) [SI §4.2]:

      ⟨P^{θ1}| O^ℓ(0) |P^{θ2}⟩ = (q^{4ℓ}; q^2)_∞ / (q^{2ℓ} e^{i(±θ1±θ2)}; q^2)_∞ ,

  and for the strange series the same with `q^{2ℓ} → q^{2ℓ} e^{iπ} = -q^{2ℓ}` [SI §4.3].

## 8. AdS and dS regions (as the sources use the words)

Recorded for use by issue #4; this file only fixes what the words mean **in the cited sources**:

- "AdS region" = lower spectral edge `θ = π - |ln q^2| k₁` where the q→1 limit is the Schwarzian
  theory on AdS₂ [SI §4.2; standard DSSYK literature].
- "dS region" = upper spectral edge `θ = |ln q^2| k₂`, **under the proposal** of
  arXiv:2411.16922 and arXiv:2505.08116 (dS-JT); SI's claim is explicitly conditional on it:
  "If this proposal were proven to be correct, the strange series would therefore have a non-zero
  amplitude between the anti-de Sitter and de Sitter regions" [SI §4.3, end].
- The rival identification (Verlinde et al.): dS at `θ = π/2` (`E(θ)=0`) [SI §4.3 footnote].

## 9. Special functions

All from SI App. A (agreeing with Koekoek–Lesky–Swarttouw (KLS) and Gasper–Rahman): q-number
`[z]_q` (§3), `μ(y) = ½(y+y^{-1})`, q-Pochhammer `(z;q)_n`, q-Gamma `Γ_q(z) =
(1-q)^{1-z}(q;q)_∞/(q^z;q)_∞`, basic hypergeometric `rφs`, continuous q-Hermite `H_k(x;q)`
(recurrence `H_{k+1} = 2xH_k - (1-q^k)H_{k-1}`), Al-Salam–Chihara, Askey–Wilson
polynomials/functions. The Casimir action in the quantum Gauss decomposition is the Askey–Wilson
second-order q-difference operator with parameters `(a,b,c,d) = (qsη, qt^{-1}ξ, qtξ, qs^{-1}η)`
[SI §4.1, discussion around the chord subsector; SI §3.5, eq. (3.85), label
`eq:CasimiractionGauss`, and the `ψ(y)` displays after it]. The DSSYK limit is `a,b,c,d → 0`
(continuous q²-Hermite).

## 10. Translation table

Verified rows are marked [num] (numerically, `src/check_conventions.py`) or [alg] (checked by
hand); rows marked [2nd] are taken from the stated source without independent verification of
that source's own printing.

| Source | Generators → ours | q | Casimir → ours | Verified |
|---|---|---|---|---|
| SI arXiv:2512.10101v2 | (identity — we adopt SI) | q | Ω, Ω̃ as §3 | — |
| GKK arXiv:0905.2830v2 | `E_SI = E_GKK`, `F_SI = -F_GKK`, `K_SI = K_GKK` (GKK have `FE-EF = (K²-K^{-2})/(q-q^{-1})`, `E* = F`) | same q | `Ω̃_SI = -Ω_GKK`, `Ω_GKK = ½((q^{-1}-q)² F_GKK E_GKK - qK² - q^{-1}K^{-2})` | [num] (C1),(C2) |
| BINN arXiv:2212.13668 | `A = K`, `B = E`, `C = F`, `D = K^{-1}` | `q̂_BINN = q`; `q_SYK = q̂² = q²` | `Ω_BINN = Ω_SI` exactly | [num] (C3) |
| MMNNSU LMP 19 (1990) | via GKK: `e = E_GKK`, `f = -F_GKK`, `k = K` [GKK §4.1] | same q | `Ω_GKK = -½(q-q^{-1})² C_MMNNSU - 1` [GKK, after label `eq:Casimir`] ⟹ `C_MMNNSU = Ω_SI` | [alg] from the two GKK rows; MMNNSU printing itself [2nd] |
| BK J.Phys.A 26 (1993) | classification cited by SI and GKK; representation list as §5 | same q | series labels as §5 | [2nd] (paper not independently obtained) |
| Koelink–Stokman (Askey–Wilson functions, e.g. math/0004053) | SI §3.5 reduce to their action at `ε=+1`, `η=-s^{-1}q^{2i}`, `ξ=t^{-1}q^{2j}` [SI §3.5, end] | same q | Askey–Wilson parameters `(a,b,c,d) = (qsη, qt^{-1}ξ, qtξ, qs^{-1}η)` | [2nd] |
| Schlösser–Isachenkov arXiv:2412.19681 (classical) | classical `H,E,F` as §2; Matsuki radial-part machinery | q = 1 side | classical `Ω_cl` [SI label `eq:classicalCasimir`] | to be fixed in issue #7 when used |
| DSSYK chord literature (arXiv:1811.02584, 2212.13668) | — | `q_SYK = q²`, `q_SYK = e^{-λ}` | `T = 2Ω̃/√(1-q²)`; `E(θ) = 2cosθ/√(1-q_SYK)` | [num] (C5) for the T-recurrence |

**Known sign/printing discrepancies found (all resolved into the conventions above):**
1. GKK's stated principal-series Casimir eigenvalue `μ(q^{2ib})` disagrees by a sign with their
   own matrix elements; resolved in favour of the matrix elements, `Ω_GKK(π_{b,ε}) =
   -cos(2b ln q)` [num, 60-digit check (C1)]. No consequence at the level of series membership.
2. SI's representation `π_θ` of A_q uses a convention differing from Koelink–Kustermans by a
   basis redefinition `δ_{±q^{2k}} → e^{ikθ} δ_{±q^{2k}}`, `x → sign(p) p^{-2}` [SI §3.1,
   footnote]. We follow SI.
3. The power-of-q mismatch between quantum-group and special-function conventions
   (`A_q(SU(1,1))` representations live on `L²_{q²}`-spaces) is historical and deliberate
   [SI §3.1, footnote]; all special functions in the DSSYK sector are in base `q²`.
