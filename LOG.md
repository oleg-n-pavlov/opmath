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
