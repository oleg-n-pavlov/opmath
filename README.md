# opmath — strange series and AdS↔dS in the DSSYK bulk dual

An attempt on an open problem in mathematical physics: whether the **strange series**
representations of $U_q(\mathfrak{su}(1,1))$ interpolate between the AdS$_2$ and dS$_2$ regions of
the bulk dual of the double-scaled SYK model.

The question is left open — argued but not proved — in

> K. Schouten, M. Isachenkov, *The von Neumann algebraic quantum group
> $\mathrm{SU}_q(1,1)\rtimes\mathbb{Z}_2$ and the DSSYK model*, [arXiv:2512.10101](https://arxiv.org/abs/2512.10101), v2 April 2026.

**Start here:** [`PROBLEM.md`](PROBLEM.md) — the formal problem statement, method and acceptance
criteria. **Then:** [`CONTRIBUTING.md`](CONTRIBUTING.md) — the working protocol this repository is
run under.

## Layout

| Path | Contents |
|---|---|
| `PROBLEM.md` | Formal problem statement, method, acceptance criteria |
| `CONTRIBUTING.md` | Working protocol: issues, gaps, commits, formality, honesty |
| `CONVENTIONS.md` | Single source of truth for notation and conventions |
| `lit/notes/` | One formal note per paper read, named by arXiv id |
| `tex/main.tex` | The mathematical record: numbered statements, each with an epistemic status tag |
| `tex/refs.bib` | Bibliography |
| `src/` | Numerical verification code |
| `results/` | Figures, tables, run logs |
| `LOG.md` | Chronological journal, including failed attempts |

## Epistemic status tags

Every mathematical statement in this repository carries exactly one of:
`PROVED`, `DERIVED`, `ARGUED`, `CONJECTURED`, `GAP`. See `CONTRIBUTING.md` §4 for the definitions.
Anything tagged `GAP` has a corresponding issue labelled `gap`.

## Status

**Verdict reached** (2026-08-11; sessions 1–4, issues #2–#27). In brief — the full
statement, with hypotheses and epistemic tags, is `tex/main.tex` (compiles to a 30+-page
PDF; see Sections 6 and 8 for the verdict and outlook):

- **The amplitude form of the claim holds**, sharpened to an exact finite-q identity (the
  strange kernel is the discrete kernel with one argument reflected through the band
  centre), conditional only on the edge-dS identification of arXiv:2411.16922 /
  arXiv:2505.08116 — and it **fails under the centre-dS identification** of
  arXiv:2310.16994.
- **The geometric form of the claim is refuted**: every normaliser-covariant AdS/dS
  criterion assigns the strange family and its fused discrete partners the same type
  (covariance obstruction, realised concretely in the radial problem); the orientation
  datum separates them but places the *strange* towers at the AdS edge; no quantity
  crosses as the label varies. The AdS↔dS bridge the framework supports is the
  normaliser grading Π itself.
- New structures along the way: the radial part of the Casimir on the reduced quantum
  AdS_{2,q} (two operators: chord wall + q-Liouville well with deficiency indices (1,1),
  per-extension spectra on two locked q²-grids realising the GKK fusion pairing); an
  exact edge-connectivity transition in the bilocal weight at the band-centre image; the
  q→1 dictionary (Schwarzian density of states, Bessel-K wavefunctions from the
  project's own connection formula, and the precise disappearance mechanism of the
  strange series); cross-literature identifications (q-Askey α<-1 bound states =
  strange-range values; sine-dilaton curvature flip = band centre = transition locus).

Open items, each with a live issue: #13 (non-spherical bilocals), and the `gap`-labelled
#14 (rescaled Plancherel), #16 (Π as a rescaled limit), #26 (extension selection); #23
(licence) awaits the repository owner. Work is tracked under the epic issue #1.
