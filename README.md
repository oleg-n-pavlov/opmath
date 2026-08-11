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

Not started. Work is tracked in the issue tree under the epic issue #1.
