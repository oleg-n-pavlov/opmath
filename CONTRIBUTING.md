# Working protocol

This repository is a research record, not a codebase. Its purpose is that a specialist — in the
first instance M. Isachenkov, whose open question this is — can read it and check every step.
The protocol below exists to make that possible.

## 1. Issues

- **One issue = one unit of work.** Nothing is worked on that does not have an issue.
- Before starting: comment on the issue saying it is in progress and what the plan is.
- On finishing: **commit first**, then close the issue with a comment linking the commit SHA and
  stating the outcome in one paragraph.
- **When a subtask appears inside a task, open a new sub-issue for it.** Do not solve it inline
  and do not silently widen the parent issue. Attach it to the correct parent via the sub-issue
  hierarchy, not a markdown checklist.
- An issue that could not be completed stays **open**, with a comment recording what was tried and
  where it broke. Closing an issue by narrowing its scope is prohibited; narrow the scope only by
  editing the issue body and saying so explicitly in a comment.

## 2. Gaps

- Every gap in an argument is recorded as its own issue labelled `gap`, and is referenced from the
  point in `tex/main.tex` where it occurs.
- Phrases of the form "it is easily seen that", "one can show", "clearly", and "it follows that"
  are not permitted as substitutes for an argument. Either give the argument or open a `gap` issue
  and cite it.

## 3. Commits

- Commit after each closed issue, not in one batch at the end.
- Commit messages reference the issue: `Closes #7`, or `Refs #7` for partial work.
- Numerical output, figures and logs are committed alongside the code that produced them.

## 4. Formality

All mathematical content lives in `tex/main.tex` and is written as formally as the material allows:

- numbered `definition`, `lemma`, `proposition`, `theorem`, `corollary`, `remark` environments;
- every statement carries its hypotheses explicitly — no hypotheses inherited silently from
  surrounding prose;
- conventions are fixed once in `CONVENTIONS.md` and in the conventions section of `main.tex`, and
  every deviation from a cited source's conventions is stated at the point of citation;
- every statement is tagged with exactly one epistemic status:

  | Tag | Meaning |
  |---|---|
  | `PROVED` | Proved here, in full, with no appeal to an unproved lemma. |
  | `DERIVED` | Obtained by a calculation that is complete but has not been independently checked. |
  | `ARGUED` | Supported by a plausibility argument, a limit, or a special case; not a proof. |
  | `CONJECTURED` | Believed, with reasons stated; no derivation. |
  | `GAP` | A step that is required and missing. Must have a corresponding `gap` issue. |

- Results quoted from the literature are cited precisely: arXiv number, version, equation number.
  Where a source is quoted in support of the claim under test, quote it **verbatim**.

## 5. Honesty requirements

- Failed checks are reported in `LOG.md` and in the verification section of `main.tex`, with the
  same prominence as successful ones.
- Numerical agreement is reported with the actual numbers and the tolerance used, never as
  "agrees well".
- If the project reaches a negative result, that is a result: state it plainly and stop, rather
  than reformulating the problem until something survives.
- Do not simplify the problem in order to make progress. If a simplification is genuinely
  necessary, open an issue for it, record what was given up, and mark all downstream statements as
  conditional on it.

## 6. Resources

There are no restrictions on sources or tools. Read any paper (arXiv full text at
`arxiv.org/abs/<id>` and `arxiv.org/html/<id>`, INSPIRE-HEP at
`inspirehep.net/api/literature?q=...` for citation lookups, journal pages, textbooks, lecture
notes). Install any package. Run any numerical experiment.
