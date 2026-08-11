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
