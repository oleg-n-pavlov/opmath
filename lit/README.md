# Literature status as of 2026-08-11 (issue #2)

Question: has the claim under test — that the strange series representations of
$U_q(\mathfrak{su}(1,1))$ interpolate between the AdS$_2$ and dS$_2$ regions of the DSSYK bulk
dual, argued in arXiv:2512.10101 — been closed, partially closed, or refuted between v2 of that
paper (27 April 2026) and today?

## Verdict

**OPEN.** The claim is argued but not proved in arXiv:2512.10101v2 itself (the authors' outlook,
item 1 of Section "Discussion and outlook", says the strange series "in some sense interpolates"
between AdS and dS regions and that "It is currently unknown what the correct physics
interpretation of such an amplitude would be"). No paper found on arXiv, INSPIRE-HEP, or Semantic
Scholar between 27 April 2026 and 11 August 2026 addresses, proves, sharpens, or refutes the
strange-series interpolation claim. No paper found even *mentions* the strange series in the DSSYK
context other than arXiv:2512.10101 itself. This is a statement about what was searched (below),
not a proof that no such work exists; in particular one recorded talk and one announced
work-in-progress (both listed below) could produce such a paper at any time.

## The reference paper

- **arXiv:2512.10101** — K. Schouten, M. Isachenkov, *The von Neumann algebraic quantum group
  $\mathrm{SU}_q(1,1)\rtimes\mathbb{Z}_2$ and the DSSYK model*. v1: 10 Dec 2025; **v2: 27 Apr
  2026** (submission history on arxiv.org/abs/2512.10101). INSPIRE recid 3092128, citation count
  10 as of 2026-08-11.
- The interpolation claim is made in the abstract ("Lastly, we make remarks on the correlation
  function related to the strange series representation, which is argued to interpolate between
  the AdS and dS regions of our $q$-homogeneous space."), in Section 1 (item on
  `sec:reductiondsSYK`), in Section `sec:strangematter` ("Strange matter"), and in the outlook.
  Verbatim quotes with locations are collected in `lit/notes/2512.10101.md` (issue #4).
- The paper's own epistemic status for the claim: **argued**. Their outlook (source line 1639):
  "the strange series gives non-zero amplitudes between AdS and dS regions, and thus in some sense
  interpolates between them. It is currently unknown what the correct physics interpretation of
  such an amplitude would be." The dS identification it relies on is itself flagged as a proposal
  (footnote to source line 1553: "If this proposal were proven to be correct, ...").

## Papers citing arXiv:2512.10101 (all of them, INSPIRE 2026-08-11, `refersto recid 3092128`)

Ten records. For each: whether it touches the strange series (full text grepped for
"strange"; the two hits found were false positives — a bibliography title "Strange Metals" and a
colloquial "looks strange").

| arXiv | date | authors | relevant content | strange series? |
|---|---|---|---|---|
| 2608.04745 | 2026-08-05 | Mariani, Mertens, Papalini, Tappeiner | ultracold RN-dS black holes from DSSYK flat-space limit; sine-dilaton | no |
| 2607.01385 | 2026-07-01 | Griguolo, Papalini, Russo, Seminara | JT in a box as Pöschl–Teller scattering; Wilson functions | no |
| 2606.26241 | 2026-06-24 | Blommaert, Tietto, Verlinde | 3d dS observer quantization; DSSYK worldline hologram | no |
| 2605.13956 | 2026-05-13 | Aguilar-Gutierrez, Kukolj, Seitz | q-Askey deformations of DSSYK; type II$_1$/I$_\infty$ algebras | no |
| 2605.03037 | 2026-05-04 | Goto, Milekhin, Verlinde, Xu | generalized free fields in dS from 1D CFT | no |
| 2604.14387 | 2026-04-15 | Rajgadia, Xu | emergent algebras from pure states in SYK double-scaling | no |
| 2602.06113 | 2026-02-05 | Aguilar-Gutierrez | $T^2$-deformed DSSYK, stretched horizon | no (false positive) |
| 2601.09801 | 2026-01-14 | Aguilar-Gutierrez, Das, Erdmenger, Xian | chaos-integrability transition | no |
| 2512.21774 | 2025-12-25 | Griguolo, Papalini, Russo, Seminara | dilaton gravity at finite cutoff | no (false positive) |
| 2511.08743 | 2025-11-11 | van der Heijden, E. Verlinde, Xu | $U_q(\mathfrak{su}(1,1))$ inside the chord algebra; JHEP 05(2026)148, DOI 10.1007/JHEP05(2026)148 | no |

Citers of arXiv:2212.13668 (Berkooz–Isachenkov–Narayan–Narovlansky, INSPIRE recid 2618661,
55 citations) with 2026 dates are a subset of the above plus one thesis-like record without arXiv
id (L. Bossi, *Dualities and Flows: Sine-Dilaton Gravity, DSSYK and TTbar-Deformed Chiral
Yang-Mills*, 2026-02-25); none engages the strange series.

## Author checks (INSPIRE + arXiv API, 2026-08-11)

- **Schouten**: two records; nothing further on DSSYK (other 2026 paper: arXiv:2604.27189, on
  long-range integrable deformations).
- **Isachenkov**: most recent record remains arXiv:2512.10101. The work in progress cited there as
  [dGIP] (J. de Groot, M. Isachenkov, H. Posthuma, on spectral triples / non-commutative geometry
  of DSSYK; bibliography key `WPJort`) **has not appeared** as of 2026-08-11. A classical
  precursor by de Groot alone exists: arXiv:2601.22171, *Pseudo-Riemannian Spectral Triples for
  $\mathrm{SU}(1,1)$* (22 Jan 2026) — classical $\mathrm{SU}(1,1)$ only, no $q$-deformation, no
  strange series.
- **van der Heijden, Blommaert, Mertens, Narovlansky, Berkooz, Verlinde, Xu, Lin, Susskind**:
  2026 output scanned via the citation lists above and arXiv searches
  (`abs:"double-scaled SYK" AND abs:"de Sitter"`); the additional hits 2604.21014 (Verlinde et
  al., 3D near-dS gravity and the DSSYK soft mode) and 2510.13986 (dS holographic complexity from
  Krylov complexity) do not discuss the strange series.
- **Groenevelt** (whose harmonic analysis of $\mathrm{SU}_q(1,1)$, arXiv:0905.2830 with Koelink
  and Kustermans, is where the strange series enters the Plancherel decomposition): 2026 papers
  are on Askey–Wilson connection formulas (arXiv:2602.15824) and Leonard trios (arXiv:2601.15052);
  no new work on the $\mathrm{SU}_q(1,1)$ Plancherel theory found. Semantic Scholar lists no
  2025–2026 citations of arXiv:0905.2830.

## arXiv listing searches (2026-08-11)

- `all:"strange series"`, most recent 40 hits: in the 2026 window only arXiv:2603.24845 (a
  $q$-analogue of Gosper's strange *number-theoretic* series — unrelated) and arXiv:2512.10101
  itself. No new paper on strange series representations of $U_q(\mathfrak{su}(1,1))$.
- `abs:"quantum group" AND abs:"DSSYK"`: only arXiv:2512.10101.
- `all:"SU_q(1,1)"`: no 2026 hits.
- `abs:"spectral triple" AND abs:"SYK"`: no hits.

## SCGP workshop (leading indicator)

*Double Scaled Sachdev-Ye-Kitaev Model: From Gravity to Many-Body Quantum Chaos*, SCGP, 11–15 May
2026 (scgp.stonybrook.edu/archives/45389; organizers Berkooz, Dietz, Jia, Lin, Verbaarschot).
M. Isachenkov spoke on **"Non-commutative geometry of the DSSYK model"**, Friday 15 May 2026,
9:30, SCGP 102. A video exists in the SCGP video portal (event id 496, video id 7653). No
associated preprint found on arXiv as of 2026-08-11; the talk title matches the announced [dGIP]
work in progress. No slides found on the workshop page.

## Adjacent 2026 results relevant to later issues (not to the open/closed status)

- arXiv:2605.13956 (q-Askey deformations of DSSYK) — the technical neighbour named in PROBLEM.md;
  cross-check target for issue #11.
- arXiv:2511.08743 = JHEP 05(2026)148 (van der Heijden–E. Verlinde–Xu) — constructs
  $U_q(\mathfrak{su}(1,1))$ generators inside the chord algebra; decomposition of the one-particle
  chord Hilbert space into **positive discrete series** only; cross-check target for issue #11.
- arXiv:2511.03779 = JHEP 05(2026)080 (Aguilar-Gutierrez; DOI 10.1007/JHEP05(2026)080) —
  cosmological entanglement entropy from the DSSYK von Neumann algebra. Note: PROBLEM.md
  attributes this arXiv id to Jiuci Xu; INSPIRE lists the author as Sergio E. Aguilar-Gutierrez.
  (Xu's type II$_1$ paper is arXiv:2403.09021, as PROBLEM.md says.)
- arXiv:2607.01385 — Pöschl–Teller radial problem for JT with Wilson-function two-point functions;
  useful classical-limit comparison for issue #7.

## Method note

Searches used: INSPIRE-HEP API (`refersto recid 3092128`, `refersto recid 2618661`, author
queries), arXiv API (`export.arxiv.org/api/query`, queries listed above), Semantic Scholar
citations API for the non-INSPIRE math literature, and the SCGP web pages. Full texts (LaTeX
source or arXiv HTML) of all ten citing papers were downloaded and grepped for "strange". Papers
without arXiv ids (e.g. the Bossi record) were not full-text checked. Coverage of pure-mathematics
venues beyond arXiv (journal-only publications in math.QA/math.OA) relies on Semantic Scholar and
is weaker; absence of evidence there is reported as exactly that.
