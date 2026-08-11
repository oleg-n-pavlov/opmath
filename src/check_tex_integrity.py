#!/usr/bin/env python3
"""Issue #22: integrity checks for tex/main.tex and tex/refs.bib.

The local TeX Live cannot compile the document (too old: mathtools/enumitem missing, tlmgr
cross-release failure; recorded in LOG.md under issue #15), so document integrity is
verified structurally:

 (T1) no duplicate \\label;
 (T2) no \\ref/\\eqref to an undefined label;
 (T3) every \\cite key exists in tex/refs.bib;
 (T4) every refs.bib key is cited somewhere (the bib policy: entries are added as sources
      are actually read and used);
 (T5) every \\begin{env} has a matching \\end{env} for the theorem-like and structural
      environments used here;
 (T6) none of the forbidden phrases of CONTRIBUTING.md Sec. 2 ("it is easily seen that",
      "one can show", "clearly", "it follows that") occurs outside comments;
 (T7) every src/... or results/... path referenced from the tex exists in the repository.

Run: python3 src/check_tex_integrity.py   (from the repository root; exit 0 iff all pass)
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def strip_comments(src):
    """Remove LaTeX comments (unescaped % to end of line)."""
    return re.sub(r'(?<!\\)%.*', '', src)


def main():
    ok = True
    with open(os.path.join(ROOT, 'tex', 'main.tex')) as f:
        raw = f.read()
    src = strip_comments(raw)
    with open(os.path.join(ROOT, 'tex', 'refs.bib')) as f:
        bib = f.read()

    def report(name, bad):
        passed = not bad
        print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {bad}" if bad else ""))
        return passed

    # (T1)/(T2)
    labels = re.findall(r'\\label\{([^}]*)\}', src)
    refs = re.findall(r'\\(?:ref|eqref)\{([^}]*)\}', src)
    ok &= report("(T1) no duplicate labels",
                 sorted({l for l in labels if labels.count(l) > 1}))
    ok &= report("(T2) no unresolved refs", sorted(set(refs) - set(labels)))

    # (T3)/(T4)
    cites = set()
    for group in re.findall(r'\\cite\{([^}]*)\}', src):
        cites.update(k.strip() for k in group.split(','))
    bibkeys = set(re.findall(r'@\w+\{([^,]+),', bib))
    ok &= report("(T3) all cite keys in refs.bib", sorted(cites - bibkeys))
    ok &= report("(T4) no uncited bib keys", sorted(bibkeys - cites))

    # (T5)
    unbalanced = []
    for env in ['definition', 'convention', 'remark', 'example', 'theorem', 'proposition',
                'lemma', 'corollary', 'claim', 'proof', 'equation', 'align', 'enumerate',
                'itemize', 'center', 'tabular', 'quote', 'abstract', 'document', 'split']:
        nb = len(re.findall(r'\\begin\{%s\}' % env, src))
        ne = len(re.findall(r'\\end\{%s\}' % env, src))
        if nb != ne:
            unbalanced.append(f"{env}: {nb} begin / {ne} end")
    ok &= report("(T5) environments balanced", unbalanced)

    # (T6)
    hits = []
    for phrase in ['easily seen', 'one can show', 'clearly', 'it follows that']:
        for m in re.finditer(r'\b' + phrase.replace(' ', r'\s+') + r'\b', src, re.I):
            ctx = ' '.join(src[max(0, m.start() - 40):m.end() + 40].split())
            hits.append(f"'{phrase}' at ...{ctx}...")
    ok &= report("(T6) no forbidden phrases (CONTRIBUTING.md Sec. 2)", hits)

    # (T7) referenced repository paths exist (paths appear with escaped underscores)
    missing = []
    for m in re.findall(r'(?:src|results)/[A-Za-z0-9_\\.]+', src):
        path = m.replace('\\_', '_').rstrip('.')
        if not os.path.exists(os.path.join(ROOT, path)):
            missing.append(path)
    ok &= report("(T7) referenced src/ and results/ paths exist", sorted(set(missing)))

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
