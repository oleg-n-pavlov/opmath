#!/usr/bin/env python3
"""Issue #4 (backbone for #9/#10): numerical check of Lemmas 1.x and Proposition 1.x of
tex/main.tex Section 1 -- the strange/discrete bilocal overlap kernels.

Checked here:
 (K1) q-Mehler: chord sum  sum_k t^k/(q^2;q^2)_k H_k(x1;q^2)H_k(x2;q^2)  equals the closed form
      (t^2;q^2)_inf / prod_{4 signs} (t e^{i(+-th1+-th2)}; q^2)_inf, for t = +q^{2l} (discrete
      kernel K^O) and t = -q^{2l} (strange kernel K^S).  [Lemma q-Mehler form]
 (K2) Reflection identity: K^S(th1, th2) == K^O(th1, pi - th2) exactly at finite q.
      [Lemma reflection]
 (K3) Limit behaviour along q -> 1 (Proposition, parts (ii)/(iv)):
      K^S(L(k1),L(k2)) / K^O(L(k1),L(k2)) -> 0   and equals  K^O(L(k1),U(k2))/K^O(L(k1),L(k2)),
      with L(k) = pi - |ln q^2| k, U(k) = |ln q^2| k.
 (K4) Limit value (Proposition, parts (i)/(iii)):
      K^S(L(k1),U(k2)) * (1-q^2)^{3-2l} (q^2;q^2)_inf^3  ->  Gamma(l+-ik1+-ik2)/Gamma(2l)
      (product over 4 signs), and the SAME renormalised limit for K^O(L(k1),L(k2)).

Run: python3 src/check_claim_kernels.py     (deterministic; output: results/check_claim_kernels.out)
"""

import mpmath as mp

mp.mp.dps = 40


def qpoch(a, q, nterms=None):
    """(a; q)_infinity by direct product (|q|<1). nterms must reach |a q^n| below the working
    precision -- n > (dps+10) ln 10 / |ln q| -- otherwise the truncated tail contributes a
    relative error ~ |a| q^n / (1-q) (this bit us: see LOG.md, issue #4)."""
    if nterms is None:
        nterms = max(200, int((mp.mp.dps + 10) * mp.log(10) / max(1e-9, -mp.log(abs(q)))) + 10)
    out = mp.mpf(1)
    aq = mp.mpc(a)
    for _ in range(nterms):
        out *= (1 - aq)
        aq *= q
        if abs(aq) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
    return out


def kernel_closed(t, th1, th2, q2):
    """(t^2;q^2)_inf / prod_{4 signs}(t e^{i(+-th1+-th2)};q^2)_inf."""
    num = qpoch(t * t, q2)
    den = mp.mpf(1)
    for s1 in (1, -1):
        for s2 in (1, -1):
            den *= qpoch(t * mp.e ** (1j * (s1 * th1 + s2 * th2)), q2)
    return num / den


def kernel_chordsum(t, th1, th2, q2, kmax=4000):
    H1p, H1 = mp.mpf(0), mp.mpf(1)  # H_{-1}, H_0 at x1
    H2p, H2 = mp.mpf(0), mp.mpf(1)
    x1, x2 = mp.cos(th1), mp.cos(th2)
    tot = mp.mpf(0)
    tk = mp.mpf(1)     # t^k
    poch = mp.mpf(1)   # (q^2;q^2)_k
    for k in range(kmax):
        term = tk / poch * H1 * H2
        tot += term
        if k > 20 and abs(term) < mp.mpf(10) ** (-mp.mp.dps - 3) * max(1, abs(tot)):
            break
        # advance recurrences: H_{k+1} = 2x H_k - (1-q^{2k}) H_{k-1}
        H1, H1p = 2 * x1 * H1 - (1 - q2 ** k) * H1p, H1
        H2, H2p = 2 * x2 * H2 - (1 - q2 ** k) * H2p, H2
        tk *= t
        poch *= (1 - q2 ** (k + 1))
    return tot


def gammaprod(l, k1, k2):
    out = mp.mpf(1)
    for s1 in (1, -1):
        for s2 in (1, -1):
            out *= mp.gamma(l + 1j * s1 * k1 + 1j * s2 * k2)
    return out / mp.gamma(2 * l)


def main():
    ok = True
    print("check_claim_kernels.py -- issue #4 kernel checks (mpmath, 40 digits)")
    print("=" * 78)

    # (K1)+(K2) at moderate q
    for q in [mp.mpf('0.6'), mp.mpf('0.9')]:
        q2 = q * q
        for l in [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)]:
            for th1, th2 in [(mp.mpf('0.7'), mp.mpf('2.1')), (mp.mpf('2.9'), mp.mpf('0.3'))]:
                t = q ** (2 * l)
                dev1 = abs(kernel_chordsum(t, th1, th2, q2) - kernel_closed(t, th1, th2, q2))
                dev2 = abs(kernel_chordsum(-t, th1, th2, q2) - kernel_closed(-t, th1, th2, q2))
                dev3 = abs(kernel_closed(-t, th1, th2, q2) - kernel_closed(t, th1, mp.pi - th2, q2))
                for name, d in [("K1 discrete (chord sum = closed form)", dev1),
                                ("K1 strange  (chord sum = closed form)", dev2),
                                ("K2 reflection K^S(th1,th2)=K^O(th1,pi-th2)", dev3)]:
                    passed = d < mp.mpf(10) ** (-25)
                    ok &= passed
                    print(f"[{'PASS' if passed else 'FAIL'}] q={float(q):.2f} l={float(l):.1f} "
                          f"th=({float(th1):.2f},{float(th2):.2f}) {name}: dev={mp.nstr(d, 3)}")

    # (K3)+(K4): q -> 1 along a sequence
    l, k1, k2 = mp.mpf(1), mp.mpf('0.8'), mp.mpf('1.3')
    target = gammaprod(l, k1, k2)
    print(f"\n  target Gamma(l+-ik1+-ik2)/Gamma(2l) = {mp.nstr(target, 12)}  (l=1, k1=0.8, k2=1.3)")
    print(f"  {'q':>8} | {'ratio (ii)=(iv) ->0':>22} | {'renorm K^S(L,U) -> target':>28} | "
          f"{'renorm K^O(L,L) -> target':>28}")
    prev_ratio = None
    rows = []
    for qv in ['0.9', '0.99', '0.999', '0.9999']:
        q = mp.mpf(qv)
        q2 = q * q
        lnq2 = mp.log(q2)
        L1 = mp.pi + lnq2 * k1     # pi - |ln q^2| k1
        L2 = mp.pi + lnq2 * k2
        U2 = -lnq2 * k2            # |ln q^2| k2
        t = q ** (2 * l)
        KO_LL = kernel_closed(t, L1, L2, q2)
        KS_LL = kernel_closed(-t, L1, L2, q2)
        KS_LU = kernel_closed(-t, L1, U2, q2)
        KO_LU = kernel_closed(t, L1, U2, q2)
        ratio_ii = abs(KS_LL / KO_LL)
        ratio_iv = abs(KO_LU / KO_LL)
        # reflection identity forces ratio_ii == ratio_iv when k2-arguments correspond;
        # RELATIVE deviation (the unrenormalised kernels grow like (q^2;q^2)_inf^{-3}, i.e.
        # beyond 1e+10000 near q=1, so an absolute deviation would be meaningless)
        dev_refl = abs(KS_LU - KO_LL) / abs(KO_LL)
        renorm = (1 - q2) ** (3 - 2 * l) * qpoch(q2, q2) ** 3
        vS = KS_LU * renorm
        vO = KO_LL * renorm
        rows.append((qv, ratio_ii, ratio_iv, vS, vO, dev_refl))
        print(f"  {qv:>8} | {mp.nstr(ratio_ii, 6):>22} | {mp.nstr(vS, 10):>28} | "
              f"{mp.nstr(vO, 10):>28}   [rel |K^S(L,U)-K^O(L,L)| = {mp.nstr(dev_refl, 3)}]")
        if prev_ratio is not None:
            ok &= ratio_ii < prev_ratio  # monotone decrease towards 0
        prev_ratio = ratio_ii
        ok &= dev_refl < mp.mpf(10) ** (-25)
    # final-q closeness to target
    final_dev_S = abs(rows[-1][3] - target) / abs(target)
    final_dev_O = abs(rows[-1][4] - target) / abs(target)
    ok &= rows[-1][1] < mp.mpf('1e-3')
    ok &= final_dev_S < mp.mpf('0.02') and final_dev_O < mp.mpf('0.02')
    print(f"\n  at q=0.9999: |ratio (ii)| = {mp.nstr(rows[-1][1], 4)} (must be < 1e-3): "
          f"{'PASS' if rows[-1][1] < mp.mpf('1e-3') else 'FAIL'}")
    print(f"  relative deviation of renormalised K^S(L,U) from Gamma-target: "
          f"{mp.nstr(final_dev_S, 4)} (must be < 2e-2): "
          f"{'PASS' if final_dev_S < mp.mpf('0.02') else 'FAIL'}")
    print(f"  relative deviation of renormalised K^O(L,L) from Gamma-target: "
          f"{mp.nstr(final_dev_O, 4)} (must be < 2e-2): "
          f"{'PASS' if final_dev_O < mp.mpf('0.02') else 'FAIL'}")

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
