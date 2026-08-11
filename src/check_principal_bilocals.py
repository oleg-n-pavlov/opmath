#!/usr/bin/env python3
"""Issue #27 (sub-issue of #13): complex-weight spherical-power bilocals and the
edge-connectivity transition.

The diagonal bilocals (2q rho_st)^{-w}, Re w > 0, have chord weights t^k with
t = q^{2w}, |t| < 1, and equal-time kernels = the q-Mehler kernel at complex t
(Lemma `lem:mehler` extends to complex t):

    K_q(t; th1, th2) = (t^2;q^2)_inf / prod_{s1,s2} (t e^{i(s1 th1 + s2 th2)}; q^2)_inf .

Edges: th^L(k) = pi - lam2 k, th^U(k) = lam2 k, lam2 = |ln q^2|.  Claims checked
(cited as (P1)-(P5) from tex Sec. 6):

 (P1) exact reflection at complex t:
        K_q(t; th^L(k1), th^U(k2)) = K_q(-t; th^L(k1), th^L(k2))   for all |t|<1.
      Consequently Rat(t) := K(t;L,U)/K(t;L,L) = prod_j R(t e^{i phi_j}),
      R(u) = (u;q^2)_inf/(-u;q^2)_inf, phi_j in {±lam2(k1+k2), ±lam2(k1-k2)},
      and Rat(t) Rat(-t) = 1 exactly.
 (P2) suppressed region |arg t| < pi/2: |Rat| -> 0 as q -> 1; the counting bound
      of prop:limits(iv) generalises: with c0 := cos(|arg t| + lam2(k1+k2)) > 0,
      every per-factor ratio with rho_n >= 1/2 is <= g_{c0}(1/2) < 1, where
      g_c(rho) = (1-2 rho c + rho^2)/(1+2 rho c + rho^2); measured |Rat| is also
      compared against the two-term asymptotic prediction
        ln|Rat| ~ sum_j Re[Li2(-t e^{i phi_j}) - Li2(t e^{i phi_j})]/lam2
                  + (1/2) sum_j ln|(1 - t e^{i phi_j})/(1 + t e^{i phi_j})|,
      with residual O(lam2) (rate measured).
 (P3) enhanced region |arg t| > pi/2: |Rat| -> infinity, exactly reciprocal to
      the mirror suppressed point (Rat(t) Rat(-t) = 1, checked exactly).
 (P4) the marginal line arg t = pi/2 (t = i rho): |Rat| = 1 EXACTLY at every q
      (the per-n factors of the symmetric phase pairs cancel identically):
      the cross- and same-edge amplitudes have equal modulus at the transition
      weight.
 (P5) the dictionary: on the principal-weight circle t = q^{1+2ib}
      (Delta = 1/2 + ib), the zero of ln|Rat| in b sits at b0(q) with
      b0/b_max -> 1/2 (b_max = pi/(2|ln q|)), i.e. at the weight whose
      theta-image is the band centre theta = pi/2; measured drift O(lam2).

Run: ./.venv/bin/python src/check_principal_bilocals.py
Output: results/check_principal_bilocals.out (deterministic).
"""

import mpmath as mp


def qpochinf(a, q):
    out = mp.mpmathify(1)
    aq = mp.mpmathify(a)
    for _ in range(200000):
        out *= (1 - aq)
        aq *= q
        if abs(aq) < mp.mpf(10) ** (-mp.mp.dps - 8):
            break
    return out


def kernel(t, th1, th2, q):
    q2 = q * q
    den = mp.mpmathify(1)
    for s1 in (1, -1):
        for s2 in (1, -1):
            den *= qpochinf(t * mp.e ** (1j * (s1 * th1 + s2 * th2)), q2)
    return qpochinf(t * t, q2) / den


def edges(q, k1, k2):
    lam2 = -2 * mp.log(q)
    return mp.pi - lam2 * k1, lam2 * k2  # th^L(k1), th^U(k2)


def ratio(t, q, k1, k2):
    lam2 = -2 * mp.log(q)
    thL1 = mp.pi - lam2 * k1
    thL2 = mp.pi - lam2 * k2
    thU2 = lam2 * k2
    return kernel(t, thL1, thU2, q) / kernel(t, thL1, thL2, q)


def main():
    print("check_principal_bilocals.py -- issue #27 (mpmath, 40 digits)")
    print("=" * 78)
    mp.mp.dps = 40
    ok = True
    k1, k2 = mp.mpf('0.8'), mp.mpf('1.3')

    # (P1) exact reflection at complex t, and Rat(t) Rat(-t) = 1
    print("(P1) K(t;L,U) == K(-t;L,L) at complex t; Rat(t)*Rat(-t) == 1")
    for (q, t) in [(mp.mpf('0.8'), mp.mpf('0.5') * mp.e ** (1j * mp.mpf('0.7'))),
                   (mp.mpf('0.9'), mp.mpf('0.9') * mp.e ** (1j * mp.mpf('2.2'))),
                   (mp.mpf('0.6'), mp.mpc('-0.3', '0.45'))]:
        lam2 = -2 * mp.log(q)
        thL1 = mp.pi - lam2 * k1
        thL2 = mp.pi - lam2 * k2
        thU2 = lam2 * k2
        lhs = kernel(t, thL1, thU2, q)
        rhs = kernel(-t, thL1, thL2, q)
        dev1 = abs(lhs - rhs) / abs(rhs)
        dev2 = abs(ratio(t, q, k1, k2) * ratio(-t, q, k1, k2) - 1)
        p = dev1 < mp.mpf(10) ** (-mp.mp.dps + 8) and dev2 < mp.mpf(10) ** (-mp.mp.dps + 8)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}, t={mp.nstr(t,4)}: "
              f"reflection rel dev = {mp.nstr(dev1,3)}; |Rat(t)Rat(-t)-1| = {mp.nstr(dev2,3)}")

    # (P2) suppressed region: counting bound + asymptotic prediction
    print("\n(P2) suppressed region arg t = pi/4 (|t| = q): counting bound and rate")
    alpha = mp.pi / 4
    for q in [mp.mpf('0.9'), mp.mpf('0.97'), mp.mpf('0.99')]:
        lam2 = -2 * mp.log(q)
        t = q * mp.e ** (1j * alpha)
        phis = [lam2 * (k1 + k2), lam2 * (k1 - k2), -lam2 * (k1 - k2), -lam2 * (k1 + k2)]
        # hypothesis and per-factor bound
        c0 = mp.cos(alpha + lam2 * (k1 + k2))
        hyp = c0 > 0
        g_half = (1 - c0 + mp.mpf('0.25')) / (1 + c0 + mp.mpf('0.25'))
        allfac = True
        count = 0
        prodbound = mp.mpf(1)
        prodactual = mp.mpf(1)
        for phi in phis:
            n = 0
            while True:
                rho = q * q ** (2 * n)  # |t| q^{2n}
                if rho < mp.mpf('0.5') or n > 100000:
                    break
                c = mp.cos(alpha + phi)
                fac = (1 - 2 * rho * c + rho * rho) / (1 + 2 * rho * c + rho * rho)
                allfac &= fac <= g_half + mp.mpf(10) ** (-30)
                prodactual *= fac
                prodbound *= g_half
                count += 1
                n += 1
        rat = ratio(t, q, k1, k2)
        # two-term asymptotic prediction for ln|Rat|
        pred = mp.mpf(0)
        for phi in phis:
            u = t * mp.e ** (1j * phi)
            pred += mp.re(mp.polylog(2, -u) - mp.polylog(2, u)) / lam2
            pred += mp.log(abs((1 - u) / (1 + u))) / 2
        meas = mp.log(abs(rat))
        resid = abs(meas - pred)
        p = hyp and allfac and abs(rat) < 1 and prodactual <= prodbound * (1 + mp.mpf(10) ** -25)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}: c0={mp.nstr(c0,4)} > 0; "
              f"{count} per-factor ratios all <= g(1/2)={mp.nstr(g_half,4)}; "
              f"|Rat|={mp.nstr(abs(rat),3)}; ln|Rat| meas/pred = "
              f"{mp.nstr(meas,6)}/{mp.nstr(pred,6)} (resid {mp.nstr(resid,3)})")
    # residual rate: O(lam2^2) -- the O(lam2) Euler-Maclaurin corrections cancel
    # in the symmetric combination sum_j [corr(-u_j) - corr(u_j)] (observed; the
    # first run of this check asserted O(lam2) and FAILED with measured exponent
    # ~2.07, recorded per CONTRIBUTING.md Sec. 5)
    resids = []
    for q in [mp.mpf('0.97'), mp.mpf('0.99')]:
        lam2 = -2 * mp.log(q)
        t = q * mp.e ** (1j * alpha)
        phis = [lam2 * (k1 + k2), lam2 * (k1 - k2), -lam2 * (k1 - k2), -lam2 * (k1 + k2)]
        pred = mp.mpf(0)
        for phi in phis:
            u = t * mp.e ** (1j * phi)
            pred += mp.re(mp.polylog(2, -u) - mp.polylog(2, u)) / lam2
            pred += mp.log(abs((1 - u) / (1 + u))) / 2
        resids.append((abs(mp.log(abs(ratio(t, q, k1, k2))) - pred), lam2))
    rate = (resids[0][0] / resids[1][0]) / (resids[0][1] / resids[1][1]) ** 2
    p = mp.mpf('0.4') < rate < mp.mpf('2.5')
    ok &= p
    print(f"  [{'PASS' if p else 'FAIL'}] residual scaling: (resid ratio)/(lam2 ratio)^2 = "
          f"{mp.nstr(rate,4)} (O(lam2^2) predicted ~ 1)")

    # (P3) enhanced region: reciprocal growth
    print("\n(P3) enhanced region arg t = 3pi/4: |Rat| = 1/|Rat(mirror)| -> infinity")
    for q in [mp.mpf('0.9'), mp.mpf('0.99')]:
        t = q * mp.e ** (1j * 3 * mp.pi / 4)
        r_enh = abs(ratio(t, q, k1, k2))
        r_sup = abs(ratio(-t, q, k1, k2))
        dev = abs(r_enh * r_sup - 1)
        p = r_enh > 10 and dev < mp.mpf(10) ** (-mp.mp.dps + 8)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}: |Rat| = {mp.nstr(r_enh,4)} "
              f"(mirror {mp.nstr(r_sup,3)}; |product-1| = {mp.nstr(dev,3)})")

    # (P4) the marginal line: |Rat(i rho)| == 1 exactly, every q
    print("\n(P4) marginal line arg t = pi/2: |Rat| == 1 EXACTLY (per-n pair cancellation)")
    for (q, rho) in [(mp.mpf('0.6'), mp.mpf('0.6')), (mp.mpf('0.9'), mp.mpf('0.9')),
                     (mp.mpf('0.97'), mp.mpf('0.5'))]:
        t = mp.mpc(0, 1) * rho
        dev = abs(abs(ratio(t, q, k1, k2)) - 1)
        p = dev < mp.mpf(10) ** (-mp.mp.dps + 8)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}, t=i*{mp.nstr(rho,2)}: "
              f"| |Rat| - 1 | = {mp.nstr(dev,3)}")

    # (P5) crossing location on the principal circle t = q^{1+2ib}: EXACT midpoint.
    # At b = b_max/2 the weight is t = q e^{-i pi/2} = -iq, i.e. ON the marginal
    # line of (P4), so |Rat| = 1 exactly at every q: the transition sits at
    # Delta = 1/2 + i b_max/2, theta-image = pi/2 (band centre), EXACTLY.
    # (The first run of this check searched for an asymptotic drift b0/b_max ->
    # 1/2 and FAILED to find one: measured drift ~ 1e-19 = bisection precision
    # at every q -- the crossing is exact, not asymptotic. Recorded per
    # CONTRIBUTING.md Sec. 5; the check now asserts exactness.)
    print("\n(P5) principal circle t = q^(1+2ib): the crossing is EXACTLY at b_max/2")
    for q in [mp.mpf('0.9'), mp.mpf('0.97'), mp.mpf('0.99')]:
        lnq = mp.log(q)
        bmax = mp.pi / (-2 * lnq)

        def F(b):
            t = q * mp.e ** (2j * b * lnq)
            return mp.log(abs(ratio(t, q, k1, k2)))

        fmid = abs(F(bmax / 2))
        fbelow = F(mp.mpf('0.4') * bmax)
        fabove = F(mp.mpf('0.6') * bmax)
        p = fmid < mp.mpf(10) ** (-mp.mp.dps + 8) and fbelow < 0 < fabove
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}: |ln|Rat|| at b_max/2 = "
              f"{mp.nstr(fmid,3)} (exact zero); sign change across: "
              f"F(0.4 b_max) = {mp.nstr(fbelow,4)} < 0 < F(0.6 b_max) = {mp.nstr(fabove,4)}")

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
