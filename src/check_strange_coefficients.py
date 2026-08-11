#!/usr/bin/env python3
"""Issue #6: radial matrix coefficients of the strange series on the chord lattice.

The radial (chord-lattice) Casimir equation at spectral value omega, in unnormalised form
(tex/main.tex Sec. 4; c_k = v_k * sqrt((q^2;q^2)_k)):

    2 omega c_k = c_{k+1} + (1 - q^{2k}) c_{k-1} ,   k >= 1,          (R)

with the k=0 boundary condition c_1 = 2 omega c_0 defining the *regular* solution
c^reg_k = H_k(omega; q^2) (continuous q^2-Hermite). Claims checked:

 (S1) closed form of the minimal solution: for omega = mu(s) = (s+1/s)/2 with 0 < s < 1,
        c^min_k = s^k * 2phi1(0,0; s^2 q^2; q^2, q^{2k+2})
      satisfies (R) for all k >= 1 (it does NOT satisfy the k=0 boundary condition unless
      omega is an eigenvalue -- it never is here; the k=0 defect is recorded).
      At the strange value omega = -mu(s) the minimal solution is (-1)^k c^min_k
      (parity map of (R)).
 (S2) growth/decay: c^min_k ~ s^k (decaying; v^min in l^2 iff sum s^{2k}/(q^2;q^2)_k < inf:
      always, since s<1); c^reg_k = H_k(mu(s);q^2) ~ s^{-k}/(s^2;q^2)_inf (growing; v^reg not
      l^2). Measured rates compared with predictions.
 (S3) Wronskian: W_k := c^reg_{k+1} c^min_k - c^min_{k+1} c^reg_k = (q^2;q^2)_k W_0,
      with W_0 = 2 omega c^min_0 - c^min_1 (the 'Jost defect'; nonzero <=> no bound state).
 (S4) The half-shift: the minimal solution at the discrete-series Casimir value
      omega = mu(q^{2l-1}) decays like q^{(2l-1)k} while the discrete bilocal weight is
      q^{2lk}; ratio q^{-k} * const -- recorded (the rho-shift), not a discrepancy.

Run: python3 src/check_strange_coefficients.py
Output: results/check_strange_coefficients.out (deterministic).
"""

import mpmath as mp

mp.mp.dps = 40


def phi_00(c, q, z, nmax=800):
    """2phi1(0,0;c;q,z) = sum_n z^n / ((q;q)_n (c;q)_n)."""
    tot = mp.mpf(0)
    term = mp.mpf(1)
    qq = mp.mpf(1)
    cc = mp.mpf(1)
    zn = mp.mpf(1)
    for n in range(nmax):
        tot += zn / (qq * cc)
        zn *= z
        qq *= (1 - q ** (n + 1))
        cc *= (1 - c * q ** n)
        if abs(zn) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
    return tot


def cmin(k, s, q2):
    return s ** k * phi_00(s * s * q2, q2, q2 ** (k + 1))


def creg_seq(omega, q2, kmax):
    c = [mp.mpf(1), 2 * omega]
    for k in range(1, kmax):
        c.append(2 * omega * c[k] - (1 - q2 ** k) * c[k - 1])
    return c


def qpochinf(a, q):
    out = mp.mpf(1)
    aq = mp.mpf(a)
    for _ in range(5000):
        out *= (1 - aq)
        aq *= q
        if abs(aq) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
    return out


def main():
    ok = True
    print("check_strange_coefficients.py -- issue #6 (mpmath, 40 digits)")
    print("=" * 78)
    for q, s in [(mp.mpf('0.6'), mp.mpf('0.35')), (mp.mpf('0.8'), mp.mpf('0.9')),
                 (mp.mpf('0.9'), mp.mpf('0.5')), (mp.mpf('0.3'), mp.mpf('0.05'))]:
        q2 = q * q
        omega = (s + 1 / s) / 2
        kmax = 60

        # (S1) recurrence for c^min at omega = +mu(s), k = 1..kmax
        devmax = mp.mpf(0)
        for k in range(1, kmax):
            lhs = 2 * omega * cmin(k, s, q2)
            rhs = cmin(k + 1, s, q2) + (1 - q2 ** k) * cmin(k - 1, s, q2)
            devmax = max(devmax, abs(lhs - rhs) / max(abs(lhs), mp.mpf(10) ** -30))
        p1 = devmax < mp.mpf(10) ** (-30)
        ok &= p1
        print(f"[{'PASS' if p1 else 'FAIL'}] q={float(q)}, s={float(s)}: (S1) minimal-solution "
              f"recurrence, max rel dev = {mp.nstr(devmax, 3)}")

        # (S1') strange value: (-1)^k c^min_k at omega = -mu(s)
        devmax = mp.mpf(0)
        for k in range(1, kmax):
            ck = (-1) ** k * cmin(k, s, q2)
            ckp = (-1) ** (k + 1) * cmin(k + 1, s, q2)
            ckm = (-1) ** (k - 1) * cmin(k - 1, s, q2)
            lhs = 2 * (-omega) * ck
            rhs = ckp + (1 - q2 ** k) * ckm
            devmax = max(devmax, abs(lhs - rhs) / max(abs(lhs), mp.mpf(10) ** -30))
        p1b = devmax < mp.mpf(10) ** (-30)
        ok &= p1b
        print(f"[{'PASS' if p1b else 'FAIL'}]            (S1') strange-value recurrence "
              f"(omega = -mu(s)), max rel dev = {mp.nstr(devmax, 3)}")

        # (S2) growth rates. Subleading corrections decay like max(s^2, q^2)^k (nearest
        # competing pole of the generating function), so the tolerance is set from that scale
        # at the comparison index K2 -- a fixed absolute tolerance at fixed k was our first,
        # WRONG attempt (see LOG.md).
        K2 = 200
        creg = creg_seq(omega, q2, K2 + 2)
        rate_reg = creg[K2] / creg[K2 - 1]
        rate_min = cmin(K2, s, q2) / cmin(K2 - 1, s, q2)
        amp_reg = creg[K2] * s ** K2 * qpochinf(s * s, q2)
        tol2 = 100 * max(s * s, q2) ** K2 + mp.mpf(10) ** (-30)
        p2 = abs(rate_reg - 1 / s) < tol2 / s and abs(rate_min - s) < tol2 \
            and abs(amp_reg - 1) < tol2 * K2
        ok &= p2
        print(f"[{'PASS' if p2 else 'FAIL'}]            (S2) rates at k={K2} "
              f"(tol {mp.nstr(tol2,2)}): c^reg ratio -> 1/s "
              f"(dev {mp.nstr(abs(rate_reg-1/s),3)}), c^min ratio -> s "
              f"(dev {mp.nstr(abs(rate_min-s),3)}); "
              f"H_k s^k (s^2;q^2)_inf -> {mp.nstr(amp_reg, 8)} (predicted -> 1)")

        # (S3) Wronskian
        W0 = 2 * omega * cmin(0, s, q2) - cmin(1, s, q2)
        devmax = mp.mpf(0)
        poch = mp.mpf(1)
        for k in range(0, kmax):
            if k > 0:
                poch *= (1 - q2 ** k)
            Wk = creg[k + 1] * cmin(k, s, q2) - cmin(k + 1, s, q2) * creg[k]
            devmax = max(devmax, abs(Wk - poch * W0) / abs(poch * W0))
        p3 = devmax < mp.mpf(10) ** (-30) and abs(W0) > mp.mpf(10) ** (-10)
        ok &= p3
        print(f"[{'PASS' if p3 else 'FAIL'}]            (S3) Wronskian W_k = (q^2;q^2)_k W_0, "
              f"max rel dev = {mp.nstr(devmax, 3)}; W_0 = {mp.nstr(W0, 8)} (nonzero: no bound "
              f"state)")

    # (S3') closed form of the Wronskian constant and (S5) the connection formula,
    # for generic (non-resonant) s: s^2 not in q^{2Z}
    print()
    for q, s in [(mp.mpf('0.6'), mp.mpf('0.35')), (mp.mpf('0.9'), mp.mpf('0.513'))]:
        q2 = q * q
        omega = (s + 1 / s) / 2
        W0 = 2 * omega * cmin(0, s, q2) - cmin(1, s, q2)
        W0_closed = (1 / s - s) / (qpochinf(s * s, q2) * qpochinf(q2, q2))
        dev = abs(W0 - W0_closed) / abs(W0)
        p = dev < mp.mpf(10) ** (-30)
        ok &= p
        print(f"[{'PASS' if p else 'FAIL'}] q={float(q)}, s={float(s)}: (S3') W_0 == "
              f"(1/s - s)/((s^2;q^2)_inf (q^2;q^2)_inf), rel dev = {mp.nstr(dev, 3)}")
        # (S5) connection: H_k(mu(s);q^2) = c^min_k(s^{-1})/(s^2;q^2)_inf
        #                                  + c^min_k(s)/(s^{-2};q^2)_inf
        K5 = 30
        creg = creg_seq(omega, q2, K5 + 2)
        devmax = mp.mpf(0)
        A1 = 1 / qpochinf(s * s, q2)
        A2 = 1 / qpochinf(1 / (s * s), q2)
        for k in range(K5):
            rhs = A1 * cmin(k, 1 / s, q2) + A2 * cmin(k, s, q2)
            devmax = max(devmax, abs(creg[k] - rhs) / abs(creg[k]))
        p = devmax < mp.mpf(10) ** (-30)
        ok &= p
        print(f"[{'PASS' if p else 'FAIL'}]            (S5) connection formula "
              f"H_k = c^min(1/s)/(s^2;q^2)_inf + c^min(s)/(s^-2;q^2)_inf, "
              f"max rel dev = {mp.nstr(devmax, 3)}")

    # (S4) the rho-shift at a discrete-series value, one sample
    q = mp.mpf('0.7'); q2 = q * q; l = mp.mpf(1)
    s = q ** (2 * l - 1)   # omega = mu(q^{2l-1})
    k = 40
    ratio = cmin(k, s, q2) / (q ** (2 * l * k))
    print(f"\n(S4) rho-shift record at q=0.7, l=1: c^min_k / q^(2lk) = "
          f"{mp.nstr(ratio, 8)} ~ q^(-k) x const (q^(-k) = {mp.nstr(q**(-k), 3)}); "
          f"const = {mp.nstr(ratio * q ** k, 8)}")

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
