#!/usr/bin/env python3
"""Issue #9: the q -> 1 limits — Schwarzian density of states and Bessel-K
wavefunctions (tex Sec. 7; checks (L1)-(L2)).

 (L1) Density of states. The chord-sector Plancherel density [SI (4.11)] is
      rho_q(theta) = (1/2pi)(q^2;q^2)_inf |(e^{2i theta};q^2)_inf|^2.  At the
      lower edge theta = pi - lam2*k (lam2 = |ln q^2|) it factorises EXACTLY
      through the q-Gamma function [SI (A.3)]:
        |(e^{2i theta};q^2)_inf|^2
          = (1-q^2)^2 (q^2;q^2)_inf^2 / |Gamma_{q^2}(2ik)|^2,
      so the normalised density rhohat(k) := 1/|Gamma_{q^2}(2ik)|^2 tends, by
      Gamma_{q^2} -> Gamma, to 2 k sinh(2 pi k)/pi — the Schwarzian density of
      states.  Checked: (a) the exact factorisation against the raw theta-
      product; (b) the limit with measured O(lam2) rates.
 (L2) Wavefunctions. From this project's own closed forms (tex Sec. 4):
      (a) the confluence of the minimal-solution kernel: for s = q^{2i stilde},
            2phi1(0,0; s^2 q^2; q^2, lam2^2 u^2/4)
              -> Gamma(1+2i stilde) (u/2)^{-2i stilde} I_{2i stilde}(u),
          (I = modified Bessel), with O(lam2) rate;
      (b) hence, via the connection formula (eq:connection), the transfer-
          matrix eigenvectors reduce to the Liouville/Schwarzian wavefunctions:
            H_k(mu(s); q^2) -> 2 K_{2i stilde}(u) / (lam2 (q^2;q^2)_inf),
          u = 2 q^{k+1}/lam2 (the renormalised-length variable e^{-ltilde/2}),
          at the upper edge; the lower (AdS) edge is the same statement times
          (-1)^k (parity).  Checked pointwise at fixed (stilde, u) with
          measured rates, in both the oscillatory (small u) and damped
          (moderate u) regions.

Run: ./.venv/bin/python src/check_q1_limits.py
Output: results/check_q1_limits.out (deterministic).
"""

import mpmath as mp


def qpochinf(a, q):
    out = mp.mpmathify(1)
    aq = mp.mpmathify(a)
    for _ in range(400000):
        out *= (1 - aq)
        aq *= q
        if abs(aq) < mp.mpf(10) ** (-mp.mp.dps - 8):
            break
    return out


def gamma_q2(x, q):
    """Gamma_{q^2}(x) = (1-q^2)^{1-x} (q^2;q^2)_inf / (q^{2x};q^2)_inf [SI (A.3)]."""
    q2 = q * q
    return (1 - q2) ** (1 - x) * qpochinf(q2, q2) / qpochinf(q ** (2 * x), q2)


def phi00(C, q, z):
    tot = mp.mpmathify(0)
    term = mp.mpmathify(1)
    n = 0
    while True:
        tot += term
        term = term * z / ((1 - q ** (n + 1)) * (1 - C * q ** n))
        n += 1
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 8) and n > 8:
            break
        if n > 100000:
            raise RuntimeError("phi00 did not converge")
    return tot


def H_seq(x, q2, kmax):
    c = [mp.mpf(1), 2 * x]
    for k in range(1, kmax):
        c.append(2 * x * c[k] - (1 - q2 ** k) * c[k - 1])
    return c


def main():
    print("check_q1_limits.py -- issue #9 (mpmath, 40 digits)")
    print("=" * 78)
    mp.mp.dps = 40
    ok = True

    # ---------------- (L1) density of states ----------------
    print("(L1) chord Plancherel density -> Schwarzian density k sinh(2 pi k)")
    for ktil in [mp.mpf('0.6'), mp.mpf('1.4')]:
        devs = []
        for q in [mp.mpf('0.9'), mp.mpf('0.99'), mp.mpf('0.999')]:
            q2 = q * q
            lam2 = -mp.log(q2)
            theta = mp.pi - lam2 * ktil
            # (a) exact factorisation: raw product vs Gamma_{q^2} form
            raw = abs(qpochinf(mp.e ** (2j * theta), q2)) ** 2
            fact = (1 - q2) ** 2 * qpochinf(q2, q2) ** 2 / abs(gamma_q2(2j * ktil, q)) ** 2
            dev_id = abs(raw - fact) / abs(fact)
            # (b) the limit
            rhohat = 1 / abs(gamma_q2(2j * ktil, q)) ** 2
            target = 2 * ktil * mp.sinh(2 * mp.pi * ktil) / mp.pi
            dev_lim = abs(rhohat - target) / target
            devs.append((dev_lim, lam2))
            p = dev_id < mp.mpf(10) ** (-mp.mp.dps + 10)
            ok &= p
            print(f"  [{'PASS' if p else 'FAIL'}] k={mp.nstr(ktil,2)}, q={mp.nstr(q,4)}: "
                  f"factorisation rel dev = {mp.nstr(dev_id,3)}; "
                  f"rhohat = {mp.nstr(rhohat,8)} vs 2k sinh(2pi k)/pi = "
                  f"{mp.nstr(target,8)} (rel dev {mp.nstr(dev_lim,3)})")
        r1 = (devs[0][0] / devs[1][0]) / (devs[0][1] / devs[1][1])
        r2 = (devs[1][0] / devs[2][0]) / (devs[1][1] / devs[2][1])
        p = mp.mpf('0.5') < r1 < 2 and mp.mpf('0.5') < r2 < 2
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}]        limit rate: (dev ratio)/(lam2 ratio) = "
              f"{mp.nstr(r1,4)}, {mp.nstr(r2,4)} (O(lam2) predicted ~ 1)")

    # ---------------- (L2a) confluence of the minimal-solution kernel ----------
    print("\n(L2a) 2phi1(0,0;q^(4i stilde) q^2;q^2, lam2^2 u^2/4) -> "
          "Gamma(1+2i stilde)(u/2)^(-2i stilde) I_(2i stilde)(u)")
    st = mp.mpf('0.7')
    for u in [mp.mpf('0.4'), mp.mpf('1.5')]:
        devs = []
        for q in [mp.mpf('0.9'), mp.mpf('0.99'), mp.mpf('0.999')]:
            q2 = q * q
            lam2 = -mp.log(q2)
            s = q ** (2j * st)
            z = lam2 ** 2 * u ** 2 / 4
            lhs = phi00(s * s * q2, q2, z)
            nu = 2j * st
            rhs = mp.gamma(1 + nu) * (u / 2) ** (-nu) * mp.besseli(nu, u)
            dev = abs(lhs - rhs) / abs(rhs)
            devs.append((dev, lam2))
        r1 = (devs[0][0] / devs[1][0]) / (devs[0][1] / devs[1][1])
        r2 = (devs[1][0] / devs[2][0]) / (devs[1][1] / devs[2][1])
        p = mp.mpf('0.5') < r1 < 2 and mp.mpf('0.5') < r2 < 2
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] stilde={mp.nstr(st,2)}, u={mp.nstr(u,2)}: "
              f"rel devs {mp.nstr(devs[0][0],3)}/{mp.nstr(devs[1][0],3)}/"
              f"{mp.nstr(devs[2][0],3)} at q=0.9/0.99/0.999; rates "
              f"{mp.nstr(r1,4)}, {mp.nstr(r2,4)} (O(lam2) ~ 1)")

    # ---------------- (L2b) H_k -> 2 K/(lam2 (q^2;q^2)_inf) ----------------
    print("\n(L2b) H_k(mu(q^(2i stilde));q^2) -> 2 K_(2i stilde)(u)/(lam2 (q^2;q^2)_inf), "
          "u = 2 q^(k+1)/lam2")
    for u0 in [mp.mpf('0.4'), mp.mpf('1.5')]:
        devs = []
        for q in [mp.mpf('0.9'), mp.mpf('0.99'), mp.mpf('0.999')]:
            q2 = q * q
            lam2 = -mp.log(q2)
            # integer k nearest to q^{k+1} = lam2 u0/2; then use the actual u
            k = int(mp.floor(mp.log(lam2 * u0 / 2) / mp.log(q) - 1 + mp.mpf('0.5')))
            u = 2 * q ** (k + 1) / lam2
            x = (q ** (2j * st) + q ** (-2j * st)) / 2  # mu(s) = cos(stilde lam2)
            Hk = H_seq(mp.re(x), q2, k + 1)[k]
            target = 2 * mp.re(mp.besselk(2j * st, u)) / (lam2 * qpochinf(q2, q2))
            dev = abs(Hk - target) / abs(target)
            devs.append((dev, lam2, k, u))
        r1 = (devs[0][0] / devs[1][0]) / (devs[0][1] / devs[1][1])
        r2 = (devs[1][0] / devs[2][0]) / (devs[1][1] / devs[2][1])
        p = mp.mpf('0.4') < r1 < mp.mpf('2.5') and mp.mpf('0.4') < r2 < mp.mpf('2.5')
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] stilde={mp.nstr(st,2)}, u0={mp.nstr(u0,2)}: "
              f"rel devs {mp.nstr(devs[0][0],3)}/{mp.nstr(devs[1][0],3)}/"
              f"{mp.nstr(devs[2][0],3)} at q=0.9/0.99/0.999 "
              f"(k = {devs[0][2]}/{devs[1][2]}/{devs[2][2]}); rates "
              f"{mp.nstr(r1,4)}, {mp.nstr(r2,4)} (O(lam2) ~ 1)")

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
