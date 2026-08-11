#!/usr/bin/env python3
"""Issue #3: numerical verification of the conventions fixed in CONVENTIONS.md.

Claims tested (all statements refer to CONVENTIONS.md and its sources):

 (C1) On each irreducible admissible *-representation of U_q(su(1,1)) as listed by
      Groenevelt-Koelink-Kustermans (GKK, arXiv:0905.2830v2, Sec. 4.4 "The decomposition of the
      GNS-space K as a U_q-module", eqs. (posdiscrserrep), (negdiscrserrep), (princunitaryserrep),
      (strangeserrep)), the Casimir Omega_GKK (GKK eq. (Casimir)) computed from the printed matrix
      elements is scalar, and equals
        D_k^+  : -mu(q^{1-2k})     [GKK's stated value; CHECK]
        D_k^-  : -mu(q^{1-2k})     [GKK's stated value; CHECK]
        pi_{b,eps}   : -cos(2 b ln q)  [GKK state +mu(q^{2ib}) = cos(2b ln q); our derived value
                                        has the opposite sign; this script decides numerically]
        pi^S_{a,eps} : +mu(q^{2a})     [GKK's stated value; CHECK]
        complementary(lambda): -mu(q^{1+2lambda})  [derived here; GKK do not state it]
 (C2) The generator map E_SI = E_GKK, F_SI = -F_GKK, K_SI = K_GKK sends the GKK relations/*
      -structure to the Schouten-Isachenkov (SI, arXiv:2512.10101v2) relations/*-structure, and
      Omega~_SI = -Omega_GKK as operators (checked on the principal series, truncated, interior).
 (C3) Omega_BINN (arXiv:2212.13668, Sec. 3.1, with A=K, B=E, C=F, D=K^{-1}, qhat=q) equals
      Omega_SI (arXiv:2512.10101v2 eq. (CasimirElement)) as operators (same truncated check).
 (C4) Omega~_SI = (1/2)((q-q^{-1})^2 Omega_SI + 2) on the same truncation.
 (C5) The SI chord-sector Casimir action Omega~|k> = (1/2)sqrt(1-q^{2k+2})|k+1> +
      (1/2)sqrt(1-q^{2k})|k-1> (SI, Sec. "The partition function and n-point functions") is
      diagonalised by v_k(theta) = H_k(cos theta; q^2)/sqrt((q^2;q^2)_k) with eigenvalue cos theta
      (continuous q-Hermite recurrence), checked pointwise.

Run:  python3 src/check_conventions.py
Output: results/check_conventions.out (committed).  Deterministic, no seeds needed.
"""

import numpy as np
import mpmath as mp

mp.mp.dps = 60  # (C1) involves exact cancellations between terms of size q^{-2N}; use 60 digits

TOL = 1e-12


def mu(y):
    return 0.5 * (y + 1.0 / y)


def report(name, maxdev, tol=TOL):
    status = "PASS" if maxdev < tol else "FAIL"
    print(f"[{status}] {name}: max deviation {maxdev:.3e} (tol {tol:.1e})")
    return maxdev < tol


# ----------------------------------------------------------------------------------------------
# GKK representations: build K (diagonal), E (one off-diagonal), F = E^dagger on a truncation,
# in GKK conventions:  (q^{-1}-q) E e_n = c_n e_{n+1} (or e_{n-1} for D^-), F = E^*.
# Casimir Omega_GKK = 1/2 ((q^{-1}-q)^2 F E - q K^2 - q^{-1} K^{-2}).   [GKK eq. (Casimir)]
# The diagonal of FE is |c_n|^2, exact for every interior n (no truncation error on the
# diagonal since FE e_n = |c_n|^2 e_n for these one-sided-shift operators).
# ----------------------------------------------------------------------------------------------

def gkk_casimir_diagonal_mp(q, weights, cE2):
    """Omega_GKK diagonal values (mpmath), from K-weights q^{w_n} and |c_n|^2 with
    (q^{-1}-q) E e_n = c_n e_{n+1}. Exact per site; 60-digit arithmetic absorbs the
    cancellation between terms of size q^{-2n-1}."""
    out = []
    for w, t2 in zip(weights, cE2):
        K2 = mp.mpf(q) ** (2 * w)
        out.append(0.5 * (t2 - mp.mpf(q) * K2 - 1 / (mp.mpf(q) * K2)))
    return out


def run_C1(q):
    ok = True
    N = 15
    qm = mp.mpf(q)
    print(f"\n== (C1) GKK Casimir eigenvalues from printed matrix elements, q = {q} "
          f"(mpmath, {mp.mp.dps} digits, sites |n| <= {N})")

    def maxdev(vals, target):
        return float(max(abs(v - target) for v in vals))

    # positive discrete series D_k^+, k in (1/2)N; K e_n = q^{k+n}, n >= 0
    for k in [0.5, 1.0, 1.5, 3.0]:
        ns = range(0, 2 * N)
        w = [k + n for n in ns]
        cE2 = [qm ** (-1 - 2 * k - 2 * n) * (1 - qm ** (2 * n + 2)) * (1 - qm ** (4 * k + 2 * n))
               for n in ns]
        om = gkk_casimir_diagonal_mp(q, w, cE2)
        dev = maxdev(om, -mu(qm ** (1 - 2 * k)))
        ok &= report(f"D^+_{k}: Omega_GKK == -mu(q^(1-2k))", dev)

    # negative discrete series D_k^-: K e_n = q^{-k-n}; E lowers: (q^{-1}-q) E e_n = c_n e_{n-1}.
    for k in [0.5, 1.0, 2.5]:
        ns = range(1, 2 * N)
        w = [-k - n for n in ns]
        cE2 = [qm ** (1 - 2 * k - 2 * n) * (1 - qm ** (2 * n)) * (1 - qm ** (4 * k + 2 * n - 2))
               for n in ns]
        om = gkk_casimir_diagonal_mp(q, w, cE2)
        dev = maxdev(om, -mu(qm ** (1 - 2 * k)))
        ok &= report(f"D^-_{k}: Omega_GKK == -mu(q^(1-2k))", dev)

    # principal unitary series pi_{b,eps}: K e_n = q^{n+eps}, n in Z
    for b, eps in [(0.3, 0.0), (0.7, 0.5), (1.1, 0.0)]:
        bmax = float(-mp.pi / (2 * mp.log(qm)))
        assert 0 <= b <= bmax, (b, bmax)
        ns = range(-N, N)
        w = [n + eps for n in ns]
        cE2 = []
        for n in ns:
            u = qm ** (2 * n + 1 + 2 * eps)
            cE2.append(qm ** (-1 - 2 * n - 2 * eps) *
                       (1 - 2 * u * mp.cos(2 * b * mp.log(qm)) + u ** 2))
        om = gkk_casimir_diagonal_mp(q, w, cE2)
        dev_stated = maxdev(om, mp.cos(2 * b * mp.log(qm)))    # GKK's stated value
        dev_derived = maxdev(om, -mp.cos(2 * b * mp.log(qm)))  # our derived value
        print(f"    pi_(b={b},eps={eps}): |Omega - cos(2b ln q)| = {dev_stated:.3e}, "
              f"|Omega + cos(2b ln q)| = {dev_derived:.3e}")
        ok &= report(f"pi_(b={b},eps={eps}): Omega_GKK == -cos(2 b ln q) (sign OPPOSITE to "
                     f"GKK's stated value)", dev_derived)

    # strange series pi^S_{a,eps}: K e_n = q^{n+eps}, n in Z
    for a, eps in [(0.5, 0.0), (1.0, 0.5), (2.5, 0.0), (0.25, 0.5)]:
        ns = range(-N, N)
        w = [n + eps for n in ns]
        cE2 = [qm ** (-1 - 2 * n - 2 * eps) *
               (1 + qm ** (2 * n + 2 * eps + 1 + 2 * a)) *
               (1 + qm ** (2 * n + 2 * eps + 1 - 2 * a)) for n in ns]
        om = gkk_casimir_diagonal_mp(q, w, cE2)
        dev = maxdev(om, mu(qm ** (2 * a)))
        ok &= report(f"pi^S_(a={a},eps={eps}): Omega_GKK == +mu(q^(2a))", dev)

    # complementary series: principal with eps=0 and 2ib -> 2*lambda+1, -1/2 < lambda < 0
    for lam in [-0.25, -0.05, -0.45]:
        ns = range(-N, N)
        w = [n + 0.0 for n in ns]
        cE2 = []
        for n in ns:
            p1 = (1 - qm ** (2 * n + 2 + 2 * lam)) * (1 - qm ** (2 * n - 2 * lam))
            assert p1 > 0
            cE2.append(qm ** (-1 - 2 * n) * p1)
        om = gkk_casimir_diagonal_mp(q, w, cE2)
        dev = maxdev(om, -mu(qm ** (1 + 2 * lam)))
        ok &= report(f"complementary(lambda={lam}): Omega_GKK == -mu(q^(1+2lambda))", dev)
    return ok


# ----------------------------------------------------------------------------------------------
# (C2)-(C4): operator identities on a truncated principal-series module.
# Truncation: indices n in [-N, N]; identities checked on the interior [-N+2, N-2] where the
# truncated products agree with the untruncated ones (band width of all operators <= 2).
# ----------------------------------------------------------------------------------------------

def run_C234(q):
    ok = True
    N, b, eps = 40, 0.4, 0.0
    n = np.arange(-N, N + 1)
    dim = len(n)
    print(f"\n== (C2)-(C4) operator identities on truncated principal series "
          f"(N={N}, b={b}, eps={eps}, q={q})")

    K = np.diag(q ** (n + eps))
    Kinv = np.diag(q ** -(n + eps))
    prod = (1 - q ** (2 * n + 1 + 2 * eps) * np.exp(2j * b * np.log(q))) * \
           (1 - q ** (2 * n + 1 + 2 * eps) * np.exp(-2j * b * np.log(q)))
    c = q ** (-0.5 - n - eps) * np.sqrt(prod.real) / (q ** -1 - q)
    E_gkk = np.zeros((dim, dim))
    for i in range(dim - 1):
        E_gkk[i + 1, i] = c[i]          # E e_n = c_n e_{n+1}
    F_gkk = E_gkk.T                     # GKK *-structure: E^* = F

    # SI generators via the map E_SI = E_GKK, F_SI = -F_GKK, K_SI = K_GKK
    E_si, F_si, K_si, Kinv_si = E_gkk, -F_gkk, K, Kinv

    interior = slice(2, dim - 2)
    scale = max(np.max(np.abs(E_gkk)), np.max(np.abs(K)), np.max(np.abs(Kinv))) ** 2

    def dev_int(X):
        """Deviation on the interior block, relative to the largest entry scale of the
        operators involved (entries grow like q^{-2N}; the identities cancel exactly)."""
        return np.max(np.abs(X[interior, interior])) / scale

    # SI relations (arXiv:2512.10101v2, eq. (quantums2lrelations))
    ok &= report("SI relation K E = q E K", dev_int(K_si @ E_si - q * E_si @ K_si))
    ok &= report("SI relation K F = q^{-1} F K", dev_int(K_si @ F_si - q ** -1 * F_si @ K_si))
    comm = E_si @ F_si - F_si @ E_si - (K_si @ K_si - Kinv_si @ Kinv_si) / (q - q ** -1)
    ok &= report("SI relation EF - FE = (K^2-K^{-2})/(q-q^{-1})", dev_int(comm))
    ok &= report("SI *-structure E^* = -F", dev_int(E_si.T.conj() + F_si))

    Om_si = (q ** -1 * Kinv_si @ Kinv_si + q * K_si @ K_si - 2 * np.eye(dim)) / (q - q ** -1) ** 2 \
        + F_si @ E_si
    Om_tilde_si = 0.5 * ((q - q ** -1) ** 2 * Om_si + 2 * np.eye(dim))
    Om_gkk = 0.5 * ((q ** -1 - q) ** 2 * F_gkk @ E_gkk - q * K @ K - q ** -1 * Kinv @ Kinv)

    ok &= report("(C2) Omega~_SI == -Omega_GKK", dev_int(Om_tilde_si + Om_gkk))
    ok &= report("(C4) Omega~_SI == (1/2)((q-q^{-1})^2 Omega_SI + 2)",
                 dev_int(Om_tilde_si - 0.5 * ((q - q ** -1) ** 2 * Om_si + 2 * np.eye(dim))))

    # BINN (arXiv:2212.13668) Casimir with A=K, B=E, C=F, D=K^{-1}, qhat = q:
    # Omega_BINN = (qhat^{-1} A^2 + qhat D^2 - 2)/(qhat^{-1}-qhat)^2 + B C
    Om_binn = (q ** -1 * K_si @ K_si + q * Kinv_si @ Kinv_si - 2 * np.eye(dim)) / (q ** -1 - q) ** 2 \
        + E_si @ F_si
    ok &= report("(C3) Omega_BINN == Omega_SI", dev_int(Om_binn - Om_si))
    return ok


# ----------------------------------------------------------------------------------------------
# (C5): continuous q-Hermite diagonalisation of the chord-sector Casimir action.
# ----------------------------------------------------------------------------------------------

def cont_q_hermite(kmax, x, q):
    """H_k(x; q) via H_{k+1} = 2 x H_k - (1 - q^k) H_{k-1}, H_0 = 1, H_1 = 2x."""
    H = [1.0, 2 * x]
    for k in range(1, kmax):
        H.append(2 * x * H[k] - (1 - q ** k) * H[k - 1])
    return np.array(H[: kmax + 1])


def run_C5(q):
    ok = True
    print(f"\n== (C5) chord-sector Casimir action diagonalised by continuous q^2-Hermite, q={q}")
    q2 = q * q
    kmax = 80
    for theta in [0.3, 1.2, np.pi / 2, 2.6]:
        x = np.cos(theta)
        H = cont_q_hermite(kmax, x, q2)
        pochh = np.cumprod(np.concatenate(([1.0], 1 - q2 ** np.arange(1, kmax + 1))))
        v = H / np.sqrt(pochh)
        # (Omega~ v)_k = (1/2) sqrt(1-q^{2k}) v_{k-1} + (1/2) sqrt(1-q^{2k+2}) v_{k+1}
        k = np.arange(1, kmax - 1)
        lhs = 0.5 * np.sqrt(1 - q2 ** k) * v[k - 1] + 0.5 * np.sqrt(1 - q2 ** (k + 1)) * v[k + 1]
        # relative deviation, normalised by the local sup of |v| (at x = cos(pi/2) ~ 0 the odd
        # components v_k vanish identically; dividing by |v_k| alone would be 0/0 noise)
        local = np.maximum.reduce([np.abs(v[k - 1]), np.abs(v[k]), np.abs(v[k + 1])])
        dev = np.max(np.abs(lhs - x * v[k]) / np.maximum(local, 1e-300))
        ok &= report(f"theta={theta:.3f}: Omega~ v = cos(theta) v (relative)", dev, tol=1e-9)
    return ok


def main():
    print("check_conventions.py -- issue #3 verification")
    print("=" * 78)
    allok = True
    for q in [0.3, 0.7, 0.95]:
        allok &= run_C1(q)
    allok &= run_C234(0.6)
    allok &= run_C5(0.7)
    print("\n" + ("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED -- see above"))
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
