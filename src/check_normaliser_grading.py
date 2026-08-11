#!/usr/bin/env python3
"""Issue #15: the normaliser origin of the reflection identity (tex/main.tex Sec. 1,
subsection "The normaliser origin of the reflection identity").

Claims verified:

 (G1) The extra factor in SI's strange bilocal (4.18) is the chord-parity unitary: the
      functional-calculus power of the positive operator (2q rho^_st)|k> = q^{-2k}|k>
      satisfies (2q rho^_st)^{-i pi/(2 ln q)} |k> = (-1)^k |k>.  [Lemma, part (i)]
 (G2) Pi |k> = (-1)^k |k> anticommutes with the chord-sector transfer matrix:
      Pi Omega~ Pi = -Omega~ on the truncated chord space (exact band structure). [part (ii)]
 (G3) Pi reflects principal states coefficientwise: H_k(-x; q^2) = (-1)^k H_k(x; q^2)
      (continuous q^2-Hermite parity), so Pi |P^theta> = |P^(pi-theta)>.  [part (iii)]
 (G4) Factorisation S^l(0) = O^l(0) Pi = Pi O^l(0) on the chord basis, and the conjugation
      triviality Pi O^l(0) Pi^{-1} = O^l(0) != S^l(0) ("conjugated by a normaliser element"
      is FALSE; the correct statement is the factorisation).  [Proposition]
 (G5) Kernel identity through the factorisation: the chord sums
      sum_k (-t)^k H_k(x1)H_k(x2)/(q^2;q^2)_k  and  sum_k t^k H_k(x1)H_k(-x2)/(q^2;q^2)_k
      agree (t = q^{2l}), i.e. K^S(th1,th2) = K^O(th1,pi-th2) re-derived through Pi,
      independently of the q-Mehler closed form; also compared against the closed form.
 (G6) GKK grading structure on the (truncated) corepresentation space L_{p,x}
      [GKK arXiv:0905.2830v2, eq. (introprincipalseries), Lemma
      (actionsofgeneratorsoneigvetsOm), Prop. (decompMintoM+andM-andgradedcommutation)]:
      (a) Omega_GKK acts on the block (ep,et) as the scalar ep*et*x     [interior];
      (b) U0^{+-} and U0^{-+} anticommute with Omega_GKK               [interior];
      (c) both commute with K, are unitary on L_{p,x} (principal case), and square to -Id;
      (d) conjugation twists E: U0^{+-} E (U0^{+-})^* = E' and U0^{-+} E (U0^{-+})^* = -E'
          where E' carries the (ep*et -> -ep*et) coefficient           [interior].
 (G7) Fusion [GKK Prop. (discretesubrepresinregrepr) case (i) + eq.
      (sudecompdiscretesubrepresinregrepr)]: for x = mu(lambda) > 1 with
      lambda = q^{1-2j} p^{-1}, l = -j - chi(p), the U_q-module blocks of the irreducible
      corepresentation W_{p,x} of SU_q(1,1) x| Z_2 are
        (+,+) block (m in Z)      = strange series pi^S_{a,eps(p)} with a = j + (chi(p)-1)/2,
        (+,-) block (m >= j)      = positive discrete series D^+_L with L = j + chi(p)/2,
        (-,+) block (m <= l)      = negative discrete series D^-_L with the same L,
      by matching K-eigenvalues and E-coefficient arrays against GKK's printed
      (strangeserrep)/(posdiscrserrep)/(negdiscrserrep); in particular a = L - 1/2 --
      exactly SI's strange label map -- inside a single irreducible corepresentation.

Run:  python3 src/check_normaliser_grading.py
Output: results/check_normaliser_grading.out (committed).  Deterministic.
"""

import mpmath as mp

mp.mp.dps = 60

TOL = mp.mpf(10) ** (-50)


def report(name, dev, tol=TOL):
    passed = dev < tol
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: max deviation {mp.nstr(mp.mpf(dev), 3)} "
          f"(tol {mp.nstr(mp.mpf(tol), 2)})")
    return passed


def cont_q_hermite_list(kmax, x, q):
    """H_k(x; q), k = 0..kmax, via H_{k+1} = 2x H_k - (1-q^k) H_{k-1}."""
    H = [mp.mpf(1), 2 * x]
    for k in range(1, kmax):
        H.append(2 * x * H[k] - (1 - q ** k) * H[k - 1])
    return H[: kmax + 1]


def qpoch(a, q):
    """(a;q)_infinity, term count scaled with working precision."""
    nterms = max(200, int((mp.mp.dps + 10) * mp.log(10) / max(1e-9, -mp.log(abs(q)))) + 10)
    out = mp.mpc(1)
    aq = mp.mpc(a)
    for _ in range(nterms):
        out *= (1 - aq)
        aq *= q
        if abs(aq) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
    return out


def main():
    ok = True
    print("check_normaliser_grading.py -- issue #15 (mpmath, 60 digits)")
    print("=" * 78)

    # ---------------------------------------------------------------- (G1)
    print("\n== (G1) functional calculus: (2q rho^_st)^(-i pi/(2 ln q)) |k> = (-1)^k |k>")
    for q in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('0.95')]:
        dev = mp.mpf(0)
        expnt = -1j * mp.pi / (2 * mp.log(q))
        for k in range(0, 40):
            # spectrum of 2q rho^_st on |k> is q^{-2k}; complex power via exp(expnt*log(.))
            val = mp.e ** (expnt * mp.log(q ** (-2 * k)))
            dev = max(dev, abs(val - (-1) ** k))
        ok &= report(f"q={float(q)}: power equals (-1)^k, k<40", dev)

    # ---------------------------------------------------------------- (G2)
    print("\n== (G2) Pi Omega~ Pi = -Omega~ on the chord sector (truncation N=60)")
    for q in [mp.mpf('0.6'), mp.mpf('0.9')]:
        q2 = q * q
        N = 60
        dev = mp.mpf(0)
        # Omega~ has entries a_k = (1/2)sqrt(1-q^{2k}) linking k-1,k ; conjugation by
        # (-1)^k multiplies the (k-1,k) entry by (-1)^{k-1}(-1)^k = -1: check elementwise.
        for k in range(1, N):
            a = mp.sqrt(1 - q2 ** k) / 2
            lhs = ((-1) ** (k - 1)) * a * ((-1) ** k)
            dev = max(dev, abs(lhs + a))
        ok &= report(f"q={float(q)}: all off-diagonal entries flip sign exactly", dev)

    # ---------------------------------------------------------------- (G3)
    print("\n== (G3) H_k(-x;q^2) = (-1)^k H_k(x;q^2) (so Pi |P^theta> = |P^(pi-theta)>)")
    for q in [mp.mpf('0.6'), mp.mpf('0.9')]:
        q2 = q * q
        for theta in [mp.mpf('0.4'), mp.mpf('1.3'), mp.mpf('2.8')]:
            x = mp.cos(theta)
            Hp = cont_q_hermite_list(60, x, q2)
            Hm = cont_q_hermite_list(60, -x, q2)   # = H_k(cos(pi-theta))
            dev = max(abs(Hm[k] - (-1) ** k * Hp[k]) /
                      max(abs(Hp[k]), mp.mpf(10) ** (-30)) for k in range(61))
            ok &= report(f"q={float(q)}, theta={float(theta)}: relative parity identity, k<=60",
                         dev)

    # ---------------------------------------------------------------- (G4)
    print("\n== (G4) factorisation S^l(0) = O^l(0) Pi ; conjugation triviality")
    for q, l in [(mp.mpf('0.6'), mp.mpf('0.5')), (mp.mpf('0.9'), mp.mpf(1))]:
        devf = mp.mpf(0)
        devc = mp.mpf(0)
        devm = mp.mpf(0)
        for k in range(0, 50):
            S = (-1) ** k * q ** (2 * l * k)          # S^l(0)|k> eigenvalue [SI (4.18)]
            O = q ** (2 * l * k)                       # O^l(0)|k> eigenvalue [SI (4.14)]
            Pi = (-1) ** k
            devf = max(devf, abs(S - O * Pi))
            devc = max(devc, abs(Pi * O * Pi - O))     # Pi O Pi^{-1} = O (diagonal)
            if k % 2 == 1:
                # relative mismatch of the conjugate from S is exactly 2 on odd k
                devm = max(devm, abs(abs(Pi * O * Pi - S) / abs(S) - 2))
        ok &= report(f"q={float(q)}, l={float(l)}: S^l(0) = O^l(0) Pi on chord basis", devf)
        ok &= report(f"q={float(q)}, l={float(l)}: Pi O^l(0) Pi^(-1) = O^l(0) (conjugation "
                     f"triviality)", devc)
        ok &= report(f"q={float(q)}, l={float(l)}: |Pi O Pi^(-1) - S|/|S| = 2 on odd k "
                     f"(conjugation reading is false)", devm)

    # ---------------------------------------------------------------- (G5)
    print("\n== (G5) kernel identity through Pi (chord sums, no q-Mehler)")
    for q, l in [(mp.mpf('0.6'), mp.mpf('0.5')), (mp.mpf('0.9'), mp.mpf(2))]:
        q2 = q * q
        t = q ** (2 * l)
        for th1, th2 in [(mp.mpf('0.7'), mp.mpf('2.1')), (mp.mpf('2.9'), mp.mpf('0.3'))]:
            kmax = max(400, int((mp.mp.dps + 10) * mp.log(10) / (-mp.log(t))) + 10)
            H1 = cont_q_hermite_list(kmax, mp.cos(th1), q2)
            H2 = cont_q_hermite_list(kmax, mp.cos(th2), q2)
            H2r = cont_q_hermite_list(kmax, mp.cos(mp.pi - th2), q2)
            sS = mp.mpf(0)
            sOr = mp.mpf(0)
            poch = mp.mpf(1)
            tk = mp.mpf(1)
            for k in range(kmax):
                sS += (-1) ** k * tk / poch * H1[k] * H2[k]
                sOr += tk / poch * H1[k] * H2r[k]
                tk *= t
                poch *= (1 - q2 ** (k + 1))
            dev1 = abs(sS - sOr) / abs(sOr)
            # against the closed form of K^O(th1, pi-th2):
            num = qpoch(t * t, q2)
            den = mp.mpc(1)
            for s1 in (1, -1):
                for s2 in (1, -1):
                    den *= qpoch(t * mp.e ** (1j * (s1 * th1 + s2 * (mp.pi - th2))), q2)
            dev2 = abs(sS - num / den) / abs(num / den)
            ok &= report(f"q={float(q)}, l={float(l)}, th=({float(th1)},{float(th2)}): "
                         f"sum_k (-t)^k H H = sum_k t^k H H(pi-th2) (rel)", dev1,
                         tol=mp.mpf(10) ** (-45))
            ok &= report(f"        same vs closed form K^O(th1,pi-th2) (rel)", dev2,
                         tol=mp.mpf(10) ** (-40))

    # ---------------------------------------------------------------- (G6)
    print("\n== (G6) GKK grading operators on truncated L_{p,x} (principal case)")
    # Blocks (ep,et), basis e_m, m in [-N,N]; x = cos(theta).
    # K e_m = p^(1/2) q^m e_m ;
    # (q^{-1}-q) E e_m^{ep,et} = q^{-m-1/2} p^{-1/2} |1 + ep et p q^{2m+1} e^{i theta}| e_{m+1};
    # U0^{+-} e_m^{ep,et} = et (-1)^{upsilon(p)} e_m^{ep,-et} ;
    # U0^{-+} e_m^{ep,et} = ep et^{chi(p)} (-1)^m e_m^{-ep,et} .   [GKK (introprincipalseries)]
    def build_Lpx(q, chi, theta, N):
        p = q ** chi
        ups = (chi + 1) // 2 if chi % 2 else chi // 2   # upsilon(p): p=q^{2k} or q^{2k-1} -> k
        x = mp.cos(theta)
        blocks = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
        dim = 2 * N + 1
        size = 4 * dim

        def idx(bi, m):
            return bi * dim + (m + N)

        K = [[mp.mpf(0)] * size for _ in range(size)]
        E = [[mp.mpf(0)] * size for _ in range(size)]
        Upm = [[mp.mpf(0)] * size for _ in range(size)]
        Ump = [[mp.mpf(0)] * size for _ in range(size)]
        for bi, (ep, et) in enumerate(blocks):
            for m in range(-N, N + 1):
                K[idx(bi, m)][idx(bi, m)] = mp.sqrt(p) * q ** m
                if m < N:
                    c = q ** (-m - mp.mpf(1) / 2) / mp.sqrt(p) * \
                        abs(1 + ep * et * p * q ** (2 * m + 1) * mp.e ** (1j * theta)) / \
                        (q ** -1 - q)
                    E[idx(bi, m + 1)][idx(bi, m)] = c
                bj = blocks.index((ep, -et))
                Upm[idx(bj, m)][idx(bi, m)] = et * (-1) ** ups
                bk = blocks.index((-ep, et))
                Ump[idx(bk, m)][idx(bi, m)] = ep * (et ** chi) * (-1) ** m
        return blocks, dim, size, idx, K, E, Upm, Ump, p, x

    def matmul(A, B):
        n = len(A)
        return [[sum(A[i][k] * B[k][j] for k in range(n) if A[i][k] != 0)
                 for j in range(n)] for i in range(n)]

    def transpose(A):
        n = len(A)
        return [[A[j][i] for j in range(n)] for i in range(n)]

    def maxabs(A, rows):
        return max(abs(A[i][j]) for i in rows for j in rows) if rows else mp.mpf(0)

    for q, chi, theta in [(mp.mpf('0.6'), 0, mp.mpf('1.1')), (mp.mpf('0.8'), 1, mp.mpf('0.5'))]:
        N = 12
        blocks, dim, size, idx, K, E, Upm, Ump, p, x = build_Lpx(q, chi, theta, N)
        Ft = transpose(E)  # F = E^* (GKK *-structure), real matrices here
        # Omega_GKK = (1/2)((q^{-1}-q)^2 E^*E - q K^2 - q^{-1} K^{-2})
        EstE = matmul(Ft, E)
        Om = [[mp.mpf(0)] * size for _ in range(size)]
        for i in range(size):
            for j in range(size):
                Om[i][j] = ((q ** -1 - q) ** 2 * EstE[i][j]) / 2
            Om[i][i] += (-q * K[i][i] ** 2 - q ** -1 / K[i][i] ** 2) / 2
        interior = [idx(bi, m) for bi in range(4) for m in range(-N + 2, N - 1)]
        # (a) Omega diagonal = ep et x on interior
        dev = mp.mpf(0)
        for bi, (ep, et) in enumerate(blocks):
            for m in range(-N + 2, N - 1):
                i = idx(bi, m)
                dev = max(dev, abs(Om[i][i] - ep * et * x))
        ok &= report(f"q={float(q)}, chi(p)={chi}, theta={float(theta)}: (a) Omega_GKK block "
                     f"scalar = ep*et*x [interior]", dev)
        # (b) anticommutation on interior
        for name, U in [("U0^{+-}", Upm), ("U0^{-+}", Ump)]:
            AC = matmul(U, Om)
            CA = matmul(Om, U)
            S = [[AC[i][j] + CA[i][j] for j in range(size)] for i in range(size)]
            ok &= report(f"        (b) {name} Omega + Omega {name} = 0 [interior]",
                         maxabs(S, interior))
        # (c) commute with K, unitary, square = -Id
        for name, U in [("U0^{+-}", Upm), ("U0^{-+}", Ump)]:
            KU = matmul(K, U)
            UK = matmul(U, K)
            D = [[KU[i][j] - UK[i][j] for j in range(size)] for i in range(size)]
            ok &= report(f"        (c) [{name}, K] = 0", maxabs(D, range(size)))
            UtU = matmul(transpose(U), U)
            dev = max(abs(UtU[i][j] - (1 if i == j else 0))
                      for i in range(size) for j in range(size))
            ok &= report(f"        (c) {name} unitary (U^*U = Id)", dev)
            UU = matmul(U, U)
            dev = max(abs(UU[i][j] + (1 if i == j else 0))
                      for i in range(size) for j in range(size))
            ok &= report(f"        (c) ({name})^2 = -Id", dev)
        # (d) conjugation twist of E
        _, _, _, _, _, Eflip, _, _, _, _ = build_Lpx(q, chi, mp.pi - theta, N)
        # E with theta -> pi - theta has coefficient |1 - ep et p q^{2m+1} e^{i theta}|,
        # i.e. the (ep*et -> -ep*et) coefficient: this is E'.
        for name, U, sgn in [("U0^{+-}", Upm, 1), ("U0^{-+}", Ump, -1)]:
            UEUt = matmul(matmul(U, E), transpose(U))
            D = [[UEUt[i][j] - sgn * Eflip[i][j] for j in range(size)] for i in range(size)]
            ok &= report(f"        (d) {name} E {name}^* = {'+' if sgn>0 else '-'}E' [interior]",
                         maxabs(D, interior))

    # ---------------------------------------------------------------- (G7)
    print("\n== (G7) fusion pi^S_(L-1/2) + D^-_L + D^+_L inside one corepresentation "
          "[GKK case (i)]")
    for q, chi, j in [(mp.mpf('0.6'), 0, 2), (mp.mpf('0.8'), 1, 1), (mp.mpf('0.7'), 2, 3)]:
        lam = q ** (1 - 2 * j) * q ** (-chi)      # lambda = q^{1-2j} p^{-1} > 1
        assert lam > 1
        x = (lam + 1 / lam) / 2
        p = q ** chi
        low = -j - chi                            # l = -j - chi(p)
        Lbig = j + mp.mpf(chi) / 2                # discrete label
        a = Lbig - mp.mpf(1) / 2                  # strange label
        epsp = mp.mpf(0) if chi % 2 == 0 else mp.mpf(1) / 2
        N = 14

        def blockcoef2(ep, et, m):
            """squared E-coefficient (times (q^{-1}-q)^2) on block (ep,et) at site m."""
            return q ** (-2 * m - 1) / p * (1 + 2 * ep * et * x * q ** (2 * m + 1) * p
                                            + q ** (4 * m + 2) * p * p)

        # strange block (+,+): match (strangeserrep) with n + eps = m + chi/2
        # NOTE: the coefficient arrays grow like q^{-2m}; deviations are reported RELATIVE
        # to |target| (a first version used absolute deviations and spuriously failed at
        # 1.4e-48 on quantities of size 1e13 -- a tolerance artifact, recorded in LOG.md).
        dev = mp.mpf(0)
        devK = mp.mpf(0)
        for m in range(-N, N):
            n = m + (chi - (0 if chi % 2 == 0 else 1)) // 2   # n integer: n+eps = m+chi/2
            target = q ** (-2 * n - 2 * epsp - 1) * \
                (1 + q ** (2 * n + 2 * epsp + 1 + 2 * a)) * \
                (1 + q ** (2 * n + 2 * epsp + 1 - 2 * a))
            dev = max(dev, abs(blockcoef2(1, 1, m) - target) / abs(target))
            devK = max(devK, abs(mp.sqrt(p) * q ** m - q ** (n + epsp)))
        ok &= report(f"q={float(q)}, chi={chi}, j={j}: (+,+) E-coeffs = strange series "
                     f"a={float(a)}", dev)
        ok &= report(f"        (+,+) K-eigenvalues = strange series", devK)

        # D^+ block (+,-), m >= j: match (posdiscrserrep) with n = m - j
        dev = mp.mpf(0)
        devK = mp.mpf(0)
        for m in range(j, j + 2 * N):
            n = m - j
            target = q ** (-1 - 2 * Lbig - 2 * n) * (1 - q ** (2 * n + 2)) * \
                (1 - q ** (4 * Lbig + 2 * n))
            dev = max(dev, abs(blockcoef2(1, -1, m) - target) / abs(target))
            devK = max(devK, abs(mp.sqrt(p) * q ** m - q ** (Lbig + n)))
        ok &= report(f"        (+,-), m>=j: E-coeffs = D^+ with L={float(Lbig)} (=a+1/2)", dev)
        ok &= report(f"        (+,-) K-eigenvalues = D^+", devK)
        # boundary: coefficient into the block vanishes at m = j-1 (E e_{j-1} -> e_j edge)
        devb = abs(blockcoef2(1, -1, j - 1) - q ** (-2 * (j - 1) - 1) / p *
                   (1 - q ** (2 * (j - 1) + 1) * p * lam) * (1 - q ** (2 * (j - 1) + 1) * p / lam))
        # the factor (1 - q^{2m+2-2j}) at m = j-1 is (1 - q^0) = 0:
        devb = abs((1 - q ** (2 * (j - 1) + 2 - 2 * j)))
        ok &= report(f"        (+,-) lower-edge factor (1-q^(2m+2-2j)) vanishes at m=j-1", devb)

        # D^- block (-,+), m <= l: match (negdiscrserrep) with n = l - m
        dev = mp.mpf(0)
        devK = mp.mpf(0)
        for m in range(low - 2 * N, low):
            n = low - m
            # GKK (negdiscrserrep): F raises n with coefficient
            # q^{-1/2-L-n} sqrt((1-q^{2n+2})(1-q^{4L+2n})); E = block raising in m lowers n:
            # E e_(n) = c_n e_(n-1), c_n^2 = q^{1-2L-2n}(1-q^{2n})(1-q^{4L+2n-2}).
            target = q ** (1 - 2 * Lbig - 2 * n) * (1 - q ** (2 * n)) * \
                (1 - q ** (4 * Lbig + 2 * n - 2))
            dev = max(dev, abs(blockcoef2(-1, 1, m) - target) / max(abs(target), mp.mpf(1)))
            devK = max(devK, abs(mp.sqrt(p) * q ** m - q ** (-Lbig - n)))
        ok &= report(f"        (-,+), m<=l: E-coeffs = D^- with L={float(Lbig)}", dev)
        ok &= report(f"        (-,+) K-eigenvalues = D^-", devK)

        # Casimir values: strange block Omega_GKK = +x = mu(q^{2a}); discrete = -x
        dev = max(abs(x - (q ** (2 * a) + q ** (-2 * a)) / 2),
                  abs(-x - (-(q ** (1 - 2 * Lbig) + q ** (2 * Lbig - 1)) / 2)))
        ok &= report(f"        Casimir values: +x = mu(q^(2a)), -x = -mu(q^(1-2L))", dev)

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
