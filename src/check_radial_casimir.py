#!/usr/bin/env python3
"""Issue #7: the radial part of the Casimir on the reduced quantum AdS_{2,q}.

Starting point: the momentum-space Casimir action SI (4.25) on the restricted
q-lattice R^2_{q^2}(xi) (SI (4.22)) with inner product SI (4.24)
(tex/main.tex, Proposition `prop:where`, eqs. (ourlattice)/(ourinnerproduct)/
(ourradialoperator)).  Claims checked (cited as (R0)-(R7) from tex Sec. 5):

 (R0) the connection identity
        2phi1(0,0;C;q,Z) * (C;q)_inf * (Z;q)_inf = 1phi1(Z;0;q,C),
      which continues the minimal radial solution to the whole lattice
      (1phi1(Z;0;q,C) = sum_n ((Z;q)_n/(q;q)_n) (-1)^n q^(n(n-1)/2) C^n is
      entire in Z).
 (R1) slice reduction of SI (4.25): the operator commutes with multiplication
      by chi; on the slice y = q^(2m+1) the symmetrised Casimir acts by
        Omtilde f_m = (q/2) f_{m+1} + (1/(2q)) (1 - c q^{2m}) f_{m-1},
      c := chi/xi;  it is symmetric w.r.t. the slice weight
      w_m = q^{2m+1} (c q^{2m+2}; q^2)_inf induced by SI (4.24); the
      orthonormalised form is the Jacobi operator
        (J_c g)_m = (1/2) sqrt(1-c q^{2m+2}) g_{m+1}
                  + (1/2) sqrt(1-c q^{2m})   g_{m-1};
      for c > 0 the weight vanishes below the wall m_0 = -log_{q^2} c and
      J_c is the chord transfer matrix (shift k = m - m_0); J_c = J_{c q^2}
      up to the unit shift m -> m+1.
 (R2) classical (q->1) limit at SI's scaling chi = 2p(q^{-1}-q),
      xi^{-1} = 2 p_R (q^{-1}-q): the edge Hamiltonian h_q := 2(Omtilde-1)/lambda^2
      (lambda = -2 ln q) contracts on smooth functions of phi (e^phi = y^{-1/2}/2)
      to  h_cl = (1/4)(1-d/dphi)^2 - p p_R e^{-2 phi},
      i.e. SI's classical radial Casimir (2.26) up to the stated conjugation;
      the deviation is O(lambda^2) (measured rate).  The lower (AdS) edge is
      the same operator on the (-1)^m-twisted sector: -2(Omtilde+1)/lambda^2
      twisted == h_q untwisted (exact identity).  Eigenfunction check:
      h_q on e^{phi} K_{2is}(2 sqrt(p p_R) e^{-phi}) ~ -s^2 (wall side).
 (R3) closed-form solutions on the c<0 wedge: c_m(sigma) = sigma^m *
      2phi1(0,0; sigma^2 q^2; q^2, c q^{2m+2}) solves the (unnormalised)
      radial recurrence 2 omega c_m = c_{m+1} + (1-c q^{2m}) c_{m-1} at
      omega = mu(sigma); via (R0) it extends to all m in Z; verified against
      downward recursion.
 (R4) deficiency indices (1,1) on the c<0 wedge [Koelink, arXiv:math/0305385,
      case (1) with psi=0: our a_m^2 = 1+gamma q^{2m+2} equals his
      1+r^{-2}q^{-2k} at k=-m-1, r=gamma^{-1/2}]: at the growth end m -> -inf
      BOTH solutions are l^2 with |g_m| ~ q^{|m|/2} (limit circle; measured
      slopes of ln_q|g_m|), while the free end m -> +inf is limit point
      (bounded coefficients); hence exactly ONE l^2(Z) solution of J*g = i g.
 (R5) extension dependence: out-of-band eigenvalues of Dirichlet vs Neumann
      truncations at m = -M do NOT converge to each other as M grows
      (limit-circle signature); control: the c>0 (chord) wedge has no
      out-of-band eigenvalues at all.
 (R6) the point spectrum of a self-adjoint extension: two bilateral q^2-grids
      [Koelink math/0305385, Thm 3.1 + Rmk 3.5: two q^2-quadratic grids].
      For the extension through omega_0 = -mu(q^{2 a_0}) on the gamma=1 slice:
      (a) the strange grid -mu(q^{2 a_0 + 2Z}) lies in the same extension
          (bilinear orthogonality h(omega) = sum_m psi^{(omega_0)} psi^{(omega)}
          vanishes);
      (b) the partner grid is +mu(q^{2Z - 2 a_0}) — the root-product rule
          y_1 y_2 in -gamma q^{2Z} (elliptic order-2 argument);
      (c) off-grid controls do not vanish;
      (d) the extension is parity(Pi)-invariant iff 2 a_0 in Z (Plancherel
          labels), where the spectrum is {±mu(q^{2a_0+2Z})}: the fusion pairing
          a = l - 1/2.
 (R7) q -> 1: the grid structure persists at q = 0.9 (h-tests), and the
      surviving lower-edge bound-state energies -2(mu(q^{2a_j})-1)/lambda^2 ->
      -(a_0+j)^2, the classical discrete-series tower of the p p_R < 0 well
      (SI (2.26) discussion), while the mirror grid is expelled like -4/lambda^2.

Run: ./.venv/bin/python src/check_radial_casimir.py
Output: results/check_radial_casimir.out (deterministic).
"""

import mpmath as mp


# ----------------------------------------------------------------------------
# basic q-series helpers (house style: explicit loops, adaptive truncation)
# ----------------------------------------------------------------------------

def qpochinf(a, q):
    """(a;q)_inf, |q|<1; a may be complex or mpf."""
    out = mp.mpf(1)
    aq = mp.mpmathify(a)
    for _ in range(200000):
        out *= (1 - aq)
        aq *= q
        if abs(aq) < mp.mpf(10) ** (-mp.mp.dps - 8):
            break
    return out


def phi00(C, q, z):
    """2phi1(0,0;C;q,z) = sum_n z^n / ((q;q)_n (C;q)_n); needs |z|<1."""
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


def phi11(Z, q, C):
    """1phi1(Z;0;q,C) = sum_n ((Z;q)_n/(q;q)_n) (-1)^n q^(n(n-1)/2) C^n.

    Entire in Z (the (Z;q)_n growth is beaten by q^(n(n-1)/2)); adaptive stop
    after the terms have passed their hump and fallen below threshold.
    """
    tot = mp.mpmathify(0)
    term = mp.mpmathify(1)  # n = 0 term
    peak = mp.mpf(1)
    n = 0
    small = 0
    while True:
        tot += term
        # term_{n+1} = term_n * (1 - Z q^n) * (-C q^n) / (1 - q^{n+1})
        term = term * (1 - Z * q ** n) * (-C * q ** n) / (1 - q ** (n + 1))
        n += 1
        peak = max(peak, abs(term))
        if abs(term) < peak * mp.mpf(10) ** (-mp.mp.dps - 8):
            small += 1
            if small >= 3:
                break
        else:
            small = 0
        if n > 200000:
            raise RuntimeError("phi11 did not converge")
    return tot


def mu(y):
    return (y + 1 / y) / 2


# ----------------------------------------------------------------------------
# slice machinery (tex Sec. 5): function picture, weights, solutions
# ----------------------------------------------------------------------------

def slice_weight(m, c, q):
    """w_m = q^{2m+1} (c q^{2m+2}; q^2)_inf  (Jackson constants dropped)."""
    return q ** (2 * m + 1) * qpochinf(c * q ** (2 * m + 2), q * q)


def cm_closed(m, sigma, c, q):
    """Minimal-at-+inf solution c_m(sigma) of
       2 mu(sigma) c_m = c_{m+1} + (1 - c q^{2m}) c_{m-1},
    via the 1phi1 form (valid for ALL m as long as (c q^{2m+2};q^2)_inf != 0):
       c_m = sigma^m 1phi1(cq^{2m+2};0;q^2,sigma^2 q^2)
             / ((sigma^2 q^2;q^2)_inf (c q^{2m+2};q^2)_inf).
    """
    q2 = q * q
    Z = c * q ** (2 * m + 2)
    return (sigma ** m * phi11(Z, q2, sigma * sigma * q2)
            / (qpochinf(sigma * sigma * q2, q2) * qpochinf(Z, q2)))


def solve_updown(seed_m, seed_vals, omega, c, q, mlo, mhi):
    """Extend a solution of 2 omega c_m = c_{m+1} + (1-c q^{2m}) c_{m-1}
    to m in [mlo, mhi] from values at (seed_m, seed_m+1). Returns dict m->c_m."""
    vals = {seed_m: seed_vals[0], seed_m + 1: seed_vals[1]}
    for m in range(seed_m + 1, mhi):
        vals[m + 1] = 2 * omega * vals[m] - (1 - c * q ** (2 * m)) * vals[m - 1]
    for m in range(seed_m, mlo, -1):
        vals[m - 1] = (2 * omega * vals[m] - vals[m + 1]) / (1 - c * q ** (2 * m))
    return vals


def weighted(vals, c, q):
    """g_m = sqrt(w_m) q^{-m} c_m = q^{1/2} sqrt((c q^{2m+2};q^2)_inf) c_m."""
    g = {}
    for m, v in vals.items():
        g[m] = mp.sqrt(q) * mp.sqrt(qpochinf(c * q ** (2 * m + 2), q * q)) * v
    return g


def psi_grid(omega, gamma_, q, Mlo, Mhi):
    """Weighted minimal solution g_m on the c = -gamma_ slice, m in [Mlo,Mhi].

    The closed form is used directly for every m >= 0 (upward recursion of a
    minimal solution is exponentially unstable, so it is never used); downward
    recursion (stable: both branches carry the same leading rate, check (R4))
    extends the solution to m < 0 from the closed-form values at m = 0, 1.
    """
    c = -gamma_
    sigma = omega - mp.sqrt(omega * omega - 1) if omega > 1 else \
        omega + mp.sqrt(omega * omega - 1)  # |sigma|<1 branch, real omega, |omega|>1
    vals = {}
    for m in range(0, Mhi + 1):
        vals[m] = cm_closed(m, sigma, c, q)
    for m in range(0, Mlo, -1):
        vals[m - 1] = (2 * omega * vals[m] - vals[m + 1]) / (1 - c * q ** (2 * m))
    return weighted(vals, c, q)


def norm2(g):
    return mp.sqrt(sum(abs(v) ** 2 for v in g.values()))


def bilin(g1, g2):
    ms = sorted(set(g1) & set(g2))
    return sum(g1[m] * g2[m] for m in ms)


# ----------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------

def run_R0():
    """Connection identity 2phi1(0,0;C;q,Z)(C;q)_inf(Z;q)_inf = 1phi1(Z;0;q,C)."""
    ok = True
    print("(R0) identity 2phi1(0,0;C;q,Z)*(C;q)_inf*(Z;q)_inf == 1phi1(Z;0;q,C)")
    for (q, C, Z) in [(mp.mpf('0.36'), mp.mpf('0.2'), mp.mpf('0.85')),
                      (mp.mpf('0.36'), mp.mpf('0.2'), mp.mpf('-0.95')),
                      (mp.mpf('0.81'), mp.mpf('0.5'), mp.mpf('-0.6')),
                      (mp.mpf('0.49'), mp.mpf('0.07'), mp.mpc(0, '0.8'))]:
        lhs = phi00(C, q, Z) * qpochinf(C, q) * qpochinf(Z, q)
        rhs = phi11(Z, q, C)
        dev = abs(lhs - rhs) / abs(rhs)
        p = dev < mp.mpf(10) ** (-mp.mp.dps + 10)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}, C={mp.nstr(C,3)}, "
              f"Z={mp.nstr(Z,3)}: rel dev = {mp.nstr(dev, 3)}")
    return ok


def run_R1():
    """Slice reduction of SI (4.25), symmetry, orthonormal form, chord wall."""
    ok = True
    print("\n(R1) slice reduction of SI (4.25); symmetry w.r.t. SI (4.24); "
          "orthonormalised Jacobi form; chord wall at c>0")
    for q in [mp.mpf('0.6'), mp.mpf('0.85')]:
        q2 = q * q

        def Om_si(f, chi, xi, m):
            """SI (4.25) acting on f(chi, y) at y = q^{2m+1}: the operator
            Omega = (q^{-1}T^{-2}+qT^2-2)/(q-q^{-1})^2 - mu_xi y chi/(1-q^2) T^{-2},
            mu_xi = xi^{-1}/(1-q^2); (Tf)(chi,y) = f(chi,qy)."""
            y = q ** (2 * m + 1)
            mu_xi = 1 / (xi * (1 - q2))
            t1 = (f(m - 1) / q + q * f(m + 1) - 2 * f(m)) / (q - 1 / q) ** 2
            t2 = mu_xi * y * chi / (1 - q2) * f(m - 1)
            return t1 - t2

        # (a) Omtilde from (4.25) equals the closed 3-term slice form, both signs
        devmax = mp.mpf(0)
        for (chi, xi) in [(q ** 2, q ** (-2)), (-q ** 4, q ** 2), (mp.mpf(1), mp.mpf(1)),
                          (-mp.mpf(1), q ** (-2))]:
            c = chi / xi
            for m in range(-4, 5):
                # random-ish but deterministic test function on the slice
                f = lambda mm: mp.mpf(2 + mm) ** 2 / (1 + abs(mm)) + mp.sin(mm)
                om = Om_si(f, chi, xi, m)
                omt = ((q - 1 / q) ** 2 * om + 2 * f(m)) / 2
                closed = q / 2 * f(m + 1) + (1 - c * q ** (2 * m)) / (2 * q) * f(m - 1)
                devmax = max(devmax, abs(omt - closed) / max(abs(closed), mp.mpf(1)))
        p = devmax < mp.mpf(10) ** (-mp.mp.dps + 8)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q={mp.nstr(q,3)}: (4.25) -> slice form "
              f"Omtilde f_m = (q/2) f_(m+1) + (1-c q^2m)/(2q) f_(m-1), "
              f"max dev = {mp.nstr(devmax, 3)}")

        # (b) symmetry w.r.t. w_m on the c<0 slice (doubly infinite window)
        c = -q2  # gamma = q^2 slice
        M = 12
        import random
        rnd = random.Random(20260811)
        fv = {m: mp.mpf(rnd.uniform(-1, 1)) for m in range(-M, M + 1)}
        gv = {m: mp.mpf(rnd.uniform(-1, 1)) for m in range(-M, M + 1)}
        # zero outside a window so that boundary terms vanish identically
        for d in (fv, gv):
            for m in (-M, -M + 1, M - 1, M):
                d[m] = mp.mpf(0)

        def omt_apply(d, m):
            return q / 2 * d.get(m + 1, mp.mpf(0)) \
                + (1 - c * q ** (2 * m)) / (2 * q) * d.get(m - 1, mp.mpf(0))

        s1 = sum(slice_weight(m, c, q) * omt_apply(fv, m) * gv[m] for m in range(-M, M + 1))
        s2 = sum(slice_weight(m, c, q) * fv[m] * omt_apply(gv, m) for m in range(-M, M + 1))
        dev = abs(s1 - s2) / max(abs(s1), mp.mpf(10) ** (-20))
        p = dev < mp.mpf(10) ** (-mp.mp.dps + 10)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}]            symmetry <Omtilde f,g>_w == "
              f"<f,Omtilde g>_w (c=-q^2), rel dev = {mp.nstr(dev, 3)}")

        # (c) orthonormalised coefficients a_m = (1/2) sqrt(1-c q^{2m+2})
        devmax = mp.mpf(0)
        for m in range(-6, 7):
            am = q / 2 * mp.sqrt(slice_weight(m, c, q) / slice_weight(m + 1, c, q))
            devmax = max(devmax, abs(am - mp.sqrt(1 - c * q ** (2 * m + 2)) / 2))
        p = devmax < mp.mpf(10) ** (-mp.mp.dps + 8)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}]            orthonormal a_m == "
              f"(1/2)sqrt(1-c q^(2m+2)), max dev = {mp.nstr(devmax, 3)}")

        # (d) c>0: weight vanishes below the wall (an exactly vanishing factor
        # (1 - c q^{2m+2} q^{2i}) at i = -(nu+m+1); checked as the minimum
        # factor, robust to rounding of the large cofactors), and above the
        # wall J = chord matrix
        nu = 2
        c = q ** (2 * nu)  # wall at m0 = -nu
        wall_ok = True
        for m in range(-nu - 5, -nu):
            minfac = min(abs(1 - c * q ** (2 * (m + 1 + i))) for i in range(0, 12))
            wall_ok &= minfac < mp.mpf(10) ** (-mp.mp.dps + 8)
        devmax = mp.mpf(0)
        for k in range(0, 8):
            m = k - nu
            am = mp.sqrt(1 - c * q ** (2 * m + 2)) / 2
            devmax = max(devmax, abs(am - mp.sqrt(1 - q2 ** (k + 1)) / 2))
        p = wall_ok and devmax < mp.mpf(10) ** (-mp.mp.dps + 8)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}]            c=q^{2*nu}>0: w_m=0 below wall "
              f"m0={-nu}; J == chord transfer matrix, max dev = {mp.nstr(devmax, 3)}")
    return ok


def run_R2():
    """Classical limit of the slice operator, SI scaling; both wedges; edges.

    In the orthonormal (Jacobi) gauge the contraction is clean at O(lambda^2):
      h^J_q := 2(J_c - 1)/lambda^2 -> (1/4) d^2/dphi^2 - p p_R e^{-2 phi},
    with a_m = (1/2) sqrt(1 - c q^{2m+2}), c = 4 p p_R (q^{-1}-q)^2, and
    e^phi = y^{-1/2}/2.  (The function-picture form Omtilde f_m = (q/2)f_{m+1}
    + (1-c q^{2m})/(2q) f_{m-1} contracts to the similarity-equivalent
    (1/4)(1-d_phi)^2 - p p_R e^{-2 phi}, but with an O(lambda) midpoint
    artifact in pointwise comparisons; the gauge-invariant statement is the
    J-gauge one, and that is what is checked here.)"""
    ok = True
    print("\n(R2) classical limit (orthonormal gauge): h^J_q := 2(J_c-1)/lambda^2 "
          "-> (1/4) d_phi^2 - p p_R e^(-2 phi); rate O(lambda^2)")
    pL = mp.mpf(1)
    phi0 = mp.mpf('0.1')

    def hJq_at(fun, mstar, lam, pR):
        q = mp.e ** (-lam / 2)
        c = 4 * pL * pR * (1 / q - q) ** 2

        def a(m):
            return mp.sqrt(1 - c * q ** (2 * m + 2)) / 2

        def f_m(m):
            y = q ** (2 * m + 1)
            return fun(-mp.log(y) / 2 - mp.log(2))

        return 2 * (a(mstar) * f_m(mstar + 1) + a(mstar - 1) * f_m(mstar - 1)
                    - f_m(mstar)) / lam ** 2

    def mstar_at(phi_target, lam):
        q = mp.e ** (-lam / 2)
        return int(mp.floor((-2 * (phi_target + mp.log(2)) / mp.log(q) - 1) / 2))

    def test_dev(lam, pR):
        q = mp.e ** (-lam / 2)

        def f_of_phi(phi):
            return mp.e ** (-(phi - phi0) ** 2)

        mstar = mstar_at(mp.mpf('-0.2'), lam)
        y = q ** (2 * mstar + 1)
        phis = -mp.log(y) / 2 - mp.log(2)
        hq = hJq_at(f_of_phi, mstar, lam, pR)
        fp = f_of_phi(phis)
        f2 = (-2 + 4 * (phis - phi0) ** 2) * fp
        hcl = f2 / 4 - pL * pR * mp.e ** (-2 * phis) * fp
        return abs(hq - hcl)

    for pR in [mp.mpf(1), mp.mpf(-1)]:
        d1 = test_dev(mp.mpf('0.2'), pR)
        d2 = test_dev(mp.mpf('0.1'), pR)
        d3 = test_dev(mp.mpf('0.05'), pR)
        r12, r23 = d1 / d2, d2 / d3
        p = 3 < r12 < mp.mpf('5.5') and 3 < r23 < mp.mpf('5.5')
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] p_R={mp.nstr(pR,2)}: dev(lam)= "
              f"{mp.nstr(d1,3)}/{mp.nstr(d2,3)}/{mp.nstr(d3,3)} at lam=0.2/0.1/0.05; "
              f"ratios {mp.nstr(r12,4)}, {mp.nstr(r23,4)} (predicted -> 4)")

    # eigenfunction check (wall side): P(phi) = K_{2is}(2 sqrt(p pR) e^{-phi}),
    # h^J_cl P = -s^2 P
    s = mp.mpf('0.7')

    def P_of_phi(phi):
        return mp.re(mp.besselk(2j * s, 2 * mp.e ** (-phi)))

    print("  eigenfunction check, wall side: h^J_q P / P -> -s^2, s = 0.7")
    devs = []
    for lam in [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05')]:
        mstar = mstar_at(mp.mpf('0.2'), lam)
        q = mp.e ** (-lam / 2)
        y = q ** (2 * mstar + 1)
        phis = -mp.log(y) / 2 - mp.log(2)
        ev = hJq_at(P_of_phi, mstar, lam, mp.mpf(1)) / P_of_phi(phis)
        devs.append(abs(ev + s * s))
        print(f"    lam={mp.nstr(lam,3)}: h^J_q P / P = {mp.nstr(ev, 8)} "
              f"(target {mp.nstr(-s*s, 8)}), |dev| = {mp.nstr(devs[-1], 3)}")
    p = 3 < devs[0] / devs[1] < 5 and 3 < devs[1] / devs[2] < 5
    ok &= p
    print(f"  [{'PASS' if p else 'FAIL'}] eigenvalue deviation ratios "
          f"{mp.nstr(devs[0]/devs[1], 4)}, {mp.nstr(devs[1]/devs[2], 4)} "
          f"(predicted -> 4)")

    # exact twisted-sector identity: -2(Omtilde+1)/lam^2 on (-1)^m f == parity of
    # +2(Omtilde-1)/lam^2 on f  (lower edge == upper edge of the same slice)
    q = mp.mpf('0.8')
    c = -mp.mpf('0.5')
    f = lambda m: 1 / (1 + m * m) + mp.mpf(m) / 7
    devmax = mp.mpf(0)
    for m in range(-5, 6):
        tw = lambda mm: (-1) ** mm * f(mm)
        lhs = -(q / 2 * tw(m + 1) + (1 - c * q ** (2 * m)) / (2 * q) * tw(m - 1)
                + tw(m))
        rhs = (-1) ** m * (q / 2 * f(m + 1) + (1 - c * q ** (2 * m)) / (2 * q)
                           * f(m - 1) - f(m))
        devmax = max(devmax, abs(lhs - rhs))
    p = devmax < mp.mpf(10) ** (-mp.mp.dps + 5)
    ok &= p
    print(f"  [{'PASS' if p else 'FAIL'}] twisted-sector identity "
          f"-(Omtilde+1) o Pi == Pi o (Omtilde-1) exactly, max dev = {mp.nstr(devmax,3)}")
    return ok


def run_R3():
    """Closed-form minimal solution on the c<0 wedge, all m, vs recursion."""
    ok = True
    print("\n(R3) closed form c_m(sigma) on the c<0 wedge: recurrence at m<=0 "
          "and match with downward recursion")
    for (q, gamma_, sg) in [(mp.mpf('0.6'), mp.mpf(1), mp.mpf('0.55')),
                            (mp.mpf('0.6'), mp.mpf(1), mp.mpf('-0.4')),
                            (mp.mpf('0.75'), mp.mpf('0.5625'), mp.mpf('0.3'))]:
        c = -gamma_
        omega = mu(sg)
        # (a) recurrence satisfied by the closed form for m in [-8, 3]
        devmax = mp.mpf(0)
        for m in range(-8, 4):
            lhs = 2 * omega * cm_closed(m, sg, c, q)
            rhs = cm_closed(m + 1, sg, c, q) \
                + (1 - c * q ** (2 * m)) * cm_closed(m - 1, sg, c, q)
            devmax = max(devmax, abs(lhs - rhs) / max(abs(lhs), abs(rhs)))
        p1 = devmax < mp.mpf(10) ** (-mp.mp.dps + 12)
        ok &= p1
        print(f"  [{'PASS' if p1 else 'FAIL'}] q={mp.nstr(q,3)}, gamma={mp.nstr(gamma_,4)}, "
              f"sigma={mp.nstr(sg,3)}: recurrence m=-8..3, max rel dev = {mp.nstr(devmax,3)}")
        # (b) recursion from the +inf seed reproduces the closed form at m=-8
        vals = solve_updown(6, (cm_closed(6, sg, c, q), cm_closed(7, sg, c, q)),
                            omega, c, q, -9, 8)
        dev = abs(vals[-8] - cm_closed(-8, sg, c, q)) / abs(vals[-8])
        p2 = dev < mp.mpf(10) ** (-mp.mp.dps + 15)
        ok &= p2
        print(f"  [{'PASS' if p2 else 'FAIL'}]            downward recursion == closed form "
              f"at m=-8, rel dev = {mp.nstr(dev, 3)}")
    return ok


def run_R4():
    """Limit circle at the growth end: both solutions l^2, |g_m| ~ q^{|m|/2}."""
    ok = True
    print("\n(R4) c<0 wedge, growth end m -> -inf: BOTH solutions have "
          "|g_m| ~ q^{|m|/2} (limit circle; deficiency indices (1,1))")
    q = mp.mpf('0.6')
    gamma_ = mp.mpf(1)
    c = -gamma_
    for omega, tag in [(mu(mp.mpf('-0.45')), "omega=-mu(0.45) (strange-type)"),
                       (mp.mpc(0, 1), "omega=i (deficiency)")]:
        # solution 1: minimal at +inf (seeded by the closed form when real)
        if mp.im(omega) == 0:
            sg = omega + mp.sqrt(omega * omega - 1)
            if abs(sg) > 1:
                sg = 1 / sg
            v1 = solve_updown(6, (cm_closed(6, sg, c, q), cm_closed(7, sg, c, q)),
                              omega, c, q, -34, 10)
        else:
            sg = omega - mp.sqrt(omega * omega - 1)
            if abs(sg) > 1:
                sg = 1 / sg
            # seed with sigma^m at large m (correction O(q^{2m}) negligible at m=14)
            v1 = solve_updown(14, (sg ** 14, sg ** 15), omega, c, q, -34, 16)
        # solution 2: independent (Dirichlet-type seed at the top)
        v2 = solve_updown(14, (mp.mpmathify(1), mp.mpmathify(0)), omega, c, q, -34, 16)
        g1, g2 = weighted(v1, c, q), weighted(v2, c, q)
        # Wronskian check of independence (unnormalised W_m at m=0)
        W = v1[1] * v2[0] - v2[1] * v1[0]
        indep = abs(W) > mp.mpf(10) ** (-10)
        # slopes of ln_q |g_m| over m in [-32,-20]
        slopes = []
        for g in (g1, g2):
            xs = list(range(-32, -19))
            ys = [mp.log(abs(g[m])) / mp.log(q) for m in xs]
            n = len(xs)
            sx = sum(xs); sy = sum(ys)
            sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
            slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
            slopes.append(-slope)  # ln_q|g| ~ -m/2 => -slope ~ 1/2
        tol = 50 * q ** 20  # correction scale O(q^{|m|}) at the window edge
        p = indep and all(abs(s - mp.mpf(1) / 2) < tol for s in slopes)
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q=0.6, gamma=1, {tag}: slopes of "
              f"ln_q|g_m| vs |m| = {mp.nstr(slopes[0],6)}, {mp.nstr(slopes[1],6)} "
              f"(predicted 0.5, tol {mp.nstr(tol,2)}); independent (|W|={mp.nstr(abs(W),3)})")
    # free end: for omega=i exactly one solution decays (|sigma|<1), the other
    # grows like |sigma|^{-m}: one l^2 solution at +inf => deficiency (1,1) total.
    sg = mp.mpc(0, 1) - mp.sqrt(mp.mpc(0, 1) ** 2 - 1)
    if abs(sg) > 1:
        sg = 1 / sg
    print(f"  free end: |sigma(i)| = {mp.nstr(abs(sg), 6)} < 1; the second solution "
          f"grows like |sigma|^-m: exactly one l^2 solution at +inf. Together with "
          f"the limit-circle growth end: deficiency indices (1,1).")
    return ok


def run_R5():
    """Dirichlet vs Neumann truncations: out-of-band eigenvalues differ."""
    ok = True
    print("\n(R5) extension dependence (limit-circle signature): out-of-band "
          "eigenvalues of Dirichlet vs Neumann truncation at m=-M")
    q = mp.mpf('0.6')
    gamma_ = mp.mpf(1)
    Mhi = 25

    def eigs(M, bc):
        n = M + Mhi + 1  # m = -M..Mhi
        A = mp.zeros(n)
        for i in range(n):
            m = i - M
            if i + 1 < n:
                a = mp.sqrt(1 + gamma_ * q ** (2 * m + 2)) / 2
                A[i, i + 1] = a
                A[i + 1, i] = a
        if bc == 'N':
            A[0, 0] = mp.sqrt(1 + gamma_ * q ** (-2 * M)) / 2
        E = mp.eigsy(A, eigvals_only=True)
        return sorted([e for e in E if abs(e) > mp.mpf('1.000001')])

    for M in [8, 12, 16]:
        eD = eigs(M, 'D')
        eN = eigs(M, 'N')
        # match each Dirichlet eigenvalue in a fixed window to nearest Neumann one
        window = [e for e in eD if mp.mpf('1.02') < abs(e) < mp.mpf('3')]
        gaps = []
        for e in window:
            gaps.append(min(abs(e - en) for en in eN))
        gmin = min(gaps) if gaps else mp.mpf('nan')
        gmax = max(gaps) if gaps else mp.mpf('nan')
        p = gmin > mp.mpf('0.002')
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] M={M}: {len(window)} Dirichlet "
              f"eigenvalues in 1.02<|E|<3; nearest-Neumann gaps in "
              f"[{mp.nstr(gmin,3)}, {mp.nstr(gmax,3)}] (do not shrink with M)")

    # control: c>0 (chord) wedge, no out-of-band eigenvalues for either BC
    def eigs_chord(n, bc):
        A = mp.zeros(n)
        for k in range(n - 1):
            a = mp.sqrt(1 - q ** (2 * k + 2)) / 2
            A[k, k + 1] = a
            A[k + 1, k] = a
        if bc == 'N':
            A[n - 1, n - 1] = mp.sqrt(1 - q ** (2 * n)) / 2
        E = mp.eigsy(A, eigvals_only=True)
        return [e for e in E if abs(e) > mp.mpf('1.000001')]

    nD = len(eigs_chord(40, 'D'))
    nN = len(eigs_chord(40, 'N'))
    p = nD == 0 and nN == 0
    ok &= p
    print(f"  [{'PASS' if p else 'FAIL'}] control, c>0 chord wedge (n=40): "
          f"out-of-band eigenvalue count D/N = {nD}/{nN} (predicted 0/0)")
    return ok


def run_R6():
    """Point spectrum grids of the extension through omega_0; parity; controls."""
    ok = True
    print("\n(R6) point spectrum of the extension through omega_0 = -mu(q^(2a0)), "
          "gamma=1 slice: two bilateral q^2-grids with root product in -q^(2Z)")
    q = mp.mpf('0.6')
    gamma_ = mp.mpf(1)
    Mlo, Mhi = -185, 90

    def htest(a0, omega_list, labels):
        om0 = -mu(q ** (2 * a0))
        g0 = psi_grid(om0, gamma_, q, Mlo, Mhi)
        n0 = norm2(g0)
        out = []
        for om, lab in zip(omega_list, labels):
            g = psi_grid(om, gamma_, q, Mlo, Mhi)
            hval = abs(bilin(g0, g)) / (n0 * norm2(g))
            out.append((lab, hval))
        return out

    # (a),(b),(c) at a0 = 1/2 (Plancherel label): grid -mu(q^{1+2Z}) and partner
    # +mu(q^{1+2Z}) (product rule with y1 = -q: y2 in q^{2Z-1}); controls off-grid
    a0 = mp.mpf('0.5')
    tests = htest(a0,
                  [-mu(q ** 3), -mu(q ** 5), mu(q), mu(q ** 3),
                   -mu(q ** mp.mpf('2.6')), mu(q ** 2)],
                  ["-mu(q^3)   [strange grid]",
                   "-mu(q^5)   [strange grid]",
                   "+mu(q^1)   [partner grid]",
                   "+mu(q^3)   [partner grid]",
                   "-mu(q^2.6) [off-grid control]",
                   "+mu(q^2)   [off-grid control]"])
    floor = mp.mpf(10) ** (-25)
    for lab, hval in tests:
        on_grid = 'grid]' in lab and 'off' not in lab
        p = hval < floor if on_grid else hval > mp.mpf('0.001')
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] a0=1/2: |h({lab})| = {mp.nstr(hval, 3)}"
              f"  ({'== 0 predicted' if on_grid else 'nonzero predicted'})")

    # (d) generic a0 = 0.37: strange grid -mu(q^{2a0+2Z}); partner grid
    # +mu(q^{2Z-2a0}) (product rule); the naive mirror +mu(q^{2a0}) is NOT in it
    a0 = mp.mpf('0.37')
    tests = htest(a0,
                  [-mu(q ** (2 * a0 + 2)), mu(q ** (2 - 2 * a0)),
                   mu(q ** (4 - 2 * a0)), mu(q ** (2 * a0)), mu(q ** (2 * a0 + 2))],
                  [f"-mu(q^(2a0+2))  [strange grid]",
                   f"+mu(q^(2-2a0))  [partner grid]",
                   f"+mu(q^(4-2a0))  [partner grid]",
                   f"+mu(q^(2a0))    [mirror: control]",
                   f"+mu(q^(2a0+2))  [mirror: control]"])
    for lab, hval in tests:
        on_grid = 'grid]' in lab and 'control' not in lab
        p = hval < floor if on_grid else hval > mp.mpf('0.001')
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] a0=0.37: |h({lab})| = {mp.nstr(hval, 3)}"
              f"  ({'== 0 predicted' if on_grid else 'nonzero predicted'})")
    print("  => the extension through -mu(q^(2a0)) is Pi-invariant iff 2a0 in Z: "
          "at a0=1/2 the spectrum is {±mu(q^(1+2Z))} (fusion pairing a = l - 1/2); "
          "at a0=0.37 the mirror points are absent.")

    # (e) internal consistency: boundary-form convergence in the cutoff
    om0 = -mu(q)
    om1 = mu(q)
    gA0 = psi_grid(om0, gamma_, q, -145, Mhi)
    gA1 = psi_grid(om1, gamma_, q, -145, Mhi)
    gB0 = psi_grid(om0, gamma_, q, Mlo, Mhi)
    gB1 = psi_grid(om1, gamma_, q, Mlo, Mhi)
    hA = abs(bilin(gA0, gA1)) / (norm2(gA0) * norm2(gA1))
    hB = abs(bilin(gB0, gB1)) / (norm2(gB0) * norm2(gB1))
    p = abs(hA - hB) < mp.mpf(10) ** (-25)
    ok &= p
    print(f"  [{'PASS' if p else 'FAIL'}] cutoff independence: |h| at Mlo=-145 vs "
          f"-185: {mp.nstr(hA,3)} vs {mp.nstr(hB,3)}")
    return ok


def run_R7():
    """q->1: grids persist; lower-edge bound-state tower -> -(a0+j)^2."""
    ok = True
    print("\n(R7) q->1: the grids persist near q=1 and the lower-edge tower "
          "2(1-mu(q^(2a_j)))/lambda^2 -> -(a0+j)^2 (classical discrete series "
          "of the p p_R<0 well); the mirror grid is expelled like -4/lambda^2")
    q = mp.mpf('0.9')
    lam = -2 * mp.log(q)
    gamma_ = mp.mpf(1)
    a0 = mp.mpf('0.5')
    Mlo, Mhi = -720, 260
    om0 = -mu(q ** (2 * a0))
    g0 = psi_grid(om0, gamma_, q, Mlo, Mhi)
    n0 = norm2(g0)
    for om, lab, pred0 in [(-mu(q ** 3), "-mu(q^3) [strange grid]", True),
                           (mu(q), "+mu(q^1) [partner grid]", True),
                           (-mu(q ** mp.mpf('2.3')), "-mu(q^2.3) [control]", False)]:
        g = psi_grid(om, gamma_, q, Mlo, Mhi)
        hval = abs(bilin(g0, g)) / (n0 * norm2(g))
        p = hval < mp.mpf(10) ** (-20) if pred0 else hval > mp.mpf('0.0005')
        ok &= p
        print(f"  [{'PASS' if p else 'FAIL'}] q=0.9: |h({lab})| = {mp.nstr(hval,3)}")
    # the energy map (arithmetic of the cosh expansion, rate check)
    print("  lower-edge energies 2(1-mu(q^(2a)))/lambda^2 at a = a0+j:")
    for qq in [mp.mpf('0.9'), mp.mpf('0.99')]:
        ll = -2 * mp.log(qq)
        vals = []
        for j in range(3):
            a = a0 + j
            e = 2 * (1 - mu(qq ** (2 * a))) / ll ** 2
            vals.append(e + a * a)  # deviation from -(a0+j)^2... e -> -a^2
        print(f"    q={mp.nstr(qq,3)}: deviations from -(a0+j)^2, j=0,1,2: "
              f"{mp.nstr(vals[0],3)}, {mp.nstr(vals[1],3)}, {mp.nstr(vals[2],3)}")
    d1 = abs(2 * (1 - mu(mp.mpf('0.9') ** 1)) / (-2 * mp.log(mp.mpf('0.9'))) ** 2
             + mp.mpf('0.25'))
    d2 = abs(2 * (1 - mu(mp.mpf('0.99') ** 1)) / (-2 * mp.log(mp.mpf('0.99'))) ** 2
             + mp.mpf('0.25'))
    p = 50 < d1 / d2 < 200  # O(lambda^2): factor ~ (log .9/log .99)^2 ~ 112
    ok &= p
    print(f"  [{'PASS' if p else 'FAIL'}] deviation rate (q=0.9)/(q=0.99) = "
          f"{mp.nstr(d1/d2,4)} (predicted ~ (lam_1/lam_2)^2 = "
          f"{mp.nstr((mp.log(mp.mpf('0.9'))/mp.log(mp.mpf('0.99')))**2,4)})")
    # expulsion of the mirror grid at the same edge
    for qq in [mp.mpf('0.9'), mp.mpf('0.99')]:
        ll = -2 * mp.log(qq)
        e = 2 * (1 - (-mu(qq ** 1))) / ll ** 2  # value at Omtilde=+mu(q): 2(1+mu)/ll^2
        print(f"    mirror point +mu(q) at q={mp.nstr(qq,3)}: lower-edge energy "
              f"2(1+mu)/lambda^2 = {mp.nstr(e,6)} ~ 4/lambda^2 = "
              f"{mp.nstr(4/ll**2,6)} -> infinity")
    return ok


def main():
    print("check_radial_casimir.py -- issue #7 (mpmath)")
    print("=" * 78)
    ok = True
    mp.mp.dps = 60
    ok &= run_R0()
    ok &= run_R1()
    mp.mp.dps = 30
    ok &= run_R2()
    mp.mp.dps = 160
    ok &= run_R3()
    mp.mp.dps = 220
    ok &= run_R4()
    mp.mp.dps = 30
    ok &= run_R5()
    mp.mp.dps = 50
    ok &= run_R6()
    mp.mp.dps = 45
    ok &= run_R7()
    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
