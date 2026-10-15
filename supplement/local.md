# Local lemmas near s = 0 and s = 1/2

**Evidence status and paths.** The mathematical derivations below form this supplement. Numerical tables, timings and validation outcomes are historical development summaries, not a fresh independent production replay. Their rounded decimal margins are diagnostics, not outward-rounded certificates. Proof-critical scripts are in `code/local/` and outputs in `certificates/local/`; current hashes, exact exported bounds and acceptance results belong to those files and the release replay receipt. Auxiliary samples, formula checks, negative controls and wider-region runs do not replace exhaustive production enclosures. No unshipped development log is a dependency.

Throughout, E, E', R, T are independent standard exponential variables,
s ∈ [0, 1/2], c = 1 − s, d = 1 − 2s, N = N(s) = s² + c², and

  Z_s = s(E − 1) − c(E' − 1) = d + sE − cE',  W = ∂Z_s/∂s = E + E' − 2,
  M(s,p) = E|Z_s|^p,  G(s,p) = (1/p) log M(s,p) − (1/2) log N(s).

For p > 1 the derivatives below are taken in s at fixed p. On the sets where it is used,
x = d/c ∈ [0,1] and y = d/s ∈ [0, ∞]. Two identities recur: 2c − x = N/c and y + 2s = N/s.

All formulas below were also checked numerically (Section 6) before anything was certified.

---

## 1. Conditioning representation

**Lemma 1 (law of (Z, W) given the side).** First let 0 < s ≤ 1/2 and X = sE − cE', so Z_s = d + X. Then P(X > 0) = s and
P(X < 0) = c. Moreover:

- conditional on {X > 0}, (Z, W) has the law of (d + sR, (R − 1) + (T − 1));
- conditional on {X < 0}, (Z, W) has the law of (d − cR, (R − 1) + (T − 1)),

where R and T are independent standard exponentials.

*Proof.* On {sE > cE'}, condition on E' = e. By memorylessness, R := E − ce/s is Exp(1) and independent
of E'. The conditional density of E' on the event is proportional to e^{−e} e^{−ce/s} = e^{−e/s}, so
E' = sT with T ~ Exp(1) independent of R, and P(sE > cE') = s/(s + c) = s. Then E = R + cT,
X = sE − cE' = sR, and W = R + cT + sT − 2 = (R − 1) + (T − 1). The case cE' > sE is the same with
the roles of E and E' exchanged: E = cT, E' = R + sT, X = −cR and W = (R − 1) + (T − 1). At s = 0, use X = −E', R = E', T = E directly; no conditioning on the null positive-side event is required. The weighted identities extend to this endpoint. ∎

**Corollary 1.** M(s,p) = s E[(d + sR)^p] + c E[|d − cR|^p]. Given the side and R,
E[W | ·] = R − 1 and E[W² | ·] = (R − 1)² + 1.

## 2. Derivative formulas

**Lemma 2 (derivatives of M).** Let p > 1, and define Z_s by the same formula for every real s. Then
M(·,p) is C¹ on [−1, 1] and C² on [0, 1]; on [1/2, 1] this follows from [0, 1/2] by the symmetry
M(1−s,p) = M(s,p) of Proposition A(2). Its derivatives are

  M_s = p E[|Z|^{p−2} Z W] = p ( s E[(d+sR)^{p−1}(R−1)] + c E[|d−cR|^{p−2}(d−cR)(R−1)] ),
  M_ss = p(p−1) E[|Z|^{p−2} W²] = p(p−1) ( s E[(d+sR)^{p−2}((R−1)²+1)] + c E[|d−cR|^{p−2}((R−1)²+1)] ).

*Proof.* Z_s is affine in s with slope W. For |s| ≤ 1 we have |Z_s| ≤ |E − 1| + 2|E' − 1| ≤ 2(E + E' + 2)
and |W| ≤ E + E' + 2. Hence ∂_s|Z_s|^p = p|Z_s|^{p−2}Z_s W is dominated by p 2^{p−1}(E+E'+2)^p, which is
integrable, and dominated convergence gives the formula for M_s.

For M_ss, put φ(z) = |z|^{p−2}z. It is absolutely continuous on ℝ with φ' = (p−1)|z|^{p−2}, and this
derivative is locally integrable because p > 1. So for each ω, u ↦ φ(Z_u)W is absolutely continuous
with derivative (p−1)|Z_u|^{p−2}W² almost everywhere. Integrating over [s₁, s₂] and using Tonelli for the
nonnegative integrand gives

  M_s(s₂) − M_s(s₁) = ∫_{s₁}^{s₂} p(p−1) E[|Z_u|^{p−2}W²] du.

The right-hand side is finite. Indeed, by Lemma 1 the integrand equals the conditional expression in the
statement, and Lemma 4 (Section 4) gives that expression in closed form, continuous in u ∈ [0, 1/2]. So M_s is C¹ with
derivative M_ss. The conditional forms follow from Lemma 1, because on each side Z is a function of R
and T is independent of R. ∎

**Lemma 3 (normalised form of G_s and G_ss).** Let p > 1, τ = N'/(2N) = −d/N and V = W − τZ. Put
A = E[|Z|^{p−2}ZV]/M and B = E[|Z|^{p−2}V²]/M. Then

  G_s = A,  G_ss = (p − 1) B − p A² + 2 (d/N) A − 1/N².

Given the side and R, E[V | ·] = c(R − 2s)/N on the + side and s(R − 2c)/N on the − side. Hence

  A = (sc/N) · a/M,  a = E[(d+sR)^{p−1}(R − 2s)] + E[|d−cR|^{p−2}(d−cR)(R − 2c)],
  B = b/M,  b = s E[(d+sR)^{p−2}(c²(R−2s)²/N² + 1)] + c E[|d−cR|^{p−2}(s²(R−2c)²/N² + 1)].

*Proof.* G_s = M_s/(pM) − τ, and M_s = pE[|Z|^{p−2}Z(V + τZ)] = pM(A + τ), so G_s = A. Next,

  M_ss/M = p(p−1)(B + 2τA + τ²),  (M_s/M)² = p²(A + τ)²,
  G_ss = (1/p)(M_ss/M − (M_s/M)²) − τ' = (p−1)B − 2τA − pA² − (τ² + τ').

Now τ' = (N''N − N'²)/(2N²) = 4s(1−s)/N² and τ² = d²/N², so τ² + τ' = (d² + 4sc)/N² = 1/N². This
gives the formula for G_ss.

For the conditional means, on the + side E[V | ·] = (R − 1) − τ(d + sR) = R(1 − τs) − (1 + τd). Since
1 − τs = (N + sd)/N = c/N and 1 + τd = (N − d²)/N = 2sc/N, this equals c(R − 2s)/N. On the − side,
1 + τc = s/N, which gives s(R − 2c)/N. Finally E[V² | ·] = E[V | ·]² + Var(T) = E[V | ·]² + 1. ∎

Two remarks on Lemma 3.

- The algebra was checked numerically against the direct derivative combination
  (M_ss M − M_s²)/(pM²) − (N''N − N'²)/(2N²), with agreement to 10^{−49} (Section 6, C3).
- This normalised form is what makes large p cheap. The O(p²) terms of that combination cancel
  analytically, and every term is O(1). At s = 0, B = F(p−2)/F(p) ~ 1/p².

## 3. Part (a): the first derivative vanishes at s = 0 and at s = 1/2

**Proposition A.** For every p > 1:

1. G_s(0, p) = 0;
2. G_s(1/2, p) = 0;
3. G_ss(0, p) = ((p−1)/p)(1/F(p) − 1), where F(p) = E|E − 1|^p;
4. G_ss(1/2, p) = 4(p − 2)²/p.

*Proof.*

(1) **s = 0.** Here c = d = 1 and the + side has weight s = 0. Lemma 2 gives
M_s(0,p) = p E[|1 − R|^{p−2}(1 − R)(R − 1)] = −p E|1 − R|^p = −pF(p), while M(0,p) = F(p) and
N(0) = 1, N'(0) = −2. Hence G_s(0,p) = M_s/(pM) − N'/(2N) = −1 + 1 = 0. Equivalently, A carries the
factor s in Lemma 3. The derivative is two-sided: M is C¹ on [−1, 1].

(2) **s = 1/2.** Exchanging E and E' gives Z_{1−s} =_d −Z_s, so M(1−s,p) = M(s,p). Also
N(1−s) = N(s). Hence G(1−s,p) = G(s,p) for all real s, and G(·,p) is C¹ on [−1, 1] by Lemma 2. Differentiating
at s = 1/2 gives G_s(1/2) = −G_s(1/2), so G_s(1/2,p) = 0. As a direct check: at s = 1/2, d = x = 0 and
the inner integral is empty, so a = 2^{1−p}(Γ(p+1) − Γ(p)) − 2^{1−p}(Γ(p+1) − Γ(p)) = 0.

(3) **G_ss at s = 0.** Here A = 0, N = 1, and V = W + Z = E − 1 is independent of Z = 1 − E'. So
B = F(p−2)/F(p) and G_ss(0,p) = (p−1)F(p−2)/F(p) − 1.

Let f(t) = |t − 1|^p. Integrating by parts twice against e^{−t} (f' is absolutely continuous for p > 1)
gives E f''(E) = E f(E) − f(0) − f'(0). That is, p(p−1)F(p−2) = F(p) − 1 + p. Substituting,

  G_ss(0,p) = 1/p + (p−1)/(pF) − 1 = ((p−1)/p)(1/F − 1).

(4) **G_ss at s = 1/2.** Here A = 0, d = 0, N = 1/2 and M = Γ(p+1)2^{−p}. Both conditional means
equal R − 1 and |Z| = R/2, so

  (p−1)B = 4(p−1)(Γ(p+1) − 2Γ(p) + 2Γ(p−1))/Γ(p+1) = 4(p − 3 + 4/p).

Then G_ss = 4(p − 3 + 4/p) − 4 = 4(p−2)²/p. ∎

**Consequences used in the main theorem.** For p > 1, G(·,p) is C² on [0,1/2] (Lemma 2, and M > 0).
By Taylor's theorem with integral remainder,
G_s(s,p) = G_s(s₀,p) + ∫_{s₀}^{s} G_ss, so:

- **L2** (G_ss > 0 on [0, 1/20]) with G_s(0) = 0 gives G_s > 0 on (0, 1/20], hence G(s,p) ≥ G0(p) there.
- **L3** (G_ss < 0) gives G(s,p) ≤ G0(p) on [0, 1/20].
- **L1** (G_ss > 0 on [9/20, 1/2]) with G_s(1/2) = 0 gives G_s < 0 on [9/20, 1/2), so G(s,p) ≥ Gh(p)
  on [9/20, 1/2].

## 4. Closed forms

Define, for x ∈ [0,1], y ≥ 0 and the indicated parameter ranges:

- J_a(x) = ∫_0^x u^a e^u du = x^{a+1}/(a+1) · ₁F₁(a+1; a+2; x), for a > −1;
- E_a(y) = e^y Γ(a, y) = ∫_0^∞ (y + t)^{a−1} e^{−t} dt, with E_a(0) = Γ(a) for a > 0;
- K_b(z) = ₁F₁(1; b; z) = Σ_{k≥0} z^k/(b)_k, for b ≥ 1;
- P_q[f] = E[(d + sR)^q f(R)];
- Q_q[f] = E[(d − cR)^q f(R); R < x] (the inner piece);
- O_q[f] = E[(cR − d)^q f(R); R > x] (the outer piece).

**Lemma 4 (moment closed forms).** For j = 0, 1, 2:

1. *(Tricomi form, d > 0, s > 0)* E[(d + sR)^q R^j] = d^q y^{j+1} j! U(j+1, q+j+2, y), with y = d/s.
   At s = 0 it equals d^q j!.
2. *(Gamma form)* Substitute v = y + R. Then P_q[f] = s^q e^y ∫_y^∞ v^q f(v − y) e^{−v} dv, so a
   polynomial f(R) = Σ_k α_k (R + y)^k gives P_q[f] = s^q Σ_k α_k E_{q+1+k}(y).
3. *(minus side)* With u = x − R inside and t = R − x outside,
   Q_q[f] = c^q e^{−x} ∫_0^x u^q f(x − u) e^u du and O_q[f] = c^q e^{−x} ∫_0^∞ t^q f(x + t) e^{−t} dt.
   For monomials these reduce to J_{q+k}(x) and Γ(q+1+k):
   - E[|d − cR|^q R^j; R < x] = c^q e^{−x} Σ_i C(j,i) x^{j−i} (−1)^i J_{q+i}(x);
   - E[|d − cR|^q R^j; R > x] = c^q e^{−x} Σ_i C(j,i) x^{j−i} Γ(q+1+i).

*Proof.* Part 1: substitute R = yt and use U(a,b,z) = Γ(a)^{−1} ∫_0^∞ e^{−zt} t^{a−1} (1+t)^{b−a−1} dt.
Parts 2 and 3 are direct substitutions. ∎

In the certifier, the polynomials are expanded in the variable that makes every coefficient explicit,
using 2c − x = N/c and y + 2s = N/s. Inside, R − 2c = −(N/c + u). Outside, R − 2c = t − N/c. On the
plus side, R − 2s = v − N/s.

**Proposition B (formulas evaluated by the certifier).** For p > 1:

  M = s·P_p[1] + c^{p+1} e^{−x} (J_p(x) + Γ(p+1)),
  a = P_{p−1}[R − 2s] − c^{p−1} e^{−x} ( (N/c) J_{p−1}(x) + J_p(x) + Γ(p+1) − (N/c) Γ(p) ),
  (p−1) b = T₊ + T₋ + S,

where

  T₋ = c^{p−1} e^{−x} [ (p−1)(2s²/(cN)) J_{p−1} + (p−1)(s²/N²) J_p + (N/c²) Γ(p)
                        − (p−1)(2s²/(cN)) Γ(p) + (p−1)(s²/N²) Γ(p+1) ],

and the plus-side terms take one of two forms.

**Mode "L1" (Gamma form).**

- s·P_p[1] = s^{p+1} E_{p+1}(y)
- P_{p−1}[R − 2s] = s^{p−1} (E_{p+1}(y) − (N/s) E_p(y))
- T₊ = s^{p−1} [ (p−1)(c²/N²) E_{p+1}(y) − (p−1)(2c²/(sN)) E_p(y) + (N/s²) e^y Γ(p) ]
- S = d^{p−1} [ (N/c²) K_p(−x) − (N/s²) K_p(y) ]

**Mode "small_s" (Tricomi form).**

- s·P_p[1] and P_{p−1}[R − 2s] = P_{p−1}[R] − 2s P_{p−1}[1] come from Lemma 4.1.
- T₊ = (p−1) s [ (c²/N²) P_{p−2}[R²] − (4sc²/N²) P_{p−2}[R] + (4s²c²/N² + 1) P_{p−2}[1] ]
- S = d^{p−1} (N/c²) K_p(−x)

*Proof.* Expand the polynomials as above and apply Lemma 4. The only terms that are singular as p → 1
are those carrying Γ(p−1, y) (plus side, exponent p − 2) and J_{p−2}(x) (inner piece). Both are
multiplied by (p − 1). Use

  (p−1) Γ(p−1, y) = Γ(p) − y^{p−1} ₁F₁(p−1; p; −y),  (p−1) J_{p−2}(x) = x^{p−1} ₁F₁(p−1; p; x),

then Kummer's transformation, which gives e^{y} ₁F₁(p−1;p;−y) = K_p(y) and e^{−x} ₁F₁(p−1;p;x) = K_p(−x).
Together with (sy)^{p−1} = (cx)^{p−1} = d^{p−1}, this gives S. In the Gamma form the plus-side
singular part enters S with a minus sign. In the Tricomi form it has no singular part, because d ≥ 9/10
there. Every term in Proposition B is finite and continuous on {p ≥ 1}. ∎

**Remark (the corner (s,p) = (1/2, 1)).** In terms of w = |Z|, the singular mass of
(p−1)E[|Z|^{p−2}W²] at w = 0 sits on both sides of Z = 0:

- the inner piece (0 < Z < d) contributes d^{p−1}(N/c²)K_p(−x) → 2d^{p−1};
- the plus piece (Z > d) contributes 2(1 − d^{p−1}), coming from the Γ(p) term minus the plus part of S;
- the outer piece (Z < 0) contributes 2.

So M_ss is continuous at the corner (G_ss(1/2,1) = 4), although the quantity d^{p−1} has no limit there.
Dropping the inner piece loses 2d^{p−1}. On the line p = 1 that
is a loss of 2 out of 4, which leaves only the bound G_ss ≥ 4 − 4 = 0 near the corner, with no margin.
The certifier therefore uses the exact cancelled form. S = d^{p−1} · [bracket], where the bracket is
O(d) and d^{p−1} ∈ [0,1] is enclosed separately, using the form d^{p−1} ₁F₁(p−1;p;b·d).

## 5. Rigour of the enclosures (`local_cert.py`)

Each box [s₀,s₁] × [p₀,p₁] has exact rational endpoints. An enclosure of G_ss valid for all (s,p) in the
box with p > 1 is built from the following ingredients.

1. **The p-dependence.** A single Arb ball P ⊇ [p₀, p₁] is formed from arb(fmpq) endpoints. All Arb
   operations (Γ, Γ(·,y), ₁F₁, U, pow, exp, log) return balls containing the exact result for every
   parameter value in the input balls. The computation runs at 128-bit precision with outward rounding.
2. **Rational functions of s.** Independent tight intervals are used for s, c, d and N (N is decreasing
   on [0,1/2]), and for x = (1−2s)/(1−s) and y = 1/s − 2, both decreasing. Each interval comes from
   exact endpoint values. Treating these as independent only enlarges the enclosure.
3. **Special-function pieces in s.** Exact monotone hulls are used: the enclosure is the union of the
   balls computed at the two (or four) exact endpoint points. The monotonicity facts are:
   - J_a(x) is increasing in x (positive integrand);
   - E_a(y) = ∫_0^∞(y+t)^{a−1}e^{−t}dt is monotone in y for each fixed a;
   - K_b(z) is increasing in z for b ≥ 1. Its derivative is (1/b)₁F₁(2;b+1;z), which equals e^z at b = 1
     and is a positive Euler integral for b > 1;
   - E[(d + sR)^q R^j] is monotone in d and, separately, in s for each fixed q. On the rectangle
     [1−2s₁, 1−2s₀] × [s₀, s₁], which contains the curve d = 1 − 2s, it therefore lies between its corner
     values. This handles s = 0 exactly, where the value is d^q j!.
4. **Powers d^{p−1} with d → 0.** For d ∈ [d₀,d₁] ⊂ [0,1] and p − 1 ∈ [e₀,e₁] ⊂ [0,∞), the power is
   increasing in d and nonincreasing in the exponent. So d^{p−1} ∈ [d₀^{e₁}, d₁^{e₀}], with 0^{e} = 0 for
   e > 0 and the value 1 at e₀ = 0. This interval is valid for every p > 1 in the box.
5. **K_b(z).** Summed explicitly to K = 48 + 8⌈|z|⌉ terms in ball arithmetic. For b > 1/2, (b)_k ≥ (k−1)!/2,
   so the tail is at most (2|z|^K/(K−1)!)·K/(K−|z|). This bound is added as a radius.
6. **Acceptance.** A box is accepted only if Arb certifies sign·G_ss > thr, so an undecided comparison
   counts as rejection. A rejected box is bisected, in s if w_s·Δs ≥ w_p·Δp and in p otherwise. Boxes
   reaching width 2^{−22} in both directions are recorded as failures, and any failure makes the run
   exit non-zero.
7. **Coverage.** Leaves tile the region by construction. The exact rational total area of the leaves is
   checked against the region area, and the sorted leaf list is hashed (sha256) into the certificate.

**Self-tests.** These run before every certification and abort it on any failure.

- (i) The seven embedded reference values of G, and G0(1.3), Gh(1.3), all to within 10^{−14}.
- (ii) The closed forms of Proposition A(3) at p ∈ {8/5, 19/10, 39/20, 41/20, 5/2, 7, 12} and of
  Proposition A(4) at p ∈ {1, 21/20, 3/2, 7/4, 5/2}. The Arb balls must overlap, and A must be exactly 0
  at s = 0.
- (iii) The Tricomi and Gamma forms agree at s = 1/4, to radius 10^{−20}.
- (iv) In 40 random sub-boxes per region, with 3 random interior points per box, all point enclosures of
  M, a, (p−1)b, A and G_ss lie inside the box enclosures (600 of 600 per region).
- (v) The same containment holds for boxes at all four corners of each region, including the singular
  corner (1/2, 1) of L1.

## 6. Historical independent-formula diagnostics (`formula_check.py`, mpmath, 50 digits)

| check | content | worst discrepancy |
|---|---|---|
| C1 | closed-form M vs direct quadrature of E\|Z_s\|^p (7 reference points) | 2.673e-51 (rel.) |
| C1 | G vs embedded reference values (incl. G0, Gh at 1.3) | 4.942e-15 |
| C2 | conditioning M_s (conditioning integrals) vs finite differences of quadrature M | 1.541e-50 (rel.) |
| C2 | conditioning M_ss vs finite differences of quadrature M (singularity removed by u = v^{1/(p−1)}) | 8.533e-51 (rel.) |
| C3 | Lemma 3 closed form G_s vs finite differences of quadrature G | 5.346e-51 |
| C3 | Lemma 3 closed form G_ss vs finite differences of quadrature G | 2.138e-49 |
| C3 | direct derivative combination (M_ss M − M_s²)/(pM²) − (N''N − N'²)/(2N²) vs Lemma 3 | 8.553e-50 |
| C3 | Gamma form (L1 mode) vs Tricomi form | 1.283e-49 |
| C4 | central FD of G_s at s = 0 (Z_s for s < 0 via hypoexponential density) | 1.928e-17 |
| C4 | FD G_ss(0,p) vs ((p−1)/p)(1/F − 1) | 2.473e-32 |
| C4 | FD G_s(1/2,p) | 0.000e+00 |
| C4 | FD G_ss(1/2,p) vs 4(p−2)²/p | 1.438e-22 |
| C5 | Arb point enclosures (local_cert.py) contain mpmath values (15 points) | 0 of 15 not contained |

C2 and C3 use the 24 points s ∈ {0.001, 0.02, 0.05, 0.3, 0.45, 0.49} × p ∈ {1.05, 1.9, 2.5, 10}. C4 uses
p ∈ {1.05, 1.3, 1.6824572, 1.9, 2.5, 5, 10, 12}. In the notation above, (1/2)G_ss(0,p) = ((p−1)/(2p))(1/F − 1)
and (1/2)G_ss(1/2,p) = 2(p−2)²/p are confirmed by FD, by the closed forms, and rigorously at
the self-test points.

## 7. Certified results

| region | box (exact) | claim | boxes evaluated | leaves | failures | min certified bound on the leaves | runtime | certificate |
|---|---|---|---|---|---|---|---|---|
| L1 | s∈[9/20,1/2], p∈[1,7/4] | d2G/ds2 > 6/125 | 50214 | 25182 | 0 | 0.048000 | 10 s | `cert_L1.json` |
| L2 | s∈[0,1/20], p∈[8/5,39/20] | d2G/ds2 > 29/2000 | 129875 | 65025 | 0 | 0.014500 | 56 s | `cert_L2.json` |
| L3 | s∈[0,1/20], p∈[41/20,7] | d2G/ds2 < -31/2000 | 57496 | 28998 | 0 | 0.015500 | 25 s | `cert_L3.json` |

The historical certificate format records a source hash, python-flint 0.9.0, 128-bit Arb, minimum box width 2^-22, and the sha256 of its sorted leaf list. The worst leaf of each run lies where |G_ss| is numerically smallest: at the corner (9/20, 7/4) for L1, and on the edges p = 39/20 and p = 41/20 for L2 and L3.

**Theorem L.** The following hold.

- **L1.** For all s ∈ [9/20, 1/2] and 1 < p ≤ 7/4: G_ss(s,p) > 6/125.
- **L2.** For all s ∈ [0, 1/20] and 8/5 ≤ p ≤ 39/20: G_ss(s,p) > 29/2000.
- **L3.** For all s ∈ [0, 1/20] and 41/20 ≤ p ≤ 7: G_ss(s,p) < −31/2000.

Each region contains the corresponding manuscript partition region: p* + 1/20 < 7/4 and p* − 1/20 > 8/5. For L3 the upper
end p = 7 joins the analytic large-p argument. The thresholds sit just
below the numerical extrema of |G_ss| on the regions:

- 0.04903 at (9/20, 7/4);
- 0.01461 at (0, 39/20);
- 0.01565 at (0, 41/20).

The certified constants are therefore essentially the true minima. Together with Proposition A, Section 3
gives the local monotonicity statements needed by the manuscript's partition.

*Reproduce from the repository root, using fresh output names.* Run `uv run --no-project --with python-flint==0.9.0 python3 code/local/local_cert.py L1 --threshold 6/125 --out NEW_L1.json`,
and the analogous commands with `L2 --threshold 29/2000 --out NEW_L2.json` and `L3 --threshold 31/2000 --out NEW_L3.json`. Each run takes about
a minute or less on one core.
