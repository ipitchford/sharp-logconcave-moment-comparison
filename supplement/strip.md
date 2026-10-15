# The strip around p = 2: Lemmas S1 and S2

**Evidence status and paths.** The mathematical derivations below form this supplement. Numerical tables, timings and validation outcomes are historical development summaries, not a fresh independent production replay. Their rounded decimal margins are diagnostics, not outward-rounded certificates. Proof-critical scripts are in `code/strip/` and outputs in `certificates/strip/`; current hashes, exact exported bounds and acceptance results belong to those files and the release replay receipt. Auxiliary samples, formula checks, negative controls and wider-region runs do not replace exhaustive production enclosures. No unshipped development log is a dependency.

Status: **both lemmas certified** with Arb (python-flint 0.9.0) on the partition values
δ = σ0 = 1/20. The shipped strip evidence consists only of the fresh production S1 and S2 certificates. Extended-strip and alternate-parameter results below are historical development summaries; their certificates are not shipped or required for the production proof.

| Certificate | Statement | Leaves | Max. upper bound | Time |
|---|---|---|---|---|
| `S1_cert.json` | H(s,p) < 0 on [1/20, 1/2] × [39/20, 41/20] | 157 | −1.539·10⁻⁵ | 3.4 s |
| `S2_cert.json` | ∂_p∂_s²G < 0 on [0, 1/20] × [39/20, 41/20] | 64 | −4.000·10⁻² | 1.0 s |

Here H(s,p) = ∂_p[G(s,p) − G0(p)]. "Max. upper bound" is the largest certified upper bound over
the leaf boxes; every leaf has a strictly negative upper bound. For comparison, diagnostic
point values on the production regions (mpmath) are H ≈ −3.53·10⁻⁴ at (s,p) = (1/20, 39/20) and
∂_p∂_s²G ≈ −0.277 at the same corner.

**Theorem S.** For every p ∈ [39/20, 41/20] and every s ∈ [0, 1/2], the difference G(s,p) − G0(p)
has the sign of 2 − p. It is strictly positive for p < 2 and strictly negative for p > 2 when
s > 0, and it vanishes when s = 0 or p = 2. Hence T_B holds for p ∈ [39/20, 2) and T_C holds for
p ∈ (2, 41/20]. The optional wider-range calculation discussed below is not part of the shipped certificate evidence.

---

## 1. Notation

E, E′, R are i.i.d. standard exponentials, X = E − 1 and Y = E′ − 1. For s ∈ [0, 1/2] let

- d = 1 − 2s, c = 1 − s, r₀ = d/c ∈ [0, 1];
- Z = Z_s = sX − cY, W = X + Y = ∂_s Z_s, N(s) = s² + c² = E Z²;
- M(s,p) = E|Z|^p, G(s,p) = p⁻¹ log M(s,p) − ½ log N(s), and G0(p) = G(0,p);
- q = p − 2 and L = log|Z|.

The three basic functions are

  Φ_A(s,q) = E|Z|^{2+q},  Φ_B(s,q) = E[|Z|^q Z W],  Φ_C(s,q) = E[|Z|^q W²].

They are used for q ∈ [−1/5, 10]; in particular q > −1 throughout.

## 2. Conditioning representation

**Lemma 1.** First let 0 < s ≤ 1/2 and D = sE − cE′, so that Z = d + D. Put σ = + and R = D/s on
{D > 0}, and σ = − and R = −D/c on {D < 0}. Then:

1. P(σ = +) = s;
2. R ~ Exp(1) and R is independent of σ;
3. conditionally on (σ, R), W has the law of (R − 1) + (T − 1) with T ~ Exp(1) independent.

In particular E[W | σ,R] = R − 1 and E[W² | σ,R] = (R − 1)² + 1. Consequently, for every
measurable φ with E|φ(Z)|(1 + W²) < ∞,

  E[φ(Z) W^j] = s ∫₀^∞ φ(d + sr) w_j(r) e^{−r} dr + c ∫₀^∞ φ(d − cr) w_j(r) e^{−r} dr,

with w₀ = 1, w₁ = r − 1 and w₂ = (r − 1)² + 1.

*Proof.* On {D > 0} substitute E = r + ct and E′ = st with r, t > 0. Then D = sr, the
Jacobian is s, and the density e^{−E−E′} becomes s e^{−r} e^{−t}, because c + s = 1. Also
W = (r − 1) + (t − 1). On {D < 0} substitute E = ct and E′ = r + st. Then D = −cr, the
Jacobian is c, the density becomes c e^{−r} e^{−t}, and W = (r − 1) + (t − 1). These two
pieces cover the full measure, which gives all three claims. For s = 0 we have D = −E′ < 0
almost surely, R = E′, T = E, and the same formulas hold. ∎

For 0 < s ≤ 1/2, Z has density f_s(z) = e^{−ρ(z)}, where ρ(z) = (z − d)/s for z > d and
ρ(z) = (d − z)/c for z < d. Given Z = z, the side's exponential is R = ρ(z).
At s = 0, use directly f_0(z) = e^{z−1} for z < 1 and zero for z > 1; no division by s or positive-side conditioning is used.

**Lemma 2 (uniform moment bound).** Let ψ ∈ {1, |r − 1|, (r − 1)² + 1}, with constants
C_ψ = 1, 1, 2 respectively. For every measurable g ≥ 0 on (0, ∞) and every s ∈ [0, 1/2],

  E[ψ(R) g(|Z|)] ≤ 2 C_ψ e^{1/2} ∫₀^∞ g(x) e^{−x/2} dx.

*Proof.* By Lemma 1, E[ψ(R)g(|Z|)] = ∫_ℝ ψ(ρ(z)) e^{−ρ(z)} g(|z|) dz. The three functions
satisfy ψ(r)e^{−r} ≤ C_ψ e^{−r/2} on r ≥ 0. To see this, note that e^{−r/2} ≤ 1, that
|r − 1| e^{−r/2} ≤ max(1, 2e^{−3/2}), and that ((r−1)²+1)e^{−r/2} is maximised at r = 0, where it
equals 2; its other critical point r = 3 + √3 gives 1.40. Since s, c ≤ 1, ρ(z) ≥ |z − d| ≥ |z| − 1. ∎

## 3. Regularity and derivative identities

**Lemma 3.** On [0, 1/2] × [−1/5, 10] the following hold:

(a) Each Φ_V is C^∞ in q, with ∂_q^j Φ_V = E[V |Z|^q L^j], where V = Z², ZW or W². These
    functions are jointly continuous in (s,q).

(b) ∂_s ∂_q^j Φ_A = ∂_q^j[(2+q)Φ_B] and ∂_s ∂_q^j Φ_B = ∂_q^j[(1+q)Φ_C], for j ≥ 0.

(c) In the conditioned form, Φ_B = E[|Z|^q Z (R−1)] and Φ_C = E[|Z|^q((R−1)²+1)].

Hence M = Φ_A, M_s = pΦ_B and M_ss = p(p−1)Φ_C, and with q = p − 2:

  G_s = Φ_B/Φ_A − N′/(2N),
  G_ss = (p−1)Φ_C/Φ_A − p(Φ_B/Φ_A)² − ½(log N)″,
  G_p = ∂_q[log Φ_A/(2+q)],
  ∂_s G_p = ∂_q(Φ_B/Φ_A),
  ∂_s² G_p = ∂_q[(1+q)Φ_C/Φ_A − (2+q)(Φ_B/Φ_A)²] =: G_pss.

*Proof.* (c) follows from Lemma 1 by the tower property.

(a) For q in a compact subinterval [q_a, q_b] ⊂ (−1, ∞), the derivatives ∂_q^j(|Z|^q V) =
|Z|^q L^j V are dominated by (|Z|^{q_a} + |Z|^{q_b}) |L|^j |V|, and the same holds with j + 1.
The expectation of this bound is finite by Lemma 2, applied after conditioning on (σ, R). Take
ψ = 1 for V = Z², and ψ = (r−1)² + 1 for V = ZW and V = W²; the latter works because
E[|W| | σ,R] ≤ (E[W² | σ,R])^{1/2} ≤ (R−1)² + 1. Take g(x) = x^m(x^{q_a} + x^{q_b})|log x|^j
with m = 2, 1, 0, and note that m + q_a > −1. So differentiation under E is justified.

For joint continuity, write each expectation in the form of Lemma 1. On the negative side,
substitute u = r − r₀(s), so that |d − cr| = c|u| with u ≥ −r₀ ≥ −1 and c ∈ [1/2, 1]. The
integrand is then continuous in (s,q) for fixed u. It is dominated, independently of s, by
C(1+u²)(|u|^{m+q_a} + |u|^{m+q_b})(1 + |log|u||)^j e^{1−u}, which is integrable on [−1, ∞).

On the positive side, for s ≥ s₀ > 0 we have s₀r ≤ x = d + sr ≤ 1 + r. The integrand is then
dominated by C(1+r²)(r^{m+q_a} + (1+r)^{m+q_b})(1 + |log r| + log(1+r))^j e^{−r}. For
s ∈ (0, s₀] with s₀ < 1/2, we have x ≥ 1 − 2s₀ > 0, so the r-integral is bounded and the
factor s makes the positive side tend to 0 as s → 0. Dominated convergence gives continuity on
[0, 1/2] × [q_a, q_b].

(b) The map s ↦ |Z_s|^{2+q} L^j is C¹ along each sample path, because 2 + q > 1. Its derivative
is W |Z|^q Z((2+q)L^j + jL^{j−1}) = W · ∂_q^j[(2+q)|Z|^q Z]. For q ∈ [−1/5, q_b],
|z|^{1+q}(|log|z||^j + |log|z||^{j−1}) ≤ C_j(1 + |z|)^{2+q_b}, and |Z_s| ≤ |X| + |Y| for
s ∈ [0,1]. So the derivative is dominated, uniformly in s, by
C_j(1 + |X| + |Y|)^{2+q_b}|W|, which is integrable. This proves the first identity.

For the second identity, z ↦ |z|^q z L^j is absolutely continuous on ℝ for q > −1, with
derivative |z|^q((1+q)L^j + jL^{j−1}), which lies in L¹_loc. Since t ↦ Z_t is affine with slope
W, for each ω and s₁ < s₂ we get

  W(|Z_{s₂}|^q Z_{s₂} L_{s₂}^j − |Z_{s₁}|^q Z_{s₁} L_{s₁}^j) = ∫_{s₁}^{s₂} W² |Z_t|^q((1+q)L_t^j + jL_t^{j−1}) dt,

where L_t = log|Z_t|. The double integral E∫ W²|Z_t|^q(…) dt is absolutely convergent, by
Lemma 2 applied with ψ = (r−1)² + 1 after conditioning, uniformly in t. So Fubini gives
∂_q^jΦ_B(s₂) − ∂_q^jΦ_B(s₁) = ∫_{s₁}^{s₂} ∂_q^j[(1+q)Φ_C](t) dt. The integrand is continuous
by (a), which gives the derivative.

The formulas for G follow by differentiating log Φ_A/(2+q) − ½ log N. For example,
∂_s[log Φ_A/(2+q)] = Φ_B/Φ_A by (b) with j = 0, and ∂_s ∂_q[log Φ_A/(2+q)] is obtained from
(b) with j = 0, 1. The resulting expression equals ∂_q(Φ_B/Φ_A). In the same way,
∂_s²∂_q[·] = ∂_q[(1+q)Φ_C/Φ_A − (2+q)(Φ_B/Φ_A)²]. All mixed derivatives are computed from (b),
so no abstract symmetry-of-derivatives theorem is needed. ∎

**Corollary 4 (values on the boundary lines).**

(i) G(s,2) = 0 for all s, because M(s,2) = E Z² = N(s); hence G_s(s,2) = G_ss(s,2) = 0.

(ii) G_s(0,p) = 0 for p ∈ [9/5, 12] (indeed for all p > 1). At s = 0 we have Z = 1 − R
     almost surely, so
     Φ_B(0,q) = E[|Z|^q Z(R−1)] = −Φ_A(0,q). Then G_s(0,p) = −1 − N′(0)/(2N(0)) = −1 + 1 = 0.

(iii) With H(s,p) := G_p(s,p) − G0′(p), we have H(0,p) = 0 and
      ∂_s H(0,p) = ∂_q(Φ_B/Φ_A)(0,q) = ∂_q(−1) = 0.

The derivative formulas were checked against independent mpmath finite differences (§7) at
s ∈ {0, 0.02, 0.2, 0.5} and p ∈ {1.96, 2, 2.04}.

## 4. Rigorous evaluation of Φ_A, Φ_B, Φ_C as power series in q

Every Φ_V(s, ·) is enclosed as a truncated Taylor series in q of length 3 (orders 0, 1, 2),
using python-flint `arb_series`. Coefficient j encloses ∂_q^jΦ_V/j! simultaneously for all
s in an s-ball and all q in a q-ball. Sums, products, quotients, log and exp of such series
again enclose the Taylor coefficients of the corresponding functions, because each coefficient
of the result is a polynomial or rational expression in the lower-order input coefficients,
evaluated in ball arithmetic. Φ_V is split into a negative side, the event {σ = −} = {Z ≤ d},
and a positive side, {σ = +} = {Z > d}.

### 4.1 Negative side: closed form, exact in q (Lemma 5)

Substitute z = d − cr in the negative-side integral of Lemma 1 and split it at z = 0.

*Part {Z < 0}.* Put u = cr − d and use r − 1 = (u − s)/c:

- A: e^{−r₀} c^{3+q} Γ(3+q);
- B: −e^{−r₀}[c^{2+q}Γ(3+q) − s c^{1+q}Γ(2+q)];
- C: e^{−r₀} c^{−2}[c^{3+q}Γ(3+q) − 2s c^{2+q}Γ(2+q) + (s²+c²)c^{1+q}Γ(1+q)].

*Part {0 < Z < d}.* Put z = dt, so r = r₀(1−t), and expand e^{−r₀(1−t)} = e^{−r₀}Σₙ r₀ⁿtⁿ/n!:

- A: e^{−r₀} d^{3+q} Σₙ (r₀ⁿ/n!) · 1/(3+q+n);
- B: e^{−r₀} d^{2+q} Σₙ (r₀ⁿ/n!) · [(r₀−1)/(2+q+n) − r₀/(3+q+n)];
- C: e^{−r₀} d^{1+q} Σₙ (r₀ⁿ/n!) · [(α²+1)/(1+q+n) − 2αr₀/(2+q+n) + r₀²/(3+q+n)],
  where α = r₀ − 1.

As a check, at s = 0 the A formula gives F(p) = e^{−1}(Γ(p+1) + ∫₀¹ u^p e^u du).

The series is truncated after n = N = 40. Since r₀ ≤ 1, and the j-th Taylor coefficient in q
of 1/(β+q+n) has modulus at most 1 once β + q + n ≥ 1, the tail adds at most
Σ|wᵢ| · 2/(N+1)! to every coefficient. Here wᵢ are the weights in the brackets, with
Σ|wᵢ| ≤ 1, 2 and 5 for A, B and C.

Powers are computed as c^{β+q} = exp((β+q) log c) on series. The Taylor coefficients of
1/(β+q+n) are (−1)^j/(β+q+n)^{j+1}.

*s-balls reaching s = 1/2.* Then d ∈ [0, D] with D ≤ 1/10, and the {0 < Z < d} part is replaced
by its bound. Its j-th q-derivative is ∫₀^d z^{m+q} log^j z · w · e^{−(d−z)/c} dz, where
|w| ≤ ψ_max = 1 for A and B. In modulus it is at most

  ψ_max ∫₀^D z^{m+q_lo}|log z|^j dz = ψ_max D^{β+1} Σ_{i=0}^{j} [j!/(j−i)!] (−log D)^{j−i}/(β+1)^{i+1},

with β = m + q_lo > −1. This uses z^q ≤ z^{q_lo} for z ≤ 1, e^{−(d−z)/c} ≤ 1, and the fact that
the integrand is nonnegative, so the integral is monotone in the upper limit. Only A and B
are needed there, because S1 uses only Φ_A and Φ_B.

### 4.2 Gamma series on a q-ball (Lemma 6)

For a q-box [q_lo, q_hi] with centre q_c and radius ρ, set z_c = β + q_c and
r = min(9/10, z_c − 1/20). The exact Taylor series of Γ(z_c + t) of length L = 24 is computed
by Arb at the exact rational centre. It is then recentred at u ∈ [−ρ, ρ] by ball Horner
evaluation.

The neglected tail Σ_{n≥L} g_n binom(n,j) u^{n−j} is bounded using the Cauchy estimate
|g_n| ≤ M/rⁿ with M = max(Γ(z_c − r), Γ(z_c + r)). This M is valid because
|Γ(w)| ≤ Γ(Re w) for Re w > 0, and Γ is log-convex on (0, ∞), so its maximum on
[z_c − r, z_c + r] ⊂ [1/20, ∞) is attained at an endpoint. With θ = ρ/r ≤ 1/10, the tail bound
is (M/r^j) binom(L,j) θ^{L−j}/(1 − 2θ). This construction replaces `arb_series.gamma` on wide
balls, which turned out to overestimate badly; `arb.digamma` on balls overestimates badly too.

### 4.3 Positive side: Taylor polynomial in q with rigorous remainder (Lemma 7)

Put x = d + sr = 1 + s(r − 2), v_A = x², v_B = x(r−1) and v_C = (r−1)² + 1. Then

  Φ_V⁺(s,q) = s ∫₀^∞ v_V(r) x^q e^{−r} dr.

On a q-box with centre q_T and half-width a ≤ 1/20, write q = q_T + u. Then

  Φ_V⁺(s, q_T + u) = Σ_{k<K} P_k u^k/k! + ρ_V(u),  where P_k = s ∫₀^∞ v_V x^{q_T} log^k x · e^{−r} dr,

and K = 12. The remainder satisfies ρ_V^{(j)}(u) = E⁺[V |Z|^{q_T} L^j T_{K−j}(uL)], where
T_n(y) = e^y − Σ_{i<n} yⁱ/i! obeys |T_n(y)| ≤ |y|ⁿ e^{|y|}/n!. Together with
e^{|u||L|} ≤ |Z|^a + |Z|^{−a} and Lemma 2 this gives

  |ρ_V^{(j)}(u)| ≤ (|u|^{K−j}/(K−j)!) · B_V,
  B_V = 2C_ψ e^{1/2} Σ_{σ=±1} [ K!/(m+q_T+σa+1)^{K+1} + ∫₁^∞ x^{m+q_T+σa} log^K x · e^{−x/2} dx ].

The integral from 1 to ∞ is computed by `acb.integral` on [1, 400]. On [400, ∞) we use
log x ≤ x, which gives the bound 2^{β+K+1} Γ(β+K+1, 200).

The moments P_k are rigorous Arb integrals:

- *Middle range.* On [r_lo, T] with T = 120 + 20⌈max(q_T, 0)⌉ we use `acb.integral`, with s
  entered as a ball. The integrand calls `log(analytic=True)`, so it returns a non-finite value
  whenever the input ball meets the branch cut x ≤ 0. This is the documented contract for
  rigorous integration. Arb's Gauss–Legendre error bounds use enclosures of the integrand on
  complex ellipses, so the output encloses the integral for every s in the ball at once.
- *Tail.* For r ≥ T ≥ 2 we have 1 ≤ x ≤ r, 0 ≤ log x ≤ r, |v_V| ≤ r² and x^{q_T} ≤ r^{max(q_T,0)}.
  So the tail is at most s·Γ(k + 3 + max(q_T,0), T).
- *Head, s-balls reaching 1/2 only.* Here r_lo = max(4·width, 10⁻⁶). On [0, r_lo] the integrand
  is at most x^{m+q_T}|log x|^k, with x ∈ (0, X] and X = (1 − 2s_lo) + s_hi r_lo. Its supremum
  g(X) is explicit, since x^μ|log x|^k increases on (0, e^{−k/μ}], has its global maximum
  (k/μ)^k e^{−k} on (0, 1], and increases on [1, ∞). The head is therefore at most
  s_hi · r_lo · g(X).

Point evaluations at exact rational s use tolerance 2^{−100}; s-ball evaluations use 2^{−40}.
Non-finite results, or imaginary parts that exclude 0, raise an error.

## 5. Box enclosures and the adaptive search

**S1 (mean value in s and in q).** Let S = [s_lo, s_hi] and Q = [q_lo, q_hi], with centres s_c
and q_c. For every (s,q) ∈ S × Q,

  H(s,q) ∈ H(s_c,q_c) + ∂_qH(s_c, Q)·(Q − q_c) + ∂_sH(S, Q)·(S − s_c).

The three ingredients are computed as follows:

- H(s_c,q_c) is coefficient 1 of the series log Φ_A/(2+q), at s = s_c and q = q_c, minus the
  same quantity at s = 0 (that is, G0′), using exact point data.
- ∂_qH(s_c, Q) is twice coefficient 2 of the same difference, taken over the q-ball Q.
- ∂_sH(S,Q) = ∂_q(Φ_B/Φ_A) (Lemma 3) is coefficient 1 of Φ_B/Φ_A over S × Q.

The mean-value form is valid because H(·,q) is C¹ in s and H(s,·) is analytic in q (Lemma 3).
This form matters: H is O(s²) while Φ_A and Φ_B vary at O(1) rates in s, so naive interval
evaluation overestimates by a factor of about 100 near s = 1/20.

**S2 (direct in s, mean value in q).** Here

  G_pss(S,Q) ⊂ G_pss(S, q_c) + G_ppss(S, Q)·(Q − q_c),

where G_pss is coefficient 1 and G_ppss/2 is coefficient 2 of the series
(1+q)Φ_C/Φ_A − (2+q)(Φ_B/Φ_A)². The term (log N)″ does not depend on q and drops out.

**Adaptive search.** A box is accepted when the Arb upper bound of its enclosure is < 0, and
undecided comparisons count as failures. Otherwise the box is bisected in s or in q, whichever
term of the enclosure is wider. The initial q-boxes have width at most 1/10, so |q − q_T| ≤ 1/20.

The program then checks, in exact rational arithmetic, that the leaves tile the region: every
leaf is contained in it, leaves have pairwise disjoint interiors, and the total area equals the
area of the region. Any failure raises `CertError` and the program exits with status 2. A
certificate file is never overwritten.

## 6. Proof of Theorem S

Let p ∈ [39/20, 41/20] and H(s,p) = ∂_p[G(s,p) − G0(p)].

1. **S2 ⇒ H < 0 on (0, 1/20].** By Lemma 3, ∂_s²H = G_pss, which is negative on
   [0, 1/20] × [39/20, 41/20] by `S2_cert.json`. By Corollary 4(iii), H(0,p) = ∂_sH(0,p) = 0.
   Taylor's formula with integral remainder therefore gives
   H(s,p) = ∫₀^s (s − t) ∂_t²H(t,p) dt < 0 for 0 < s ≤ 1/20.
2. **S1.** `S1_cert.json` gives H < 0 on [1/20, 1/2] × [39/20, 41/20]. So H < 0 on
   (0, 1/2] × [39/20, 41/20].
3. **Sign.** By Corollary 4(i), G(s,2) − G0(2) = 0, so G(s,p) − G0(p) = ∫₂^p H(s,t) dt. For
   s > 0 this is strictly positive when p < 2 and strictly negative when p > 2. For s = 0,
   G(0,p) = G0(p) by definition. ∎

*Second route (equivalent integration order).* Since G_ss(s,2) = 0 (Corollary 4(i)),
G_ss(s,p) = ∫₂^p G_pss(s,t) dt, which has the sign of 2 − p on [0, 1/20] by S2. Since
G_s(0,p) = 0 (Corollary 4(ii)),

  G(s,p) − G0(p) = ∫₀^s (s − v) G_ss(v,p) dv,

which has the sign of 2 − p for 0 < s ≤ 1/20. On [1/20, 1/2], S1 gives the sign through
∫₂^p H dt. Both routes need only S1, S2 and Corollary 4.

*Historical extension (not shipped).* Development records report the same argument on
p ∈ [19/10, 3] using optional extended certificates. Those outputs are excluded from this release and are not a reader dependency. This does not change the production partition, whose strip boundaries remain 39/20 and 41/20; a reader wishing to check the extension must generate fresh outputs with the optional arguments below.

*Further range (numerical, not certified).* Numerically, H < 0 holds on (0, 1/2] for every
tested p ∈ [1.85, 12]: H/s² ranges from about −0.004 at p = 12 to −0.29. G_pss < 0 also holds on
[0, 1/20] across the same range. Certifying p near 12 with the present code is expensive.
The H_q enclosure suffers from cancellation in (log Φ_A)″ when Φ_A ≈ 10⁸, and the S2 enclosure
suffers from cancellation of order p between the two terms of G_ss. The left limit is real:
H(1/2, p) changes sign near p ≈ 1.84, and the zero set of H crosses s ∈ [1/20, 1/2] for
smaller p. See the negative control in §7.

## 7. Historical validation diagnostics

1. **Finite differences** (`code/strip/strip_fdcheck.py`; mpmath, 40 digits). The
   formulas for M_p, M_s, M_ss, M_sp, M_ssp, G_p, G_s, G_ss and G_pss were compared against
   4th-order finite differences of M(s,p) = E|Z|^p and of G alone, at s ∈ {0, 0.02, 0.2, 0.5}
   and p ∈ {1.96, 2, 2.04}. Agreement is at the 10⁻¹⁵ to 10⁻²¹ level for the first and
   mixed derivatives. For G_pss from nested finite differences of G, agreement is 3·10⁻⁹ to
   6·10⁻¹⁰ in the interior. At s = 1/2 the one-sided nested finite differences are only
   accurate to 2·10⁻⁴, but the formula agrees there with the exact value
   G_pss(1/2,p) = 4(2q/p − q²/p²) to 12 digits. M itself was compared against an independent 2-D
   quadrature of E|s(E−1) − c(E′−1)|^p, with no conditioning, at four points, agreeing to
   10⁻¹³.
2. **embedded reference values.** All seven values of G pass to 10⁻¹³. They are computed by the
   certifier's closed forms at arbitrary p plus a direct Arb integral of the positive side.
3. **mpmath cross-checks** inside the certifier: G_p and G_pss at 14 points with p from 1.9 to 12
   agree to 10⁻¹⁰.
4. **Identities inside the certifier.** The certifier checks G_ss(s,2) = 0 at 4 points and
   G_s(0,p) = 0 at 3 points, both to 10⁻²⁵. It also checks Taylor-centre consistency: the
   enclosures over a q-box contain the point enclosures.
5. **Second parameter set.** With 192 bits, K = 16, N = 60 and L = 32, both lemmas certify again
   with the same leaf counts and the same worst bounds to 7 digits in the historical records. These alternate-parameter outputs are not shipped or part of the current replay claim.
6. **Negative control.** The same program, asked to certify S1 on p ∈ [9/5, 19/10] where H
   changes sign, stops with status 2 and writes no certificate:
   `CERTIFICATION FAILURE: S1: cannot certify box s=[2018083/5242880,504521/1310720] …
   enc=[+/- 1.12e-8]`. That box lies on the zero set of H near (s,p) = (0.385, 1.8). The
   program also stops immediately on any box whose whole enclosure is ≥ 0.
7. **Sampling check** (`code/strip/strip_sample_check.py`; historical logs not shipped). Random points in every
   leaf were evaluated by plain mpmath quadrature of the derivative formulas, with no closed
   forms and no Taylor expansion. Every value lay inside its leaf's certified enclosure:
   S1 had 471/471 points inside, S2 192/192, S1ext 1448/1448 and S2ext 817/817.

## 8. Rigour notes and dependencies

- **Arithmetic.** Region endpoints are exact rationals, converted to Arb balls containing them.
  All arithmetic is Arb ball arithmetic, so it rounds outward. Comparisons use `upper() < 0`,
  and an undecided comparison is treated as a failure. A box whose enclosure is entirely ≥ 0
  aborts the run, because the inequality is then false there.
- **Trusted components.** The proof relies on Arb's enclosure semantics for:
  - `acb.integral` with ball parameters and the `analytic` flag;
  - `arb_series` arithmetic, `log` and `exp`;
  - `arb_series.gamma` at an exact point;
  - `gamma`, `gamma_upper` and `exp`/`log` on real balls.
- **Truncations.** Every truncation has an explicit bound: the r₀-series tail (§4.1), the Gamma
  Cauchy tail (§4.2), the q-Taylor remainder (§4.3), the r-integral tail (§4.3), and the head
  piece near d = 0 (§4.3).
- **Code path for the analytic inputs.** The lemmas use only Lemma 3 (derivative identities),
  Lemma 1 (conditioning) and Corollary 4. These were proved above, and also checked numerically
  against finite differences of E|Z|^p.

## 9. Files and reproduction

Paths relative to the repository root:

- `code/strip/strip_cert.py`: the certifier; current source hashes are recorded in regenerated production certificates.
- `certificates/strip/S1_cert.json`, `certificates/strip/S2_cert.json`: production certificates.
- Alternate-parameter and extended-strip certificates are historical, omitted outputs, not files in the public package.
- `code/strip/strip_sample_check.py`, `code/strip/strip_fdcheck.py`, `code/strip/mpref.py`: optional mpmath checks. Historical logs are not shipped.

Commands from the repository root (use distinct fresh output names):

```
uv run --no-project --with python-flint==0.9.0 python3 code/strip/strip_cert.py --lemma S1 --out NEW_S1.json
uv run --no-project --with python-flint==0.9.0 python3 code/strip/strip_cert.py --lemma S2 --out NEW_S2.json
# Optional new experiments, not shipped outputs: --p-lo 19/10 --p-hi 3
# Alternate parameters, also requiring a fresh run: --prec 192 --k-taylor 16 --n-ser 60 --l-gamma 32
uv run --no-project --with mpmath python3 code/strip/strip_sample_check.py certificates/strip/S1_cert.json --per-leaf 3
uv run --no-project --with mpmath python3 code/strip/strip_fdcheck.py
```

Each certificate JSON records the script hash, the parameters, the validation results, the
statement, the region, and every leaf as exact rational boxes with an enclosure [lower, upper].
The program exits with status 0 only if certification succeeds; on failure it exits with
status 2 and writes no certificate.
