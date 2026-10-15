# LARGE-P: the regime $p \ge 7$ of $T_C$

**Evidence status and paths.** The mathematical derivations below form this supplement. Numerical tables, timings and validation outcomes are historical development summaries, not a fresh independent production replay. Their rounded decimal margins are diagnostics, not outward-rounded certificates. Proof-critical scripts are in `code/largep/` and outputs in `certificates/largep/`; current hashes, exact exported bounds and acceptance results belong to those files and the release replay receipt. Auxiliary samples, formula checks, negative controls and wider-region runs do not replace exhaustive production enclosures. No unshipped development log is a dependency.

Supplement to the manuscript: the notation and proof are given below.

$E, E', R$ are i.i.d. standard exponential random variables. For $s\in[0,\tfrac12]$ put
$c = 1-s$, $d = 1-2s$, $Z_s = s(E-1) - c(E'-1)$, $N(s) = s^2 + c^2 = \mathbb E Z_s^2$,
$M(s,p) = \mathbb E|Z_s|^p$, $F(p) = \mathbb E|E-1|^p = M(0,p)$, and
$G(s,p) = \tfrac1p\log M(s,p) - \tfrac12\log N(s)$, $G_0(p) = G(0,p) = \tfrac1p\log F(p)$.

Throughout, $\Gamma_p := \Gamma(p+1) = \mathbb E R^p$ and
$$J(x) := \int_0^x u^p e^u\,du \qquad (0\le x\le 1).$$

## Lemma LP

Let $p \ge 7$ be real and $s \in [0,\tfrac12]$. Then
$$M(s,p) \le N(s)^{p/2}\,F(p), \qquad\text{equivalently}\qquad G(s,p) \le G_0(p),$$
with equality only at $s = 0$. Quantitatively, for $s\in(0,\tfrac12]$,
$$N(s)^{p/2}F(p) - M(s,p) \;\ge\; \frac{s^2 (1-s)^{p-1}}{e}\,K(p) > 0,
\qquad
K(p) := \Big(\frac p2 - 2e + 2\Big)\Gamma_p - \frac{e(p-2)}{2} - \frac{e(e-2)}{p+1},$$
and consequently
$G_0(p) - G(s,p) \ge \dfrac{s^2(1-s)^{p-1}K(p)}{e\,p\,N(s)^{p/2}F(p)}$.

In particular the upper moment comparison $T_C$ holds for every $p\ge7$. The production partition uses $P=7$.

The proof is self-contained. It does not use $\partial_s G(0,p)=0$ or any other fact from the
partition; the first-order cancellation at $s=0$ is carried exactly through Steps 3–6.

## Proof

**Step 1 (conditioning representation).** Let $s\in(0,\tfrac12]$ and $X = sE - cE'$. For $t\ge 0$,
$$\mathbb P(X>t) = \mathbb E\big[e^{-(t+cE')/s}\big] = e^{-t/s}\,\frac{1}{1+c/s} = s\,e^{-t/s},
\qquad
\mathbb P(X<-t) = \mathbb E\big[e^{-(t+sE)/c}\big] = e^{-t/c}\,\frac{1}{1+s/c} = c\,e^{-t/c},$$
using $s + c = 1$. The mixture $Y$ that equals $sR$ with probability $s$ and $-cR$ with probability $c$
has the same two tails. Since both $X$ and $Y$ have densities, these tails determine their distribution
functions, so $X \overset{d}{=} Y$. Since $Z_s = X + (c - s) = X + d$,
$$M(s,p) = s\,\mathbb E(d+sR)^p + c\,\mathbb E|d - cR|^p. \tag{1}$$
(Here $d + sR \ge 0$ because $d\ge0$.) For $s = 0$, $M(0,p) = F(p)$ by definition.

**Step 2 (two-sided moment).** For $x \ge 0$,
$$\mathbb E|x-R|^p = \int_x^\infty (r-x)^p e^{-r}dr + \int_0^x (x-r)^p e^{-r}dr
= e^{-x}\Gamma_p + e^{-x}\!\int_0^x u^p e^{u}du = e^{-x}\big(\Gamma_p + J(x)\big), \tag{2}$$
substituting $r = x+v$ in the first integral and $u = x - r$ in the second. In particular
$F(p) = \mathbb E|1-R|^p = e^{-1}(\Gamma_p + J(1))$.

**Step 3 (rescaling).** Let $s\in(0,\tfrac12]$ and $\tau := s/c \in (0,1]$. Then $c = 1/(1+\tau)$,
$s = \tau c$, $d/c = 1-\tau$, and $N(s) = c^2(1+\tau^2)$. By (1) and (2),
$$\mathbb E(d+sR)^p = c^p\,\mathbb E(1-\tau+\tau R)^p,\qquad
\mathbb E|d-cR|^p = c^p\,\mathbb E|1-\tau-R|^p = c^p e^{\tau-1}\big(\Gamma_p + J(1-\tau)\big),$$
so
$$M(s,p) = \frac{c^p}{1+\tau}\Big[\tau\,\mathbb E(1-\tau+\tau R)^p + e^{\tau-1}\big(\Gamma_p + J(1-\tau)\big)\Big],
\qquad N(s)^{p/2}F(p) = c^p (1+\tau^2)^{p/2} e^{-1}\big(\Gamma_p + J(1)\big).$$
Put $A(\tau) := (1+\tau)(1+\tau^2)^{p/2}$. Multiplying by $e(1+\tau)c^{-p} > 0$,
$$e(1+\tau)c^{-p}\big[N(s)^{p/2}F(p) - M(s,p)\big] = S(\tau),$$
$$S(\tau) := A(\tau)\big(\Gamma_p + J(1)\big) - e\tau\,\mathbb E(1-\tau+\tau R)^p - e^{\tau}\big(\Gamma_p + J(1-\tau)\big). \tag{3}$$
It suffices to show $S(\tau) \ge \tau^2 K(p) > 0$ for $\tau\in(0,1]$ and $p\ge 7$; then
$N^{p/2}F - M \ge c^p\tau^2K(p)/(e(1+\tau)) = s^2c^{p-1}K(p)/e$, which is the quantitative claim.

**Step 4 (Jensen).** For $p \ge 1$ the map $t\mapsto t^p$ is convex on $[0,\infty)$, and
$1-\tau+\tau R$ is the convex combination of $1$ and $R$ with weights $1-\tau,\tau\in[0,1]$. Hence
$(1-\tau+\tau R)^p \le (1-\tau) + \tau R^p$ pointwise, and
$$\mathbb E(1-\tau+\tau R)^p \le (1-\tau) + \tau\Gamma_p. \tag{4}$$
Inserting (4) into (3):
$$S(\tau) \ge \Delta(\tau) := \big[A(\tau) - e^\tau - e\tau^2\big]\Gamma_p
+ \Big\{A(\tau)J(1) - e^{\tau}J(1-\tau) - e\tau(1-\tau)\Big\}. \tag{5}$$

**Step 5 (elementary inequalities).** For $\tau\in[0,1]$ and $p \ge 2$:

- (a) $A(\tau) \ge (1+\tau)\big(1 + \tfrac p2\tau^2\big) \ge 1 + \tau + \tfrac p2\tau^2$, and $A(\tau)\ge 1+\tau$.
  Bernoulli: $(1+y)^q \ge 1+qy$ for $y\ge0$, $q = p/2 \ge 1$.
- (b) $0 \le e^\tau - 1 - \tau \le (e-2)\tau^2$, because
  $e^\tau - 1-\tau = \tau^2\sum_{k\ge2}\tau^{k-2}/k! \le \tau^2\sum_{k\ge2}1/k! = (e-2)\tau^2$.
- (c) $0 \le J(1) \le e\int_0^1 u^p\,du = e/(p+1)$.
- (d) $\int_{1-\tau}^1 u^p e^u\,du \ge e^{1-\tau}\int_{1-\tau}^1 u^p\,du = e^{1-\tau}\,\dfrac{1-(1-\tau)^{p+1}}{p+1}$,
  because $e^u \ge e^{1-\tau}$ on $[1-\tau,1]$.
- (e) $1-(1-\tau)^{p+1} \ge (p+1)\tau - \tfrac{(p+1)p}{2}\tau^2$. Put $n = p+1\ge3$ and
  $g(\tau) = (1-\tau)^n$. Then $g'''(\tau) = -n(n-1)(n-2)(1-\tau)^{n-3} \le 0$ on $[0,1]$, so
  Taylor's theorem with Lagrange remainder gives $g(\tau) \le 1 - n\tau + \tfrac{n(n-1)}2\tau^2$.

**Step 6 (the bound $\Delta \ge \tau^2K$).** First bracket of (5). By (a) and (b),
$$A(\tau) - e^\tau - e\tau^2 \ge \Big(1+\tau+\frac p2\tau^2\Big) - \big(1+\tau+(e-2)\tau^2\big) - e\tau^2
= \Big(\frac p2 - 2e + 2\Big)\tau^2,$$
and $\Gamma_p > 0$. Second bracket of (5). Write $J(1-\tau) = J(1) - \int_{1-\tau}^1u^pe^u\,du$. Using
$A(\tau)J(1) \ge (1+\tau)J(1)$ (by (a) and $J(1)\ge0$), then (d), (b), (c) and (e):
$$\begin{aligned}
A(\tau)J(1) - e^\tau J(1-\tau)
&\ge (1+\tau-e^\tau)J(1) + e^{\tau}\int_{1-\tau}^1u^pe^u\,du\\
&\ge -(e-2)\tau^2\cdot\frac{e}{p+1} + e\,\frac{1-(1-\tau)^{p+1}}{p+1}
\;\ge\; -\frac{e(e-2)}{p+1}\tau^2 + e\Big(\tau - \frac p2\tau^2\Big).
\end{aligned}$$
Subtracting $e\tau(1-\tau) = e\tau - e\tau^2$, the terms linear in $\tau$ cancel exactly:
$$A(\tau)J(1) - e^\tau J(1-\tau) - e\tau(1-\tau) \ge -\frac{e(e-2)}{p+1}\tau^2 - \frac{e(p-2)}{2}\tau^2.$$
Adding the two brackets gives $\Delta(\tau) \ge \tau^2 K(p)$ for all $\tau\in[0,1]$, $p\ge2$.

**Step 7 ($K(p) > 0$ for $p \ge 7$).** We use $e < 87/32$ and $e - 2 < 23/32$, both immediate
from $e = 2.71828\ldots$ and certified in Arb by `largep_certify.py`, Part A.

- For $p\ge7$: $\frac p2 - 2e + 2 \ge \frac{11}{2} - 2e > \frac{11}{2} - \frac{87}{16} = \frac1{16}$.
- $\Gamma$ is increasing on $[2,\infty)$: $\psi = \Gamma'/\Gamma$ is strictly increasing
  ($\psi'(x) = \sum_{k\ge0}(x+k)^{-2} > 0$), and $\Gamma(1)=\Gamma(2)$ gives a zero of $\psi$ in $(1,2)$
  by Rolle's theorem. Hence for $p\ge7$, $\Gamma_p = p\,\Gamma(p) \ge p\,\Gamma(7) = 720p$.
- $\frac{e(p-2)}{2} \le \frac e2 p < \frac{87}{64}p$, and for $p\ge7$,
  $\frac{e(e-2)}{p+1} < \frac{87}{32}\cdot\frac{23}{32}\cdot\frac18 = \frac{2001}{8192}$.

Therefore, for $p \ge 7$,
$$K(p) > \frac{720p}{16} - \frac{87}{64}p - \frac{2001}{8192} = \frac{2793}{64}p - \frac{2001}{8192}
\ge \frac{2793\cdot7}{64} - \frac{2001}{8192} = \frac{2500527}{8192} > 0.$$

**Conclusion.** For $p\ge7$ and $\tau\in(0,1]$, $S(\tau) \ge \Delta(\tau) \ge \tau^2K(p) > 0$. By Step 3
this is $N(s)^{p/2}F(p) - M(s,p) \ge s^2(1-s)^{p-1}K(p)/e > 0$ for every $s\in(0,\tfrac12]$. Taking
logarithms, $G(s,p) - G_0(p) = \tfrac1p\log\big(M/(N^{p/2}F)\big) < 0$. The bound on $G_0 - G$ follows from
$-\log(1-y) \ge y$ for $y\in[0,1)$ applied to $y = (N^{p/2}F - M)/(N^{p/2}F)$. At $s=0$ both sides are
equal. $\blacksquare$

## Remarks

1. **Second-order structure.** In (3), $S(0) = 0$ and the $O(\tau)$ terms cancel exactly
   (Step 6). This is the identity $\partial_sG(0,p) = 0$ in the present variables. The only loss that
   is not lower order in $p$ is the Jensen chord (4), which costs $e\tau^2\Gamma_p$. It is beaten by the
   gain $A(\tau) - e^\tau \ge (\tfrac p2 - (e-2))\tau^2$ once $\tfrac p2 > 2e-2$. With exact Taylor
   coefficients the gain is $\tfrac{p-1}{2}\tau^2 + O(\tau^3)$, matching
   $\tfrac12\partial_s^2G(0,p) = \tfrac{p-1}{2p}(1/F - 1)\to -\tfrac12$.
2. **Threshold of this argument.** $K(p) > 0$ requires $\tfrac p2 - 2e + 2 > 0$, that is
   $p > 4e - 4 \approx 6.873$. The Arb audit confirms $K(6.8) < 0$ and $K(6.9) > 0$. Going below
   about $6.87$ requires replacing the chord bound (4) by a sharper bound on
   $\mathbb E(1-\tau+\tau R)^p = \tau^p e^{(1-\tau)/\tau}\Gamma(p+1,(1-\tau)/\tau)$. This is not attempted
   here, because the BULK and LOCAL regions already reach $P$.
3. **Partition.** Lemma LP covers $p \in [7,\infty)$. The production regions B4 and L3 end at $7$.
   Lemma LP holds for all $s\in[0,\tfrac12]$ including the corner $s\le\sigma_0$, so no local lemma is
   needed for $p\ge7$.

## Certificates and audit

- `code/largep/largep_certify.py` → `certificates/largep/largep_cert.json`, the sole shipped large-p certificate, using python-flint 0.9.0 Arb at 320 bits.
  - **Part A, used by the proof:** certifies $e<87/32$, $e-2<23/32$ and $e>5/2$ in Arb, and
    $\Gamma(7)=720$. It then redoes the Step 7 arithmetic in exact rationals: $\kappa = 1/16$,
    $720\kappa = 45$, linear coefficient $2793/64$, constant $2001/8192$, and
    $K\text{-lower}(7) = 2500527/8192$.
  - **Part B, audit only:** checks the evaluator against the nine reference values embedded in the script to $10^{-13}$.
    It then evaluates, at $p\in\{7, 15/2, 8, 10, 12, 20, 50, 100, 1000\}$ and $\tau = k/64$
    ($k=1,\dots,64$), the following rigorous balls:
    - $S$, both from the $\tau$-form (3) and directly from (1) at $s=\tau/(1+\tau)$, requiring
      agreement to relative $10^{-40}$;
    - the Jensen gap $j = (1-\tau)+\tau\Gamma_p-\mathbb E(1-\tau+\tau R)^p$, with the identity
      $S-\Delta = e\tau j$, $j>0$ for $\tau<1$ and $j=0$ at $\tau=1$;
    - the inequalities $\Delta > \tau^2K > 0$ and $G(s,p) < G_0(p)$.

    An undecided ball comparison raises, and the script exits with status 2.
  - **Part C, negative control:** the Part B audit is run at $p=6$ and $p=34/5$ and is required to
    fail there, since $K<0$. The script also certifies $K(34/5) < 0 < K(69/10)$, with
    $K(6.8)\approx-130.2$ and $K(6.9)\approx48.5$.
- `code/largep/largep_numeric_scan.py` is an optional mpmath, non-rigorous exploration script. Its historical JSON output is not shipped or required. Historical development notes report: for
  $p\in\{7,12,20,50,100,1000\}$ on $s = i/800$, the maximum of $G-G_0$ is about $-7\times10^{-7}$ at the
  smallest grid point $s = 1/800$. The ratio $(G-G_0)/s^2$ at $s=10^{-4}$ matches
  $\tfrac{p-1}{2p}(1/F-1)$ to about $10^{-4}$, for example $-0.45843$ against $-0.45833$ at $p=12$.

Reproduce from the repository root, choosing new output paths:
```
uv run --no-project --with python-flint==0.9.0 python3 code/largep/largep_certify.py --out NEW_largep_cert.json
uv run --no-project --with mpmath python3 code/largep/largep_numeric_scan.py --out NEW_largep_numeric_scan.json
```
