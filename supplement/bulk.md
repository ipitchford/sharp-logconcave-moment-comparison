# Appendix B. Certified inequalities in the bulk (regions B1–B4) and the crossing lemma

**Evidence status and paths.** The mathematical derivations below form this supplement. Numerical tables, timings and validation outcomes are historical development summaries, not a fresh independent production replay. Their rounded decimal margins are diagnostics, not outward-rounded certificates. Proof-critical scripts are in `code/bulk/` and outputs in `certificates/bulk/`; current hashes, exact exported bounds and acceptance results belong to those files and the release replay receipt. Auxiliary samples, formula checks, negative controls and wider-region runs do not replace exhaustive production enclosures. No unshipped development log is a dependency.

This appendix proves the four bulk inequalities B1–B4 and the crossing lemma X by
interval arithmetic. Every numerical step is an evaluation in Arb ball arithmetic
(python-flint 0.9.0) whose output ball provably contains the exact value. The analysis
below shows which exact statement each computed ball encloses.

## B.1 Notation and statements

Let $E, E', R$ be independent standard exponential random variables. For $s\in[0,\tfrac12]$ put

$$
Z_s=s(E-1)-(1-s)(E'-1),\qquad c=1-s,\qquad d=1-2s,\qquad x=\frac dc,\qquad z=\frac ds ,
$$

$$
N(s)=\mathbb E Z_s^2=s^2+(1-s)^2,\qquad M(s,p)=\mathbb E|Z_s|^p,\qquad
G(s,p)=\frac1p\log M(s,p)-\frac12\log N(s).
$$

Write $F(p)=\mathbb E|E-1|^p$, $G_0(p)=G(0,p)=\frac1p\log F(p)$ and
$G_h(p)=G(\tfrac12,p)=\frac1p\log\Gamma(p+1)-\frac12\log2$. Let $H=G_0-G_h$. The parameters are
$\eta=\sigma=\sigma_0=\delta=\tfrac1{20}$. The two rationals

$$
p^*_{lo}=\frac{1682457208}{10^9},\qquad p^*_{hi}=\frac{1682457209}{10^9}
$$

bracket the crossing point $p^*$ (Lemma X).

**Theorem B.** The following hold.

* **(X)** $H>0$ on $[1,p^*_{lo}]$, $H<0$ on $[p^*_{hi},\tfrac{39}{20}]$ and $H'<0$ on $[p^*_{lo},p^*_{hi}]$. Hence $H$ has exactly one zero $p^*$ in $[1,\tfrac{39}{20}]$, and $p^*_{lo}<p^*<p^*_{hi}$. In particular $G_h-G_0>0$ on $[p^*+\eta,\,2-\delta]$.
* **(B1)** $G(s,p)-G_h(p)>0$ for $s\in[0,\tfrac9{20}]$ and $p\in[1,\,p^*_{hi}-\eta]\supset[1,p^*-\eta]$.
* **(B2)** $G(s,p)-\min(G_0(p),G_h(p))>0$ for $s\in[\tfrac1{20},\tfrac9{20}]$ and $p\in[p^*_{lo}-\eta,\,p^*_{hi}+\eta]\supset[p^*-\eta,p^*+\eta]$.
* **(B3)** $G(s,p)-G_0(p)>0$ for $s\in[\tfrac1{20},\tfrac12]$ and $p\in[p^*_{lo}+\eta,\,\tfrac{39}{20}]\supset[p^*+\eta,2-\delta]$.
* **(B4)** $G_0(p)-G(s,p)>0$ for $s\in[\tfrac1{20},\tfrac12]$ and $p\in[\tfrac{41}{20},\,7]$.

B4 stops at $p=7$ because the analytic large-$p$ lemma covers $p\ge7$; the original partition had $P=12$. Historical development records also report region B4ext on $p\in[7,12]$ (§B.7), duplicating part of that lemma. Its optional certificate is not shipped or required.

The regions are slightly enlarged so that every endpoint is an explicit rational. The enlargement uses only $p^*\in(p^*_{lo},p^*_{hi})$ from (X).

## B.2 Closed forms

**Lemma B.1 (conditioning representation).** For $s\in(0,\tfrac12]$ and $p>-1$,

$$
M(s,p)=s\,J_p(s,d)+c^{p+1}\varphi(x,p),
$$

where, for $q\ge0$, $s>0$ and $d\ge0$,

$$
J_q(s,d)=\mathbb E(d+sR)^q=s^q e^{z}\,\Gamma(q+1,z)\quad(z=d/s),\qquad J_q(s,0)=s^q\Gamma(q+1),
$$

$$
\varphi(x,p)=\mathbb E|x-R|^p=e^{-x}\bigl(\Gamma(p+1)+I(x,p)\bigr),\qquad
I(x,p)=\int_0^x u^pe^u\,du=\frac{x^{p+1}}{p+1}\,{}_1F_1(p+1;p+2;x).
$$

For $s=0$ we have $M(0,p)=\varphi(1,p)=F(p)$. Here $\Gamma(a,z)=\int_z^\infty t^{a-1}e^{-t}dt$.

*Proof.* $Y=sE-cE'$ has density $f(y)=e^{-y/s}$ for $y>0$ and $f(y)=e^{y/c}$ for $y<0$; the normalisation is $s+c=1$. Hence

$$
M=\int_0^\infty (y+d)^p e^{-y/s}dy+\int_0^\infty|d-y|^pe^{-y/c}dy=s\,\mathbb E(d+sR)^p+c\,\mathbb E|d-cR|^p .
$$

Substituting $t=z+r$ gives $\mathbb E(d+sR)^q=s^q\int_0^\infty(z+r)^qe^{-r}dr=s^qe^z\Gamma(q+1,z)$. Next, $c\,\mathbb E|d-cR|^p=c^{p+1}\mathbb E|x-R|^p$. Splitting at $r=x$,

$$
\int_0^x(x-r)^pe^{-r}dr=e^{-x}\int_0^xu^pe^udu ,\qquad \int_x^\infty(r-x)^pe^{-r}dr=e^{-x}\Gamma(p+1).
$$

Termwise integration of $e^u=\sum u^k/k!$ gives $I(x,p)=\sum_{k\ge0}\frac{x^{p+k+1}}{k!(p+k+1)}$. Since $(p+1)_k/(p+2)_k=(p+1)/(p+k+1)$, this series equals the ${}_1F_1$ expression. $\square$

At $s=\tfrac12$ we have $d=x=0$ and $M=2^{-p}\Gamma(p+1)$, so $G(\tfrac12,p)=G_h(p)$.

**Lemma B.2 ($s$-derivative).** For $p\ge1$, the function $s\mapsto M(s,p)$ is $C^1$ on $[0,\tfrac12]$, and

$$
\partial_sM=(p+1)J_p-p\,J_{p-1}-c^{p-1}\bigl((p+1)c-1\bigr)\varphi(x,p)-\frac{d^p}{c},
$$

with $J_q=J_q(s,1-2s)$. At $s=0$, read $J_q=1$.

*Proof.* Write $T_A=sJ_p(s,1-2s)$ and $T_B=c^{p+1}\varphi(x,p)$. Differentiation under $\mathbb E$ is legitimate by dominated convergence for $q=p-1\ge0$. It gives $\partial_sJ_q=q\,\mathbb E[(d+sR)^{q-1}R]$ and $\partial_dJ_q=q\,\mathbb E[(d+sR)^{q-1}]$. Hence

$$
\frac{d}{ds}T_A=J_p+s\,p\,\mathbb E\bigl[(d+sR)^{p-1}(R-2)\bigr]=J_p+p\,\mathbb E\bigl[(d+sR)^{p-1}\bigl((d+sR)-(d+2s)\bigr)\bigr]=(p+1)J_p-pJ_{p-1},
$$

using $d+2s=1$. From Lemma B.1, $\partial_x\varphi=x^p-\varphi$, and $dx/ds=-1/c^2$. Therefore

$$
\frac{d}{ds}T_B=-(p+1)c^p\varphi-c^{p-1}(x^p-\varphi)=-c^{p-1}\bigl((p+1)c-1\bigr)\varphi-\frac{d^p}{c},
$$

because $c^{p-1}x^p=d^p/c$. As $s\downarrow0$ we have $J_q(s,1-2s)\to1$ and $T_A/s\to1$, so the formula holds at $s=0$ as a one-sided derivative and is continuous there. $\square$

As a check, at $s=0$ the formula gives $\partial_sM=-pF(p)$, i.e. $\partial_sG(0,p)=0$. At $s=\tfrac12$ it gives $0$.

**Lemma B.3 ($p$-derivative).** For $s\in(0,\tfrac12]$ and $p>-1$,

$$
\partial_pM=\log s\cdot sJ_p+s^{p+1}e^{z}K(p+1,z)+\log c\cdot c^{p+1}\varphi+c^{p+1}e^{-x}\bigl(\Gamma(p+1)\psi(p+1)+I_p(x,p)\bigr),
$$

where

$$
I_p(x,p)=\partial_pI=\log x\cdot I(x,p)-\frac{x^{p+1}}{(p+1)^2}\,{}_2F_2(p+1,p+1;p+2,p+2;x),
$$

$$
K(a,z)=\partial_a\Gamma(a,z)=\Gamma(a)\psi(a)-\Bigl[\log z\cdot\gamma(a,z)-\frac{z^a}{a^2}\,{}_2F_2(a,a;a+1,a+1;-z)\Bigr],
$$

with $\gamma(a,z)=\Gamma(a)-\Gamma(a,z)$. When $z=0$ (that is, $s=\tfrac12$), read $K(a,0)=\Gamma(a)\psi(a)$ and $I_p=0$. For $s=0$, only the $c$-terms with $c=1,x=1$ remain. In every case $\partial_pM=\mathbb E\bigl[|Z_s|^p\log|Z_s|\bigr]$.

*Proof.* Differentiate the closed form of Lemma B.1 in $p$; note that $z$ and $x$ do not depend on $p$. Termwise differentiation of the series for $I$ gives

$$
\partial_p I=\sum_k\frac{x^{p+k+1}}{k!}\Bigl[\frac{\log x}{p+k+1}-\frac1{(p+k+1)^2}\Bigr].
$$

Because $(p+1)_k^2/(p+2)_k^2=(p+1)^2/(p+k+1)^2$, the second part is the ${}_2F_2$ term. The same computation for

$$
\gamma(a,z)=\sum_k(-1)^k\frac{z^{a+k}}{k!(a+k)}
$$

gives $\partial_a\gamma$. The identity with $\mathbb E|Z|^p\log|Z|$ follows by differentiating $\mathbb E|Z_s|^p$ under the expectation, which dominated convergence justifies because $|Z_s|$ has a bounded density and exponential tails. $\square$

The formula for $K$ cancels about $2.9z$ bits when $z$ is large. The certificate therefore evaluates $\partial_pM$ only at exact rational points, at working precision $192+3z$ bits. The inputs are converted to balls after the precision is raised. A result with relative radius above $10^{-20}$ aborts the run.

The derivatives of the references are

$$
G_0'(p)=\frac{F'(p)}{pF(p)}-\frac{\log F(p)}{p^2},\qquad F'(p)=e^{-1}\bigl(\Gamma(p+1)\psi(p+1)+I_p(1,p)\bigr),
$$

$$
G_h'(p)=\frac{\psi(p+1)}p-\frac{\log\Gamma(p+1)}{p^2}.
$$

## B.3 Monotonicity facts

**Lemma B.4 (log-convexity).** Let $Y\ge0$ be a random variable with $\mathbb EY^{q}<\infty$ for $q$ in an open interval $\mathcal I$, and let $m(p)=\mathbb EY^p$. Then $m$ is convex on $\mathcal I$, and $m'(p)=\mathbb E[Y^p\log Y]$ is nondecreasing there.

*Proof.* On compact subsets of $\mathcal I$, both $Y^p|\log Y|$ and $Y^p\log^2Y$ are dominated by integrable functions of the form $C(Y^{q_0}+Y^{q_1})$ with $q_0<q_1$ in $\mathcal I$. Hence $m''(p)=\mathbb E[Y^p\log^2Y]\ge0$. $\square$

Applied to $Y=|Z_s|$, $Y=|E-1|$ and $Y=d+sR$, this gives:

* $\partial_pM(s,\cdot)$ is nondecreasing;
* $F'$ is nondecreasing;
* $q\mapsto J_q(s,d)$ is convex.

Consequently, for every interval $[p_0,p_1]$ and every $[q_0,q_1]$,

$$
\partial_pM(s,[p_0,p_1])\subset[\partial_pM(s,p_0),\,\partial_pM(s,p_1)],\qquad F'([p_0,p_1])\subset[F'(p_0),F'(p_1)],
$$

$$
0<J_q(s,d)\le\max\bigl(J_{q_0}(s,d),J_{q_1}(s,d)\bigr)\quad(q\in[q_0,q_1]).
$$

**Lemma B.5 (monotone enclosures).** Let $0\le s_0\le s\le s_1\le\tfrac12$ and $1\le p_0\le p\le p_1$.

1. **(in $s$ and $d$)** For $q\ge0$, $J_q(s,d)=\mathbb E(d+sR)^q$ is nondecreasing in $s\ge0$ and in $d\ge0$. Since $d=1-2s\in[1-2s_1,1-2s_0]$,
   $$J_q(s_0,1-2s_1)\le J_q(s,1-2s)\le J_q(s_1,1-2s_0),\qquad s_0J_p(s_0,1-2s_1)\le T_A(s)\le s_1J_p(s_1,1-2s_0),$$
   with $J_q(0,d)=d^q$.
2. **(in $q$)** For $s,d>0$ and $z=d/s$,
   $$J_q(s,d)=d^q\,\Phi_q(z),\qquad \Phi_q(z)=\mathbb E(1+R/z)^q=z^{-q}e^z\Gamma(q+1,z),$$
   and $\Phi_q(z)$ is nondecreasing in $q$ because $1+R/z\ge1$. Hence for $q\in[q_0,q_1]$,
   $$J_q\in d^{[q_0,q_1]}\cdot[\Phi_{q_0}(z),\Phi_{q_1}(z)].$$
3. **($x$)** $x(s)=(1-2s)/(1-s)$ is decreasing, because $x'(s)=-1/(1-s)^2$. Hence $x\in[x(s_1),x(s_0)]\subset[0,1]$.
4. **($I$)** $I(x,p)=\int_0^xu^pe^udu$ is nondecreasing in $x\ge0$. For $x\le1$ it is nonincreasing in $p$, because $u^p$ is nonincreasing in $p$ for $u\in[0,1]$. So $I(x,p)\in[I(x(s_1),p_1),\,I(x(s_0),p_0)]$ on the box.
5. **(elementary)** On $[p_0,p_1]\subset[1,\infty)$, $\Gamma(p+1)$, $\psi(p+1)$ and $\log\Gamma(p+1)$ are nondecreasing; $\Gamma$ increases on $[1.4617,\infty)$, and $\psi$ is increasing. Also $e^{-x}$ is decreasing in $x$, and $c\mapsto c^{p\pm1}$ is monotone on $c\in[1-s_1,1-s_0]$. Finally $(2s-1)/N(s)$ is nondecreasing on $[0,\tfrac12]$ (its derivative is $4s(1-s)/N(s)^2$), and $\partial_s\bigl(\tfrac12\log N\bigr)=(2s-1)/N$.

*Proof.* (1) The integrand $(d+sr)^q$ is nondecreasing in $d$ and in $s$ for $r\ge0$. (2) $J_q=\mathbb E\bigl(d(1+R/z)\bigr)^q$, and the incomplete-gamma form follows as in Lemma B.1. (3)–(5) follow by direct computation. $\square$

**Enclosures on a box $B=[s_0,s_1]\times[p_0,p_1]$.** Every factor of Lemma B.1 and Lemma B.2 is enclosed from evaluations at exact endpoints:

* **$J_p$, $J_{p-1}$ and $T_A$.** Item 1 reduces them to the two corners $(s_0,1-2s_1)$ and $(s_1,1-2s_0)$. At each corner, the dependence on $q\in[p_0,p_1]$ (respectively $[p_0-1,p_1-1]$) is enclosed by intersecting three valid enclosures:
  * the upper bound of Lemma B.4;
  * item 2;
  * when $z<1$, the plain Arb evaluation of $s^qe^z\Gamma(q+1,z)$ on the ball $q\in[q_0,q_1]$.
* **$\Gamma(p+1)+I(x,p)$.** By items 4–5.
* **$c^{p+1}$, $c^{p-1}$, $e^{-x}$, $d^p$.** By items 3 and 5, with the exponent ball $P\supseteq[p_0,p_1]$. These are elementary functions, which Arb evaluates on balls without significant over-estimation.

This interval-endpoint construction matters for efficiency, not for rigour. Arb's $\Gamma(a,z)$ with a wide ball $a$ is rigorous but can over-estimate enormously. For example, $a\in[1.32,1.63]$ and $z=44$ return $[\pm189]$ for a quantity of size $10^{-18}$. Such balls would stall the subdivision.

Inclusion isotonicity of interval arithmetic gives enclosures $[M](B)\ni M(s,p)$ and $[M_s](B)\ni\partial_sM(s,p)$ for all $(s,p)\in B$. When $[M](B)$ has a positive lower bound,

$$
[G_s](B)=\frac{[M_s](B)}{P\,[M](B)}-\Bigl[\frac{2s_0-1}{N(s_0)},\frac{2s_1-1}{N(s_1)}\Bigr]\ \ni\ \partial_sG(s,p)\qquad\text{for all }(s,p)\in B.
$$

In the same way, $M(s_c,p)$ for $p\in[p_0,p_1]$ is enclosed using items 2, 4 and 5, and $F(p)=\varphi(1,p)$ using items 4 and 5.

## B.4 The box enclosure

**Proposition B.6.** Let $R\in\{G_0,G_h\}$ and $D(s,p)=G(s,p)-R(p)$. Let $B=[s_0,s_1]\times[p_0,p_1]$ be a box with centre $(s_c,p_c)$ and half-widths $h_s,h_p$. Suppose $[D_p]\supseteq\{\partial_pD(s_c,p):p\in[p_0,p_1]\}$ and $[G_s]\supseteq\{\partial_sG(s,p):(s,p)\in B\}$. Then for all $(s,p)\in B$,

$$
D(s,p)\in D(s_c,p_c)+[D_p]\cdot[-h_p,h_p]+[G_s]\cdot[-h_s,h_s].
$$

*Proof.* Since $R$ does not depend on $s$,

$$
D(s,p)-D(s_c,p_c)=\bigl[D(s_c,p)-D(s_c,p_c)\bigr]+\bigl[G(s,p)-G(s_c,p)\bigr].
$$

Apply the mean value theorem in $p$ to the first bracket and in $s$ to the second. $M>0$ and $M$ is $C^1$ (Lemmas B.2 and B.3), so $G$ is $C^1$ on $B$. $\square$

The certificate obtains the three ingredients as follows.

* **Centre value.** $D(s_c,p_c)$ is computed from Lemma B.1 at the exact rational point. Its radius is about $10^{-35}$.
* **$[D_p]$ ("mv-p" mode, used when $s_c\ge\tfrac1{64}$).**
  $$\partial_pG(s_c,p)=\frac{\partial_pM}{pM}-\frac{\log M}{p^2}.$$
  The numerator $\partial_pM$ is enclosed by $[\partial_pM(s_c,p_0),\partial_pM(s_c,p_1)]$ (Lemma B.4), evaluated by Lemma B.3 at exact rational points with raised precision. $M(s_c,[p_0,p_1])$ is enclosed as in §B.3. For the reference, $G_0'$ uses $[F'(p_0),F'(p_1)]$ and $F([p_0,p_1])$; $G_h'=\psi(p+1)/p-\log\Gamma(p+1)/p^2$ uses item 5.
* **"ball-p" mode ($s_c<\tfrac1{64}$, region B1 only).** The term $D(s_c,p_c)+[D_p][-h_p,h_p]$ is replaced by the interval enclosure of $D(s_c,[p_0,p_1])$. That enclosure is $\log M(s_c,[p_0,p_1])/P-\tfrac12\log N(s_c)-R([p_0,p_1])$, and it needs no derivative in $p$.
* **$[G_s]$.** As in §B.3.

A box is *certified* when the lower endpoint of the resulting ball is $>0$; for B4 the upper endpoint must be $<0$. Arb guarantees that every evaluated ball contains the exact value for all points of the input balls. An undecided comparison counts as "not certified", and the box is then bisected.

Proposition B.6 is a mixed first-order (centred) form, which explains its efficiency. The over-estimation of $[D_p]$ and $[G_s]$ is $O(h)$, so the enclosure width is

$$
2|\partial_pD|h_p+2|\partial_sG|h_s+O(h^2).
$$

The first two terms are the true variation of $D$ on the box. Near $p=2$ we have $D(s,p)\approx(2-p)\,\kappa(p)s^2$, with both $D$ and its gradient of size $|2-p|$. A plain ball evaluation in $p$ would instead need $h_p\lesssim D$, i.e. $h_p\approx10^{-5}$ at the corner $(s,p)=(\tfrac1{20},\tfrac{39}{20})$.

**Regions with two references (B2).** $G>\min(G_0,G_h)$ holds at a point if $G>G_0$ or $G>G_h$ there. A box of B2 is certified when either $D=G-G_0$ or $D=G-G_h$ has a positive certified lower bound on the whole box.

## B.5 Branch and bound

Each region starts as one box with exact rational endpoints. Boxes are processed depth first. An uncertified box is bisected in $s$ if $|[G_s]|h_s>|[D_p]|h_p$, and in $p$ otherwise. When an enclosure is not finite, the relatively wider side is bisected. The ratio of the normalised side lengths is capped at $2^{12}$.

A box narrower than $2^{-34}$ in both directions counts as a failure. So does a box whose centre value has the certified wrong sign, which is a counterexample. On either kind of failure the run stops and the script exits with a non-zero status. Every certified leaf is hashed (SHA-256 of the list of rational leaf boxes in processing order), and the hash is recorded in the certificate.

**Lemma X (crossing).** This lemma is one-dimensional. $H$ is enclosed on $[p_0,p_1]$ by $H(p_c)+[H'](P)[-h,h]$, where

$$
H'=G_0'-G_h' ,
$$

with $F'(P)\subset[F'(p_0),F'(p_1)]$ (Lemma B.4) and all other quantities evaluated on the ball $P$. X3 is the direct enclosure of $H'$ on $[p^*_{lo},p^*_{hi}]$. Uniqueness of the zero follows from X1–X3: $H$ is positive on $[1,p^*_{lo}]$, strictly decreasing on $[p^*_{lo},p^*_{hi}]$, and negative on $[p^*_{hi},\tfrac{39}{20}]$, so by continuity it has exactly one zero there.

## B.6 Validation of the evaluator

Before any certification, `bulk_cert.py` runs the following self-checks and aborts on any mismatch.

* It evaluates the embedded reference values; the tolerance is $10^{-14}$. All seven $G$ values agree with the 40-digit mpmath values to every printed digit, and so do $G_0(1.3)$ and $G_h(1.3)$. See the table below.
* It compares $\partial_sM$ and $\partial_pM$ with central differences of the closed form for $M$, at 256 bits with $h=10^{-12}$, at four points; the tolerance is $10^{-15}$.
* It runs the containment test of §B.7.

Separately, `check_formulas.py` compares $M$, $\partial_sM$ and $\partial_pM$ from the certified code with independent mpmath quadrature of the defining integrals and with `mp.diff` of the quadrature. It also checks the conditioning form of $M_s$,

$$
M_s=p\bigl(s\,\mathbb E[(d+sR)^{p-1}(R-1)]+c\,\mathbb E[|d-cR|^{p-2}(d-cR)(R-1)]\bigr).
$$

At 15 points (7 fixed, 8 random) the worst relative discrepancy is $2.9\times10^{-25}$; this historical diagnostic output is not a shipped dependency.

| $(s,p)$ | $G$ (Arb, closed form) | Reference (mpmath) |
|---|---|---|
| (0.1, 1.3) | −0.213049905296182717 | −0.213049905296183 |
| (0.3, 1.9) | −0.0294149150315885971 | −0.0294149150315886 |
| (0.25, 2.5) | 0.143192764694942902 | 0.143192764694943 |
| (0.4, 5) | 0.641602497665675327 | 0.641602497665675 |
| (0.05, 12) | 1.58099627165212506 | 1.58099627165213 |
| (0.5, 1.2) | −0.265784034621107008 | −0.265784034621107 |
| (0, 3.3) | 0.375483768956575059 | 0.375483768956575 |
| $G_0(1.3)$ | −0.213229500079512965 | −0.213229500079513 |
| $G_h(1.3)$ | −0.227966317234102977 | −0.227966317234103 |

**Negative controls.** The certifier must fail on false statements, and it does (`bulk_cert.py CTRL_*`). It was run on four:

* B3's inequality on $p\in[\tfrac{41}{20},\tfrac{21}{10}]$;
* B4's inequality on $p\in[\tfrac{19}{10},\tfrac{39}{20}]$;
* B1's inequality on $p\in[1.7,1.8]$;
* B3's inequality on $p\in[1.9,2.01]$, which is true below $p=2$ and false above.

Each run exits with status 1 and reports a certified counterexample at a box centre. In the last case the counterexample appears after 2506 boxes, in the box $[\tfrac1{20},\tfrac{137}{2560}]\times[\tfrac{1597}{800},\tfrac{201}{100}]$, whose centre has $p\approx2.0031$. Historical control outputs are not shipped or required.

## B.7 Results

All runs used python-flint 0.9.0 at 128-bit working precision, with $\partial_pM$ evaluated at $192+3z$ bits. Historical runs used a single thread. "Min. certified bound" is the smallest lower bound of $\pm D$ over the certified leaves; it is a rigorous lower bound for the margin on the whole region. "Approx. true min." is the minimum of $\pm D$ at the region's worst corner, from point evaluation; it is not rigorous.

| Region | $s$ | $p$ | Claim | Boxes evaluated / leaves | Failures | Min. certified bound | Approx. true min. (location) | Runtime |
|---|---|---|---|---|---|---|---|---|
| B1 | $[0,\tfrac9{20}]$ | $[1,\ 1.632457209]$ | $G-G_h>0$ | 1653 / 827 | 0 | $3.24\times10^{-6}$ | $3.5\times10^{-4}$ $(0.45,\,p^*-\eta)$ | 1.0 s |
| B2 | $[\tfrac1{20},\tfrac9{20}]$ | $[1.632457208,\ 1.732457209]$ | $G-\min(G_0,G_h)>0$ | 1127 / 564 (455 via $G_0$, 109 via $G_h$) | 0 | $4.03\times10^{-8}$ | $7.8\times10^{-5}$ $(0.05,\,p^*+\eta)$ | 1.0 s |
| B3 | $[\tfrac1{20},\tfrac12]$ | $[1.732457208,\ \tfrac{39}{20}]$ | $G-G_0>0$ | 4413 / 2207 | 0 | $1.43\times10^{-8}$ | $1.8\times10^{-5}$ $(0.05,\,1.95)$ | 3.8 s |
| B4 | $[\tfrac1{20},\tfrac12]$ | $[\tfrac{41}{20},\ 7]$ | $G_0-G>0$ | 8635 / 4318 | 0 | $8.02\times10^{-8}$ | $2.0\times10^{-5}$ $(0.05,\,2.05)$ | 7.6 s |
| B4ext (optional) | $[\tfrac1{20},\tfrac12]$ | $[7,\ 12]$ | $G_0-G>0$ | 4107 / 2054 | 0 | $4.70\times10^{-5}$ | $\approx1.2\times10^{-3}$ $(0.05,\,7)$ | 3.0 s |

| Lemma X part | Interval | Boxes | Failures | Certified bound |
|---|---|---|---|---|
| X1: $H>0$ | $[1,\ 1.682457208]$ | 73 | 0 | $H\ge3.5\times10^{-12}$ |
| X2: $H<0$ | $[1.682457209,\ \tfrac{39}{20}]$ | 65 | 0 | $-H\ge6.4\times10^{-12}$ |
| X3: $H'<0$ | $[1.682457208,\ 1.682457209]$ | 1 | 0 | $H'\le-0.016038$; enclosure $[-0.01603876\pm9\times10^{-9}]$ |

The certified bracket is therefore $p^*\in(1.682457208,\,1.682457209)$. Bisection of the certified sign of $H$ (informational) gives $p^*=1.68245720845516070\ldots$

**Historical higher-precision diagnostic (outputs not shipped).** Development records report successful B1–B4 runs at 192-bit working precision. The box counts are 1655, 1127, 4415 and 8635. The B2 and B4 leaf hashes were identical to those of the historical 128-bit runs. These optional outputs are not current release bindings or proof dependencies.

**Containment test.** Each certificate run first evaluates $M$, $\partial_sM$, $\partial_sG$, $M(s_c,\cdot)$ and $G-R$ at 96 exact random interior points of 24 random boxes. It checks that every point value lies in the corresponding box enclosure.

**Historical partition trials (not production evidence).** Development notes report the following; no exploration directory is shipped or required:

* the B3/B4 inequalities with $\delta=\tfrac1{100}$ (8595 and 12179 boxes);
* the B3/B4 inequalities with $\delta=\tfrac1{1000}$ (14885 and 18773 boxes);
* B2–B4 with $\sigma_0=\tfrac1{40}$;
* B1–B2 with $\sigma=\tfrac1{100}$.

Each run takes under 20 s. With $\sigma_0=\tfrac1{100}$ and $\delta=\tfrac1{100}$ together, a cap of $4\times10^5$ boxes was not enough. These trials do not change the production partition.

## B.8 Reproduction

Run from the repository root with fresh output names. The shipped bulk inventory is exactly `certificates/bulk/cert_B1.json`, `cert_B2.json`, `cert_B3.json`, `cert_B4.json` and `cert_X.json`. B4ext and 192-bit outputs are omitted historical diagnostics; their optional commands generate new results rather than replay a shipped certificate.

```
uv run --no-project --with python-flint==0.9.0 python3 code/bulk/bulk_cert.py B1 --out NEW_cert_B1.json   # likewise B2, B3, B4
uv run --no-project --with python-flint==0.9.0 python3 code/bulk/bulk_cert.py B3 --prec 192           # optional new diagnostic
uv run --no-project --with python-flint==0.9.0 python3 code/bulk/bulk_cert.py CTRL_B3_overshoot      # must exit 1
uv run --no-project --with python-flint==0.9.0 python3 code/bulk/crossing_cert.py --out NEW_cert_X.json
uv run --no-project --with python-flint==0.9.0 --with mpmath python3 code/bulk/check_formulas.py
```

Each certificate records:

* the SHA-256 of the script that produced it;
* the exact rational region and all parameters;
* the python-flint version and the precision;
* the box counts;
* the minimum certified margin and the box where it is attained;
* the leaf hash and the self-test output.

Current source hashes are recorded in production certificates under `certificates/bulk/`; historical hashes are not current release bindings.
