# Introduction

For a random variable $X$ write $\|X\|_p=(\mathbb E|X|^p)^{1/p}$. Let
$\mathcal L$ be the class of nondegenerate mean-zero log-concave random
variables, that is, those with a density $e^{-V}$ for a proper convex
$V:\mathbb R\to(-\infty,+\infty]$ (so supports may be half-lines or
intervals; $E-1\in\mathcal L$). Melbourne, Roysdon, Tang and Tkocz
[@MRTT] determined the sharp constants comparing $\|X\|_p$ with
$\|X\|_1$ on $\mathcal L$, finding a phase transition in the extremiser.
In their concluding remarks [@MRTT Sec. 3.1] they ask for the sharp
constants in the variance-normalised problem,
$$m_p=\inf_{X\in\mathcal L}\frac{\|X\|_p}{\|X\|_2}\quad(1<p<2),\qquad
M_p=\sup_{X\in\mathcal L}\frac{\|X\|_p}{\|X\|_2}\quad(p>2),$$ and, based
on numerical experiments, conjecture that the extremiser is the
symmetric exponential for $1<p<p^*\approx1.68$ and the one-sided
exponential for $p>p^*$. The phase transition and the value
$p^*\approx1.68$ were thus found numerically by Melbourne, Roysdon, Tang
and Tkocz; our contribution is a proof. Eitan [@Eitan Theorem 1.4]
(Theorem 4 in the arXiv version) proved the case of even integers $p>2$.

At $p=0$ we use the limiting convention
$\|X\|_0=\exp(\mathbb E\log|X|)$. The logarithm is integrable for the
laws considered here: the density is bounded near zero and has
exponential tails. The new results below concern $p>1$.

Let $E,E'$ be independent standard exponential random variables,
$L=E-E'$ (symmetric exponential) and $F(p)=\mathbb E|E-1|^p$. Then
$\|L\|_p/\|L\|_2=\Gamma(p+1)^{1/p}/\sqrt2$ and
$\|E-1\|_p/\|E-1\|_2=F(p)^{1/p}$.

::: {#thm:main .theorem}
**Theorem 1**. *There is a unique $p^*\in(1,2)$ with
$\Gamma(p^*+1)^{1/p^*}/\sqrt2=F(p^*)^{1/p^*}$; numerically
$p^*=1.6824572084\ldots$. Moreover:*

- *for $1<p\le p^*$, $m_p=\Gamma(p+1)^{1/p}/\sqrt2$, attained by the
  symmetric exponential;*

- *for $p^*\le p<2$, $m_p=F(p)^{1/p}$, attained by the one-sided
  exponential $E-1$;*

- *for $p>2$, $M_p=F(p)^{1/p}$, attained by $E-1$.*
:::

With [@MRTT Theorem 1], which treats $-1<p\le1$ (symmetric exponential
extremal), Theorem [1](#thm:main){reference-type="ref"
reference="thm:main"} gives the sharp constants comparing $L_p$ with
$L_2$ on $\mathcal L$ for all $p>-1$; for even $p>2$ it recovers [@Eitan
Theorem 1.4]. We also classify the cases of equality. Here
$X\overset{d}{=}Y$ means equality in law.

::: {#thm:U .theorem}
**Theorem 2** (Equality cases). *Let $X\in\mathcal L$, and let $p^*$,
$m_p$, $M_p$ be as in Theorem [1](#thm:main){reference-type="ref"
reference="thm:main"}.*

1.  *If $1<p<p^*$, then $\|X\|_p=m_p\|X\|_2$ if and only if
    $X\overset{d}{=}\lambda(E-E')$ for some $\lambda>0$.*

2.  *If $p=p^*$, then $\|X\|_p=m_p\|X\|_2$ if and only if
    $X\overset{d}{=}\lambda(E-E')$ or $X\overset{d}{=}\pm\lambda(E-1)$
    for some $\lambda>0$.*

3.  *If $p^*<p<2$, then $\|X\|_p=m_p\|X\|_2$ if and only if
    $X\overset{d}{=}\pm\lambda(E-1)$ for some $\lambda>0$.*

4.  *If $p>2$, then $\|X\|_p=M_p\|X\|_2$ if and only if
    $X\overset{d}{=}\pm\lambda(E-1)$ for some $\lambda>0$.*
:::

#### Contribution and proof architecture.

The new optimisation and complete classification in
Theorems [1](#thm:main){reference-type="ref"
reference="thm:main"}--[2](#thm:U){reference-type="ref"
reference="thm:U"} concern the remaining real exponents. Eitan's
reduction already covers real $1<a<b$; neither that scope nor the
crossing method is claimed as new.
Section [2](#sec:selfred){reference-type="ref" reference="sec:selfred"}
supplies a self-contained, equality-sensitive proof of the reduction,
following Eitan [@Eitan Sections 2--3] and the crossing methods of
Eskenazis, Nayar and Tkocz [@ENT].
Sections [4](#sec:largep){reference-type="ref"
reference="sec:largep"}--[8](#sec:assembly){reference-type="ref"
reference="sec:assembly"} establish the reduced optimisation, and
Section [9](#sec:equality){reference-type="ref"
reference="sec:equality"} extracts the strictness needed for the
equality classification. All arguments and proof-critical computational
material are included in this paper and its accompanying public package;
no earlier draft is required.

#### Related work.

Murawski [@Murawski] proved the mean-zero bound
$\|X\|_p\le\frac pq\|X\|_q$ ($p\ge q\ge1$) and sharp shifted-exponential
results for general (not necessarily centred) log-concave variables.
Eitan conjectured [@Eitan Conjecture 16 in arXiv v1] that for even $p<q$
the extremal $Z_s$ is $s\in\{0,1\}$;
Theorem [1](#thm:main){reference-type="ref" reference="thm:main"} shows
that this endpoint behaviour fails for real exponents, since for $(p,2)$
with $p<p^*$ the extremal parameter is $s=\frac12$.
Theorem [3](#thm:R2){reference-type="ref" reference="thm:R2"} reduces
that conjecture, with its equality case, to an inequality within the
family $Z_s$; we do not address it. Since $Z_s/\|Z_s\|_2$ is a two-term
centred sum of exponentials, the problem also belongs to the study of
sharp moment bounds for sums of exponential random variables
[@BTT; @BP].

#### Method.

By Theorem [3](#thm:R2){reference-type="ref" reference="thm:R2"}
(Section [2](#sec:selfred){reference-type="ref"
reference="sec:selfred"}) the extremal ratio is attained within the
family $Z_s=s(E-1)-(1-s)(E'-1)$, $0\le s\le1$, and only there, so the
theorem is equivalent to a two-parameter inequality for
$G(s,p)=\frac1p\log\mathbb E|Z_s|^p-\frac12\log\mathbb EZ_s^2$
(Section [3](#sec:reduction){reference-type="ref"
reference="sec:reduction"}). The inequality is degenerate along $s=0$
and $s=1/2$ (where $\partial_sG$ vanishes), along $p=2$ (where
$G\equiv0$), and at the crossing $p=p^*$. We prove it analytically for
$p\ge7$ (Section [4](#sec:largep){reference-type="ref"
reference="sec:largep"}). For $1<p<7$ we partition
$[0,\frac12]\times(1,7]$ into regions with explicit margins, certified
by ball arithmetic, and regions around the degenerate lines, where
second-order lemmas reduce the claim to certified sign conditions on
derivatives (Sections [5](#sec:local){reference-type="ref"
reference="sec:local"}--[7](#sec:bulk){reference-type="ref"
reference="sec:bulk"}). All the certified inequalities are strict, and
this gives the equality cases
(Section [9](#sec:equality){reference-type="ref"
reference="sec:equality"}).

#### A concrete sharp bound.

At $p=4$, direct integration gives
$F(4)=\mathbb E(E-1)^4=24-24+12-4+1=9$. Thus every $X\in\mathcal L$
satisfies $\mathbb EX^4\le9(\mathbb EX^2)^2$, with equality exactly for
a scaled or reflected centred exponential. This even-integer case was
already known from Eitan; it illustrates the normalisation rather than a
new fourth-moment inequality. In contrast, at $p=3/2<p^*$ the sharp
lower ratio is $\Gamma(5/2)^{2/3}/\sqrt2$ and only the Laplace family
attains it. At $p=2$ the ratio equals one for every admissible law.
Section [10](#sec:stability){reference-type="ref"
reference="sec:stability"} also shows that a nearly extremal ratio at a
fixed $p>1$, $p\ne2$, forces closeness in total variation and every
finite Wasserstein distance to the corresponding variance-normalised
equality set. This qualitative consequence has no explicit
deficit-to-distance rate.

# A self-contained reduction (after Eitan) {#sec:selfred}

Throughout, $E,E'$ are independent standard exponential random variables
and $$Z_s=s(E-1)-(1-s)(E'-1),\qquad 0\le s\le1.$$ Thus
$Z_{1/2}=\frac12(E-E')$, $Z_1=E-1$, $Z_0=1-E'\overset{d}{=}-(E-1)$, and
$Z_{1-s}\overset{d}{=}-Z_s$ (exchange $E$ and $E'$).

#### Relation to Eitan's work.

The two theorems below give self-contained versions of Theorems 12
and 14 of Eitan's arXiv:2110.14382v1 [@Eitan]; Theorem 12 corresponds to
Theorem 3.1 of the published article. The real-exponent reduction is
already established there. Although its arXiv Theorem 12 does not state
equality explicitly, Theorem 14 does, and the proof of Theorem 4 uses an
equality characterisation for the reduction. Thus we do not claim
priority for either the reduction or equality-sensitive reasoning. Our
detailed comparison uses the arXiv text; the published full text was not
inspected. The proofs below follow Eitan's crossing method and its
antecedents [@ENT], using Rolle's theorem, determinants, convexity,
dominated convergence and the intermediate value theorem. Their role is
to make the logical route to the complete extremiser classification
explicit.

::: {#thm:R2 .theorem}
**Theorem 3** (Reduction with equality). *Let $1<a<b$ be real and
$X\in\mathcal L$. There exist $s_0\in[0,1]$ and $\lambda>0$ such that
$\mathbb E|X|^a=\mathbb E|\lambda Z_{s_0}|^a$ and
$$\mathbb E|X|^b\le\mathbb E|\lambda Z_{s_0}|^b,$$ with equality if and
only if $X\overset{d}{=}\lambda Z_{s_0}$. Consequently the maximum below
exists and
$$\frac{\|X\|_b}{\|X\|_a}\le\max_{0\le s\le1}\frac{\|Z_s\|_b}{\|Z_s\|_a},$$
with equality if and only if $X\overset{d}{=}\lambda Z_{s_1}$ for some
$\lambda>0$ and some $s_1$ at which the maximum is attained.*
:::

::: {#thm:R1 .theorem}
**Theorem 4** (Odd moments). *Let $1<a<r$ be real and $X\in\mathcal L$.
Then
$$\frac{\mathbb E[|X|^r\operatorname{sgn}X]}{\|X\|_a^{\,r}}\le\frac{\mathbb E[|E-1|^r\operatorname{sgn}(E-1)]}{\|E-1\|_a^{\,r}},$$
with equality if and only if $X\overset{d}{=}\lambda(E-1)$ for some
$\lambda>0$.*
:::

## Tools

::: {#lem:A .lemmaL}
**Lemma 1** (log-concave preliminaries). *Let $X\in\mathcal L$ have
density $e^{-V}$, where $V:\mathbb R\to(-\infty,\infty]$ is convex. Then
$I:=\{V<\infty\}$ is an interval with nonempty interior; $V$ is finite
and convex on $\operatorname{int}I$, and
$f:=e^{-V}\mathbf 1_{\operatorname{int}I}$ is a version of the density;
and there are $C,\kappa>0$ with $f(x)\le Ce^{-\kappa|x|}$ for all $x$.
In particular $\mathbb E|X|^t<\infty$ for every $t\ge0$.*
:::

::: proof
*Proof.* $I$ is convex and has positive Lebesgue measure because
$\int e^{-V}=1$, so it is an interval with nonempty interior; its
boundary has at most two points, a null set. Fix
$x_0\in\operatorname{int}I$ and $m\in\partial V(x_0)$ (the
subdifferential is nonempty at interior points). Then
$V(x)\ge V(x_0)+m(x-x_0)$ for all $x\in\mathbb R$, so
$f(x)\le e^{-V(x_0)}e^{-m(x-x_0)}$ everywhere; in particular $f$ is
bounded on every bounded set, and on a bounded side of $I$ it vanishes
beyond the endpoint.

Suppose $I$ is unbounded above. If $V(x)\le V(x_0)$ for all $x>x_0$,
then $e^{-V}\ge e^{-V(x_0)}$ on $(x_0,\infty)$, contradicting
integrability. So there is $x_1>x_0$ with $V(x_1)>V(x_0)$, and by
convexity, for $x\ge x_1$,
$$V(x)\ge V(x_1)+\kappa_+(x-x_1),\qquad \kappa_+=\frac{V(x_1)-V(x_0)}{x_1-x_0}>0.$$
The same argument applies below if $I$ is unbounded below. Hence on each
unbounded side $f$ decays exponentially beyond some point; on each
bounded side $f$ vanishes beyond the endpoint; and in between $f$ is
bounded by the affine-minorant bound. Combining these gives $C,\kappa>0$
with $f(x)\le Ce^{-\kappa|x|}$. ◻
:::

::: {#lem:B .lemmaL}
**Lemma 2** (Descartes' rule for generalised polynomials). *Let
$\lambda_0<\dots<\lambda_n$ be real and $c=(c_0,\dots,c_n)\neq0$. Then
$P(t)=\sum_ic_it^{\lambda_i}$ has at most $\mathrm V(c)$ distinct zeros
in $(0,\infty)$, where $\mathrm V(c)$ is the number of sign changes in
$c$ after zero entries are deleted.*
:::

::: proof
*Proof.* Induction on $\mathrm V(c)$. If $\mathrm V(c)=0$, all nonzero
terms $c_it^{\lambda_i}$ have the same sign for $t>0$, and at least one
is nonzero, so $P(t)\ne0$. If $\mathrm V(c)\ge1$, choose consecutive
nonzero entries $c_j,c_{j'}$ ($j<j'$) of opposite signs and
$\mu\in(\lambda_j,\lambda_{j'})$. The function $\psi(t)=t^{-\mu}P(t)$
has the same zeros as $P$ on $(0,\infty)$, and
$$\psi'(t)=\sum_ic_i(\lambda_i-\mu)t^{\lambda_i-\mu-1}.$$ Its exponents
are distinct, and its coefficients are nonzero exactly where $c_i\ne0$.
Passing from $c$ to these coefficients flips the sign of every entry
with $i\le j$ and keeps every entry with $i\ge j'$; this preserves the
sign changes inside each group and removes the one between $c_j$ and
$c_{j'}$. So the coefficients of $\psi'$ have $\mathrm V(c)-1$ sign
changes. If $P$ had $N$ distinct zeros, Rolle's theorem would give
$\psi'$ at least $N-1$ distinct zeros, and by induction
$N-1\le\mathrm V(c)-1$. ◻
:::

This is the classical rule of Laguerre [@PolyaSzego Part V, Ch. 1,
Sec. 6]; the proof is included for completeness.

::: {#lem:C .lemmaL}
**Lemma 3**. *For $c\in\mathbb R^{n+1}\setminus\{0\}$ let
$c^*_i=(-1)^ic_i$. Then
$\mathrm V(c)+\mathrm V(c^*)\le n-\mathbf 1[c_0=0]$.*
:::

::: proof
*Proof.* Consider consecutive nonzero entries at indices $i<i'$, with
gap $g=i'-i\ge1$. This pair contributes
$\mathbf 1[c_ic_{i'}<0]+\mathbf 1[(-1)^gc_ic_{i'}<0]$ to the left side.
If $g$ is odd the contribution is exactly $1\le g$; if $g$ is even it is
$0$ or $2\le g$. Summing over pairs gives at most
$i_{\rm last}-i_{\rm first}\le n-\mathbf 1[c_0=0]$. ◻
:::

::: {#prop:D .propL}
**Proposition 4** (Chebyshev property). *Let $1<a<r<b$ and set
$$u_0=1,\quad u_1=x,\quad u_2=|x|^a,\quad u_3=|x|^r\operatorname{sgn}x,\quad u_4=|x|^b .$$
These are continuous on $\mathbb R$. For each $n\in\{0,\dots,4\}$, every
nontrivial combination $\varphi=\sum_{i\le n}c_iu_i$ has at most $n$
distinct zeros in $\mathbb R$.*
:::

::: proof
*Proof.* Let $(\lambda_0,\dots,\lambda_4)=(0,1,a,r,b)$. For $t>0$,
$u_i(t)=t^{\lambda_i}$ and $u_i(-t)=(-1)^it^{\lambda_i}$. So on
$(0,\infty)$, $\varphi=\sum c_it^{\lambda_i}$ has at most $\mathrm V(c)$
distinct zeros by Lemma [2](#lem:B){reference-type="ref"
reference="lem:B"}; $t\mapsto\varphi(-t)=\sum c_i^*t^{\lambda_i}$ has at
most $\mathrm V(c^*)$ distinct zeros in $(0,\infty)$; and $0$ is a zero
if and only if $c_0=0$. Lemma [3](#lem:C){reference-type="ref"
reference="lem:C"} bounds the total by $n$. ◻
:::

The same proof applies to any number of functions with alternating
parity and increasing exponents $0<1<\lambda_2<\cdots$.

::: {#lem:E .lemmaL}
**Lemma 5** (alternating kernels). *Let $u_0,\dots,u_n$ be continuous on
$\mathbb R$, and suppose every nontrivial combination of them has at
most $n$ distinct zeros. Given $z_1<\dots<z_n$, there is a nontrivial
combination $\varphi$ of $u_0,\dots,u_n$ with
$\{\varphi=0\}=\{z_1,\dots,z_n\}$ that has constant nonzero sign on each
interval $(z_j,z_{j+1})$ ($z_0=-\infty$, $z_{n+1}=+\infty$), these signs
alternating in $j$.*
:::

::: proof
*Proof.* For $t_0<\dots<t_n$ let
$D(t_0,\dots,t_n)=\det[u_i(t_j)]_{i,j=0}^n$. It never vanishes:
otherwise some $c\ne0$ satisfies $\sum_ic_iu_i(t_j)=0$ for all $j$,
giving a nontrivial combination with $n+1$ zeros. Being continuous on
the connected set $\{t_0<\dots<t_n\}$, $D$ has a constant sign $\sigma$.
Put $\varphi(x)=\det[u_i(z_1),\dots,u_i(z_n),u_i(x)]$, a combination of
the $u_i$ by expansion along the last column. It vanishes at each $z_j$.
For $x\in(z_j,z_{j+1})$, moving the last column to position $j$ takes
$n-j$ transpositions, so
$\varphi(x)=(-1)^{n-j}D(z_1,\dots,z_j,x,z_{j+1},\dots,z_n)$, whose sign
is $(-1)^{n-j}\sigma\ne0$. In particular $\varphi$ is nontrivial. ◻
:::

::: {#lem:F .lemmaL}
**Lemma 6** (sign patterns and moment constraints). *Let $1<a<r<b$, let
$u_i$ and $\lambda_i$ be as in
Proposition [4](#prop:D){reference-type="ref" reference="prop:D"}, and
let $n\in\{3,4\}$. Let $h$ be measurable and not a.e. zero, with
$\int|h|(1+|x|^{\lambda_n})<\infty$ and $\int hu_i=0$ for
$i=0,\dots,n-1$. Suppose there are
$-\infty=y_0\le y_1\le\dots\le y_n\le y_{n+1}=+\infty$ such that
$(-1)^{n-j}h\ge0$ a.e. on $(y_j,y_{j+1})$ for $j=0,\dots,n$ (alternating
signs, the last one $+$). Then $\int hu_n>0$.*
:::

::: proof
*Proof.* *Step 1: a minimal alternating pattern.* Let $J=\{j:h$ is not
a.e. zero on $(y_j,y_{j+1})\}=\{j_1<\dots<j_m\}$, which is nonempty, and
let $e_l=(-1)^{n-j_l}$. Group the indices of $J$ into maximal runs of
constant $e_l$; if $j_l$ is the last index of a run that is not the
final run, put a change point at $y_{j_l+1}$. These points
$z_1<\dots<z_k$ are finite and strictly increasing:
$(y_{j_l},y_{j_l+1})$ has positive measure, so
$y_{j_l+1}>y_{j_l}\ge-\infty$, and
$y_{j_l+1}\le y_{j_{l+1}}<y_{j_{l+1}+1}\le+\infty$. Moreover
$k\le m-1\le n$. Let $\omega$ be the step function equal to $e_1$ on
$(-\infty,z_1)$ that switches sign at each $z_i$. Every interval
$(y_j,y_{j+1})$ with $j\in J$ lies between the two change points that
bound its run, so $\omega h\ge0$ a.e., and $\omega h$ is not a.e. zero.

*Step 2: $k<n$ is impossible.* If $k\le n-1$, then $u_0,\dots,u_k$
satisfy the hypothesis of Lemma [5](#lem:E){reference-type="ref"
reference="lem:E"} by Proposition [4](#prop:D){reference-type="ref"
reference="prop:D"}. Lemma [5](#lem:E){reference-type="ref"
reference="lem:E"} gives
$\varphi\in\operatorname{span}\{u_0,\dots,u_k\}$ with, after a choice of
sign, $\varphi\omega>0$ off $\{z_i\}$. Then $\varphi h\ge0$ a.e. and is
not a.e. zero, so $\int\varphi h>0$; but $\int\varphi h=0$ by the
constraints.

*Step 3: $k=n$.* Then $m=n+1$ and $J=\{0,\dots,n\}$, so the final sign
is $e_m=+1$. Take $\varphi$ from Lemma [5](#lem:E){reference-type="ref"
reference="lem:E"} for $u_0,\dots,u_n$ with sign chosen so that
$\varphi>0$ on $(z_n,\infty)$. Then $\varphi$ has the sign of $\omega$
off $\{z_i\}$, so $\int\varphi h>0$. Writing
$\varphi=\sum_{i\le n}c_iu_i$, the constraints give
$\int\varphi h=c_n\int hu_n$, so $c_n\ne0$. As $x\to+\infty$,
$\varphi(x)/x^{\lambda_n}\to c_n$, because $u_n(x)=x^{\lambda_n}$ for
$x>0$ has the largest exponent; since $\varphi>0$ on $(z_n,\infty)$,
$c_n>0$. Hence $\int hu_n>0$. ◻
:::

::: {#lem:G .lemmaL}
**Lemma 7** (crossing pattern against $Z_s$). *Let $X\in\mathcal L$ have
the density $f=e^{-V}\mathbf 1_{\operatorname{int}I}$ of
Lemma [1](#lem:A){reference-type="ref" reference="lem:A"}. Let
$\lambda>0$, $s\in[0,1]$, let $g$ be the density of $\lambda Z_s$, and
let $h=g-f$.*

1.  *There are $-\infty\le y_1\le y_2\le y_3\le y_4\le+\infty$ such
    that, a.e., $h$ is $\ge0$, $\le0$, $\ge0$, $\le0$, $\ge0$ on the
    five successive intervals; this is the pattern of
    Lemma [6](#lem:F){reference-type="ref" reference="lem:F"} with
    $n=4$.*

2.  *If $s=1$, there are $y_1\le y_2\le y_3$ with $h\le0$, $\ge0$,
    $\le0$, $\ge0$ a.e. on the four successive intervals; this is the
    pattern of Lemma [6](#lem:F){reference-type="ref" reference="lem:F"}
    with $n=3$.*
:::

::: proof
*Proof.* Put $c=1-s$ and $d=1-2s$. For $0<s<1$, $Z_s=d+(sE-cE')$, and
$sE-cE'$ has density $e^{-t/s}$ for $t>0$ and $e^{t/c}$ for $t<0$ (the
tails are given in Section [3](#sec:reduction){reference-type="ref"
reference="sec:reduction"}). So $g>0$ everywhere and $\log g$ is affine
on $L=(-\infty,\lambda d]$ and on $R=[\lambda d,\infty)$. For $s=0$,
$Z_0=1-E'$ and $g$ is positive and log-affine on $(-\infty,\lambda]$ and
zero after; for $s=1$, $Z_1=E-1$ and $g$ is positive and log-affine on
$[-\lambda,\infty)$ and zero before.

*Key observation.* Let $P$ be an interval on which $g=e^{\ell}$ with
$\ell$ affine. On $P\cap\operatorname{int}I$, $h<0$ if and only if
$\ell+V<0$; as $\ell+V$ is convex there, this set is an interval $N_P$
(possibly empty, possibly all of $P\cap\operatorname{int}I$, possibly
half-infinite). On $P\setminus\overline I$, $h=g>0$, and $\partial I$ is
null. Hence, a.e. on $P$, $h<0$ exactly on the interval $N_P$ and
$h\ge0$ on the rest of $P$, which consists of at most two intervals, one
on each side of $N_P$. Where $g=0$, $h=-f\le0$. Whether the endpoints of
the intervals below are included is immaterial, since all statements are
a.e.

The patterns follow.

- $0<s<1$: let $(y_1,y_2)$ be the endpoints of $N_L$ and $(y_3,y_4)$
  those of $N_R$, so $y_2\le\lambda d\le y_3$. If $N_L=\emptyset$ take
  $y_1=y_2=\lambda d$; if $N_R=\emptyset$ take $y_3=y_4=\lambda d$.

- $s=0$: let $(y_1,y_2)$ be the endpoints of $N_{(-\infty,\lambda]}$,
  $y_3=\lambda$ and $y_4=+\infty$ (the fifth interval is empty); if
  $N_{(-\infty,\lambda]}=\emptyset$ take $y_1=y_2=\lambda$.

- $s=1$, part 1: $y_1=-\infty$, $y_2=-\lambda$, and $(y_3,y_4)$ the
  endpoints of $N_{[-\lambda,\infty)}$; if it is empty take
  $y_3=y_4=-\lambda$.

- $s=1$, part 2: $y_1=-\lambda$ and $(y_2,y_3)$ the endpoints of
  $N_{[-\lambda,\infty)}$; if it is empty take $y_2=y_3=-\lambda$.

In each case the stated sign holds a.e. on each interval. ◻
:::

## Proofs of Theorems [3](#thm:R2){reference-type="ref" reference="thm:R2"} and [4](#thm:R1){reference-type="ref" reference="thm:R1"} {#proofs-of-theorems-thmr2-and-thmr1}

::: proof
*Proof of Theorem [4](#thm:R1){reference-type="ref"
reference="thm:R1"}.* Let $Y=\lambda(E-1)=\lambda Z_1$ with $\lambda>0$
chosen so that $\mathbb E|Y|^a=\mathbb E|X|^a$, let $g$ be its density
and $h=g-f$. By Lemma [1](#lem:A){reference-type="ref"
reference="lem:A"} all moments are finite, so $\int|h|(1+|x|^r)<\infty$.
We have $\int h=0$ (both are densities), $\int xh=0$ (both are centred)
and $\int|x|^ah=0$ (choice of $\lambda$); these are $\int hu_i=0$ for
$i=0,1,2$, with $u_i$ as in
Proposition [4](#prop:D){reference-type="ref" reference="prop:D"} for
any $b>r$. By Lemma [7](#lem:G){reference-type="ref"
reference="lem:G"}(2), $h$ has the $n=3$ pattern. If $h$ is not
a.e. zero, Lemma [6](#lem:F){reference-type="ref" reference="lem:F"}
with $n=3$ gives $\int|x|^r\operatorname{sgn}(x)\,h>0$, that is,
$\mathbb E[|Y|^r\operatorname{sgn}Y]>\mathbb E[|X|^r\operatorname{sgn}X]$.
Dividing by $\|X\|_a^r=\|Y\|_a^r$ and using the scale invariance of the
ratio gives strict inequality unless $X\overset{d}{=}Y$. ◻
:::

::: proof
*Proof of Theorem [3](#thm:R2){reference-type="ref"
reference="thm:R2"}.* Fix $r=(a+b)/2$ and let
$\rho(Y)=\mathbb E[|Y|^r\operatorname{sgn}Y]/(\mathbb E|Y|^a)^{r/a}$. It
is scale invariant and $\rho(-Y)=-\rho(Y)$.
Theorem [4](#thm:R1){reference-type="ref" reference="thm:R1"} applied to
$X$ and to $-X$ (both in $\mathcal L$) gives $\rho(X)\le\rho(Z_1)$ and
$-\rho(X)\le\rho(Z_1)$. Since $Z_0\overset{d}{=}-Z_1$,
$\rho(Z_0)=-\rho(Z_1)$, so $$\rho(Z_0)\le\rho(X)\le\rho(Z_1).$$ The map
$s\mapsto Z_s$ is pointwise continuous and $|Z_s|\le E+E'+2$. For each
$t\ge0$, $(E+E'+2)^t$ is integrable, so by dominated convergence
$s\mapsto\mathbb E|Z_s|^t$ and
$s\mapsto\mathbb E[|Z_s|^r\operatorname{sgn}Z_s]$ are continuous on
$[0,1]$; moreover $\mathbb E|Z_s|^t>0$. Hence $s\mapsto\rho(Z_s)$ is
continuous, and the intermediate value theorem gives $s_0\in[0,1]$ with
$\rho(Z_{s_0})=\rho(X)$.

Let $\lambda=\|X\|_a/\|Z_{s_0}\|_a$ and $Y=\lambda Z_{s_0}$, with
density $g_Y$. Then $\mathbb E|Y|^a=\mathbb E|X|^a$ and
$\mathbb E[|Y|^r\operatorname{sgn}Y]=\rho(X)(\mathbb E|X|^a)^{r/a}=\mathbb E[|X|^r\operatorname{sgn}X]$.
With $h=g_Y-f$ we have $\int|h|(1+|x|^b)<\infty$ and $\int hu_i=0$ for
$i=0,\dots,3$. By Lemma [7](#lem:G){reference-type="ref"
reference="lem:G"}(1), $h$ has the $n=4$ pattern, and
Lemma [6](#lem:F){reference-type="ref" reference="lem:F"} gives
$\mathbb E|Y|^b-\mathbb E|X|^b=\int|x|^bh>0$ unless $h=0$ a.e., that is,
unless $X\overset{d}{=}Y$.

*Consequences.* Let $Q(s)=\|Z_s\|_b/\|Z_s\|_a$. By the continuity just
proved (with $t=a$ and $t=b$), $Q$ is continuous on the compact interval
$[0,1]$ and attains its maximum. Then
$$\frac{\|X\|_b}{\|X\|_a}\le\frac{\|Y\|_b}{\|Y\|_a}=Q(s_0)\le\max Q.$$
If $\|X\|_b/\|X\|_a=\max Q$, both inequalities are equalities: the first
forces $X\overset{d}{=}\lambda Z_{s_0}$, and the second shows that $s_0$
is a maximiser. Conversely, if $X\overset{d}{=}\lambda Z_{s_1}$ with
$s_1$ a maximiser, equality is immediate. ◻
:::

::: remark
*Remark 5*. If $s_0\in\{0,1\}$ in the proof of
Theorem [3](#thm:R2){reference-type="ref" reference="thm:R2"}, then $h$
has at most three sign changes, and
Lemma [6](#lem:F){reference-type="ref" reference="lem:F"} forces $h=0$
a.e. For $s_0=1$, the pattern of Lemma [7](#lem:G){reference-type="ref"
reference="lem:G"}(2) ($-,+,-,+$) is the $n=4$ pattern of
Lemma [6](#lem:F){reference-type="ref" reference="lem:F"} with an empty
first interval, and Step 2 of the proof of
Lemma [6](#lem:F){reference-type="ref" reference="lem:F"} ($k\le3<4$)
rules out $h\ne0$. For $s_0=0$ the pattern is $+,-,+,-$, so
Lemma [6](#lem:F){reference-type="ref" reference="lem:F"} is applied in
the same way to $-h$, which satisfies the same four constraints. This is
consistent with the equality case of
Theorem [4](#thm:R1){reference-type="ref" reference="thm:R1"}.
:::

# Reduction to a two-parameter inequality {#sec:reduction}

Each $Z_s$ is itself in $\mathcal L$ (for $0<s<1$ it has the density of
an asymmetric Laplace law shifted to mean zero), so the maximum in
Theorem [3](#thm:R2){reference-type="ref" reference="thm:R2"} is the
sharp constant. Applying Theorem [3](#thm:R2){reference-type="ref"
reference="thm:R2"} with $(a,b)=(p,2)$ for $1<p<2$ and $(2,p)$ for
$p>2$, and using $Z_{1-s}\overset d=-Z_s$:
$$\log m_p=\min_{0\le s\le1/2}G(s,p)\quad(1<p<2),\qquad \log M_p=\max_{0\le s\le1/2}G(s,p)\quad(p>2),$$
where $G(s,p)=\frac1p\log M(s,p)-\frac12\log N(s)$,
$M(s,p)=\mathbb E|Z_s|^p$ and $N(s)=\mathbb EZ_s^2=s^2+(1-s)^2$. Since
$G(0,p)=\frac1p\log F(p)=:G_0(p)$ and
$G(\frac12,p)=\frac1p\log\Gamma(p+1)-\frac12\log2=:G_h(p)$,
Theorem [1](#thm:main){reference-type="ref" reference="thm:main"} is
equivalent to:

- $G(s,p)\ge G_h(p)$ for $1<p\le p^*$, $0\le s\le\frac12$;

- $G(s,p)\ge G_0(p)$ for $p^*\le p<2$, $0\le s\le\frac12$;

- $G(s,p)\le G_0(p)$ for $p>2$, $0\le s\le\frac12$;

together with the existence and uniqueness of $p^*$.

#### A representation.

Put $c=1-s$, $d=1-2s$, and let $R$ be standard exponential. For
$0<s\le\frac12$, $sE-cE'$ has tails $\mathbb P(sE-cE'>t)=se^{-t/s}$ and
$\mathbb P(sE-cE'<-t)=ce^{-t/c}$ ($t\ge0$), the same as the mixture
equal to $sR$ with probability $s$ and $-cR$ with probability $c$. Since
$Z_s=sE-cE'+d$, $$\begin{equation}
\label{eq:rep}
M(s,p)=s\,\mathbb E(d+sR)^p+c\,\mathbb E|d-cR|^p,
\end{equation}$$ and
$\mathbb E|x-R|^p=e^{-x}\bigl(\Gamma(p+1)+\int_0^xu^pe^u\,du\bigr)$ for
$x\ge0$. In particular
$F(p)=e^{-1}\bigl(\Gamma(p+1)+\int_0^1u^pe^u\,du\bigr)$.

# Large exponents {#sec:largep}

::: {#prop:largep .proposition}
**Proposition 6**. *For every real $p\ge7$ and $0<s\le\frac12$,
$$\begin{gathered}
N(s)^{p/2}F(p)-M(s,p)\;\ge\;\frac{s^2(1-s)^{p-1}}{e}\,K(p)>0,\\
K(p)=\Bigl(\frac p2-2e+2\Bigr)\Gamma(p+1)-\frac{e(p-2)}2-\frac{e(e-2)}{p+1}.
\end{gathered}$$ In particular $G(s,p)<G_0(p)$, so (T$_C$) holds for
$p\ge7$.*
:::

::: proof
*Proof.* Write $\Gamma_p=\Gamma(p+1)$ and $J(x)=\int_0^xu^pe^u\,du$. Put
$\tau=s/c\in(0,1]$, so that $c=1/(1+\tau)$, $d/c=1-\tau$ and
$N(s)=c^2(1+\tau^2)$. By [\[eq:rep\]](#eq:rep){reference-type="eqref"
reference="eq:rep"}, $$\begin{aligned}
M(s,p)&=\frac{c^p}{1+\tau}\Bigl[\tau\,\mathbb E(1-\tau+\tau R)^p+e^{\tau-1}\bigl(\Gamma_p+J(1-\tau)\bigr)\Bigr],\\
N(s)^{p/2}F(p)&=c^p(1+\tau^2)^{p/2}e^{-1}\bigl(\Gamma_p+J(1)\bigr).
\end{aligned}$$ With $A(\tau)=(1+\tau)(1+\tau^2)^{p/2}$ this gives
$e(1+\tau)c^{-p}\bigl[N^{p/2}F-M\bigr]=S(\tau)$, where
$$S(\tau)=A(\tau)\bigl(\Gamma_p+J(1)\bigr)-e\tau\,\mathbb E(1-\tau+\tau R)^p-e^\tau\bigl(\Gamma_p+J(1-\tau)\bigr).$$
Since $t\mapsto t^p$ is convex and $1-\tau+\tau R$ is a convex
combination of $1$ and $R$,
$\mathbb E(1-\tau+\tau R)^p\le(1-\tau)+\tau\Gamma_p$. Hence $S\ge\Delta$
with
$$\Delta(\tau)=\bigl[A(\tau)-e^\tau-e\tau^2\bigr]\Gamma_p+\Bigl\{A(\tau)J(1)-e^\tau J(1-\tau)-e\tau(1-\tau)\Bigr\}.$$
For $\tau\in[0,1]$ and $p\ge2$ we use:

- $A(\tau)\ge(1+\tau)(1+\frac p2\tau^2)\ge1+\tau+\frac p2\tau^2$, by
  Bernoulli's inequality with exponent $p/2\ge1$;

- $0\le e^\tau-1-\tau\le(e-2)\tau^2$;

- $0\le J(1)\le e/(p+1)$;

- $\int_{1-\tau}^1u^pe^u\,du\ge e^{1-\tau}\bigl(1-(1-\tau)^{p+1}\bigr)/(p+1)$;

- $1-(1-\tau)^{p+1}\ge(p+1)\tau-\frac{(p+1)p}2\tau^2$, by Taylor's
  theorem, since the third derivative of $(1-\tau)^{p+1}$ is $\le0$ on
  $[0,1]$.

By (a) and (b), the first bracket is at least $(\frac p2-2e+2)\tau^2$.
For the second, write $J(1-\tau)=J(1)-\int_{1-\tau}^1u^pe^u\,du$ and use
(a), (b), (c), (d), (e):
$$A(\tau)J(1)-e^\tau J(1-\tau)\ge(1+\tau-e^\tau)J(1)+e\,\frac{1-(1-\tau)^{p+1}}{p+1}\ge-\frac{e(e-2)}{p+1}\tau^2+e\Bigl(\tau-\frac p2\tau^2\Bigr).$$
Subtracting $e\tau(1-\tau)$, the terms linear in $\tau$ cancel, and the
second bracket is at least
$-\bigl(\frac{e(e-2)}{p+1}+\frac{e(p-2)}2\bigr)\tau^2$. Hence
$S(\tau)\ge\Delta(\tau)\ge\tau^2K(p)$, which is the stated bound after
multiplying by $c^p/(e(1+\tau))$.

Finally let $p\ge7$. Then
$\frac p2-2e+2>\frac{11}2-\frac{87}{16}=\frac1{16}$ (as
$e<\frac{87}{32}$). Since $\Gamma$ is increasing on $[2,\infty)$,
$\Gamma_p=p\,\Gamma(p)\ge720p$. Also $\frac{e(p-2)}2<\frac{87}{64}p$ and
$\frac{e(e-2)}{p+1}<\frac{87}{32}\cdot\frac{23}{32}\cdot\frac18=\frac{2001}{8192}$.
Thus
$$K(p)>45p-\frac{87}{64}p-\frac{2001}{8192}=\frac{2793}{64}p-\frac{2001}{8192}>0. \qedhere$$ ◻
:::

::: remark
*Remark 7*. The argument needs $\frac p2>2e-2$,
i.e. $p>4e-4\approx6.87$; indeed $K(6.8)<0<K(6.9)$. The loss comes from
the chord bound for $\mathbb E(1-\tau+\tau R)^p$.
:::

# Derivatives in $s$ and the endpoints {#sec:local}

Let $W=E+E'-2=\partial_sZ_s$.

::: {#lem:cond .lemma}
**Lemma 8**. *Let $0<s<1$ and $X=sE-cE'$, so that $Z_s=d+X$. Then
$\mathbb P(X>0)=s$, and conditionally on $\{X>0\}$ the pair $(Z_s,W)$
has the law of $(d+sR,\,(R-1)+(T-1))$. Conditionally on $\{X<0\}$ it has
the law of $(d-cR,\,(R-1)+(T-1))$. Here $R,T$ are independent standard
exponential random variables. At $s=0,1$ use the direct formulas for
$Z_s,W$, or their moment limits; no conditioning on a zero-probability
event is intended.*
:::

::: proof
*Proof.* On $\{sE>cE'\}$, condition on $E'$. By memorylessness
$R:=E-cE'/s$ is standard exponential and independent of $E'$, while $E'$
has conditional density proportional to $e^{-u}e^{-cu/s}=e^{-u/s}$ at
$u$, so $E'=sT$. Then $X=sR$ and $W=R+cT+sT-2$. The other side is
symmetric. ◻
:::

::: {#lem:deriv .lemma}
**Lemma 9**. *For $p>1$, $G(\cdot,p)$ is $C^2$ on $[0,\frac12]$. With
$\tau=-d/N$ and $V=W-\tau Z_s$, put $A=\mathbb E[|Z_s|^{p-2}Z_sV]/M$ and
$B=\mathbb E[|Z_s|^{p-2}V^2]/M$. Then
$$\partial_sG=A,\qquad \partial_s^2G=(p-1)B-pA^2+\frac{2d}{N}A-\frac1{N^2}.$$*
:::

::: proof
*Proof.* The map $s\mapsto Z_s$ is affine with slope $W$.
Differentiation under the expectation is justified by dominated
convergence for the first derivative. For the second, use the absolute
continuity of $z\mapsto|z|^{p-2}z$, whose derivative $(p-1)|z|^{p-2}$ is
locally integrable for $p>1$, together with Tonelli; finiteness comes
from the closed forms in supplement `local.md`. This gives
$M_s=p\mathbb E[|Z|^{p-2}ZW]$ and
$M_{ss}=p(p-1)\mathbb E[|Z|^{p-2}W^2]$. Substituting $W=V+\tau Z$ and
using $\tau'+\tau^2=1/N^2$ (as $d^2+4sc=1$) gives the formulas. ◻
:::

::: {#prop:endpoints .proposition}
**Proposition 10**. *For every $p>1$:*

- *$\partial_sG(0,p)=0$ and $\partial_sG(\frac12,p)=0$;*

- *$\partial_s^2G(0,p)=\frac{p-1}p\bigl(\frac1{F(p)}-1\bigr)$ and
  $\partial_s^2G(\frac12,p)=\frac{4(p-2)^2}p$.*
:::

::: proof
*Proof.* At $s=0$, $c=d=1$ and $V=W+Z=E-1$ is independent of $Z=1-E'$,
so $\mathbb E[|Z|^{p-2}ZV]=0$. Hence $A=0$ and $B=F(p-2)/F(p)$. Two
integrations by parts against $e^{-t}$ give $p(p-1)F(p-2)=F(p)-1+p$,
which yields the first formula in (ii).

At $s=\frac12$, $Z_{1-s}\overset d=-Z_s$ (exchange $E$ and $E'$) and
$N(1-s)=N(s)$, so $G(\cdot,p)$ is even about $\frac12$, and
$\partial_sG(\frac12,p)=0$. There $d=0$, $N=\frac12$,
$M=2^{-p}\Gamma(p+1)$, $|Z|=R/2$ on both sides, and
$\mathbb E[V\mid\cdot]=R-1$. So
$(p-1)B=4(p-1)\bigl(\Gamma(p+1)-2\Gamma(p)+2\Gamma(p-1)\bigr)/\Gamma(p+1)=4(p-3+4/p)$,
and $\partial_s^2G=4(p-3+4/p)-4=4(p-2)^2/p$. ◻
:::

::: {#prop:local .proposition}
**Proposition 11** (computer-assisted).

- *$\partial_s^2G>\frac{6}{125}$ on
  $[\frac9{20},\frac12]\times(1,\frac74]$.*

- *$\partial_s^2G>\frac{29}{2000}$ on
  $[0,\frac1{20}]\times[\frac85,\frac{39}{20}]$.*

- *$\partial_s^2G<-\frac{31}{2000}$ on
  $[0,\frac1{20}]\times[\frac{41}{20},7]$.*

*Consequently, $G\ge G_h$ on $[\frac9{20},\frac12]\times(1,\frac74]$,
$G\ge G_0$ on $[0,\frac1{20}]\times[\frac85,\frac{39}{20}]$, and
$G\le G_0$ on $[0,\frac1{20}]\times[\frac{41}{20},7]$.*
:::

::: proof
*Certificate.* By Lemma [8](#lem:cond){reference-type="ref"
reference="lem:cond"}, the quantities $M$, $\mathbb E[|Z|^{p-2}ZV]$ and
$(p-1)\mathbb E[|Z|^{p-2}V^2]$ are explicit combinations of
$\Gamma(p+1)$, $\Gamma(p)$, upper incomplete gamma functions,
$\int_0^xu^ae^u\,du=\frac{x^{a+1}}{a+1}{}_1F_1(a+1;a+2;x)$ and
${}_1F_1(1;b;z)$ (supplement `local.md`, Lemma 4 and Proposition B). The
terms that are singular as $p\to1$ appear with a factor $p-1$ and are
rewritten as
$d^{p-1}\bigl[\frac N{c^2}{}_1F_1(1;p;-x)-\frac N{s^2}{}_1F_1(1;p;d/s)\bigr]$,
which is continuous at the corner $(\frac12,1)$, where
$\partial_s^2G=4$.

Each region is covered by boxes with exact rational endpoints. On each
box, $p$ is a ball, the $s$-dependence is enclosed through monotone
hulls of the special-function pieces, and $d^{p-1}$ is enclosed in
$[d_0^{e_1},d_1^{e_0}]$. Precision is $128$ bits, boxes are bisected
down to width $2^{-22}$, and an undecided comparison counts as a
failure. Region L1 used $50{,}214$ boxes, L2 $129{,}875$ and L3
$57{,}496$, all without failure.

The consequences follow from
Proposition [10](#prop:endpoints){reference-type="ref"
reference="prop:endpoints"}(i) and Taylor's theorem: $\partial_sG$ has a
fixed sign on each strip away from the endpoint. ◻
:::

# The neighbourhood of $p=2$ {#sec:strip}

At $p=2$ we have $M(s,2)=N(s)$, so $G(\cdot,2)\equiv0$ and every
comparison degenerates. Let
$\mathcal H(s,p)=\partial_p\bigl[G(s,p)-G_0(p)\bigr]$.

::: {#prop:strip .proposition}
**Proposition 12** (computer-assisted).

- *$\mathcal H<0$ on
  $[\frac1{20},\frac12]\times[\frac{39}{20},\frac{41}{20}]$.*

- *$\partial_p\partial_s^2G<0$ on
  $[0,\frac1{20}]\times[\frac{39}{20},\frac{41}{20}]$.*

*Consequently, for $p\in[\frac{39}{20},\frac{41}{20}]$ and
$0<s\le\frac12$, $G(s,p)-G_0(p)$ has the sign of $2-p$.*
:::

::: proof
*Proof.* *Consequence.* Since $G(0,\cdot)=G_0$, $\mathcal H(0,p)=0$. The
mixed partial derivatives of $G$ up to $\partial_p\partial_s^2$ exist,
are continuous and commute on the strip (supplement `strip.md`,
Lemma 3), so $\partial_s\mathcal H(0,p)=\partial_p\partial_sG(0,p)=0$ by
Proposition [10](#prop:endpoints){reference-type="ref"
reference="prop:endpoints"}(i), and
$\partial_s^2\mathcal H=\partial_p\partial_s^2G<0$ on $[0,\frac1{20}]$
by (S2), so Taylor's formula gives $\mathcal H<0$ on $(0,\frac1{20}]$.
Together with (S1), $\mathcal H<0$ for $0<s\le\frac12$. Since
$G(s,2)-G_0(2)=0$, we have $G(s,p)-G_0(p)=\int_2^p\mathcal H(s,t)\,dt$,
which has the sign of $2-p$.

*Certificate.* Write $q=p-2$, $\Phi_A=\mathbb E|Z_s|^{2+q}$,
$\Phi_B=\mathbb E[|Z_s|^qZ_s(R-1)]$ and
$\Phi_C=\mathbb E[|Z_s|^q((R-1)^2+1)]$, the expectations being taken in
the conditioned form of Lemma [8](#lem:cond){reference-type="ref"
reference="lem:cond"}. Then $M_s=p\Phi_B$, $M_{ss}=p(p-1)\Phi_C$,
$\partial_s\partial_pG=\partial_q(\Phi_B/\Phi_A)$ and
$\partial_p\partial_s^2G=\partial_q\bigl[(1+q)\Phi_C/\Phi_A-(2+q)(\Phi_B/\Phi_A)^2\bigr]$
(with the $N$-terms independent of $p$).

On the negative side ($Z\le d$) the moments have closed forms in $q$ via
Gamma functions and a rapidly convergent series with bounded tail. On
the positive side they are enclosed by order-$12$ Taylor polynomials in
$q$, whose coefficients are rigorous integrals (Arb's `acb.integral`
with explicit tail and endpoint bounds), plus an explicit remainder
bound.

S1 is certified with a mean-value form in $s$ and $q$ ($157$ boxes,
worst certified upper bound below $-1.5\times10^{-5}$). S2 is certified
with direct evaluation in $s$ and a mean-value form in $q$ ($64$ boxes,
worst bound below $-0.040$). The shipped outputs retain exact rational
enclosure endpoints. Both certificates are regenerated in normal and
optimised replay; the regression suite rejects altered signs, empty or
overlapping coverage, and malformed domains. ◻
:::

# The crossing point and the bulk {#sec:bulk}

Let $H(p)=G_0(p)-G_h(p)$, where $G_0(p)=\frac1p\log F(p)$ and
$G_h(p)=\frac1p\log\Gamma(p+1)-\frac12\log2$.

::: {#prop:crossing .proposition}
**Proposition 13** (computer-assisted). *Let $p^*_{\rm lo}=1.682457208$
and $p^*_{\rm hi}=1.682457209$. Then:*

- *$H>0$ on $[1,p^*_{\rm lo}]$;*

- *$H<0$ on $[p^*_{\rm hi},\frac{39}{20}]$;*

- *$H'<-0.016$ on $[p^*_{\rm lo},p^*_{\rm hi}]$.*

*Hence $H$ has exactly one zero $p^*$ in $[1,\frac{39}{20}]$, $H>0$ on
$[1,p^*)$, and $H<0$ on $(p^*,\frac{39}{20}]$.*
:::

::: {#prop:bulk .proposition}
**Proposition 14** (computer-assisted). *Let $\eta=\frac1{20}$. The
following hold:*

- *$G(s,p)>G_h(p)$ for $s\in[0,\frac9{20}]$,
  $p\in[1,p^*_{\rm hi}-\eta]$;*

- *$G(s,p)>\min\{G_0(p),G_h(p)\}$ for $s\in[\frac1{20},\frac9{20}]$,
  $p\in[p^*_{\rm lo}-\eta,p^*_{\rm hi}+\eta]$;*

- *$G(s,p)>G_0(p)$ for $s\in[\frac1{20},\frac12]$,
  $p\in[p^*_{\rm lo}+\eta,\frac{39}{20}]$;*

- *$G(s,p)<G_0(p)$ for $s\in[\frac1{20},\frac12]$,
  $p\in[\frac{41}{20},7]$.*
:::

::: proof
*Certificates.*
$M(s,p)=s\,J_p+c^{p+1}e^{-x}\bigl(\Gamma(p+1)+\frac{x^{p+1}}{p+1}{}_1F_1(p+1;p+2;x)\bigr)$,
where $x=d/c$ and $J_q=s^qe^{d/s}\Gamma(q+1,d/s)$. The partial
derivatives $M_s$ and $M_p=\mathbb E|Z_s|^p\log|Z_s|$ have explicit
forms (supplement `bulk.md`).

On a box with centre $(s_0,p_0)$ and half-widths $(h_s,h_p)$, the
difference $D$ to be certified is enclosed by
$D(s_0,p_0)+[\partial_pD]\,[-h_p,h_p]+[\partial_sD]\,[-h_s,h_s]$. Here
$[\partial_pD]$ is bounded through the convexity of
$p\mapsto\mathbb E|Z_s|^p$: $M_p=\mathbb E|Z_s|^p\log|Z_s|$ and $F'$ are
nondecreasing in $p$ (their derivatives are
$\mathbb E|Z_s|^p\log^2|Z_s|\ge0$), so their values at the exact
rational endpoints bound them; $M$ and $F$ are enclosed separately and
interval division is used. $[\partial_sD]$ is bounded from monotone
enclosures of the closed-form pieces in $s$. No incomplete gamma
function is evaluated at a wide parameter ball.

Branch and bound at $128$ bits gives:

- B1: $1653$ boxes, margin $\ge3.2\times10^{-6}$;

- B2: $1127$ boxes, margin $\ge4.0\times10^{-8}$;

- B3: $4413$ boxes, margin $\ge1.4\times10^{-8}$;

- B4: $8635$ boxes, margin $\ge8.0\times10^{-8}$;

- the crossing lemma: $73$, $65$ and $1$ boxes.

There were no failures. B1--B4 are regenerated at $128$ bits in both
normal and optimised replay. The four deliberately false bulk statements
are rejected in both modes. Optional floating-point quadrature checks
are diagnostic only and do not replace the interval certificates. ◻
:::

# Proof of Theorem [1](#thm:main){reference-type="ref" reference="thm:main"} {#sec:assembly}

By Section [3](#sec:reduction){reference-type="ref"
reference="sec:reduction"} it suffices to prove (T$_A$), (T$_B$),
(T$_C$) and the uniqueness of $p^*$. Throughout, $s\in[0,\frac12]$.

#### The crossing.

By Proposition [13](#prop:crossing){reference-type="ref"
reference="prop:crossing"}, $H=G_0-G_h$ has a unique zero $p^*$ in
$[1,\frac{39}{20}]$, with $H>0$ before it and $H<0$ after it. On
$[\frac{39}{20},2)$, Proposition [12](#prop:strip){reference-type="ref"
reference="prop:strip"} at $s=\frac12$ gives $G_h-G_0>0$, so $p^*$ is
the unique zero in $(1,2)$.

#### (T$_A$): $1<p\le p^*$.

- For $p\le p^*_{\rm hi}-\eta$: B1 covers $s\le\frac9{20}$ and L1 covers
  $s\ge\frac9{20}$.

- For $p\in[p^*_{\rm hi}-\eta,p^*]$:

  - for $s\le\frac1{20}$, L2 gives $G\ge G_0$, and $G_0\ge G_h$ because
    $p\le p^*$;

  - for $s\in[\frac1{20},\frac9{20}]$, B2 gives $G>\min\{G_0,G_h\}=G_h$;

  - for $s\ge\frac9{20}$, L1 applies.

#### (T$_B$): $p^*\le p<2$.

- For $p\in[p^*,p^*_{\rm lo}+\eta]$:

  - L2 covers $s\le\frac1{20}$;

  - B2 gives $G>\min\{G_0,G_h\}=G_0$ on $[\frac1{20},\frac9{20}]$;

  - L1 gives $G\ge G_h\ge G_0$ on $[\frac9{20},\frac12]$.

- For $p\in[p^*_{\rm lo}+\eta,\frac{39}{20}]$: L2 covers
  $s\le\frac1{20}$ and B3 the rest.

- For $p\in[\frac{39}{20},2)$:
  Proposition [12](#prop:strip){reference-type="ref"
  reference="prop:strip"}.

#### (T$_C$): $p>2$.

- $(2,\frac{41}{20}]$:
  Proposition [12](#prop:strip){reference-type="ref"
  reference="prop:strip"}.

- $[\frac{41}{20},7]$: L3 for $s\le\frac1{20}$ and B4 for the rest.

- $[7,\infty)$: Proposition [6](#prop:largep){reference-type="ref"
  reference="prop:largep"}. 0◻

# Equality cases: proof of Theorem [2](#thm:U){reference-type="ref" reference="thm:U"} {#sec:equality}

#### Reduction to $G$.

Recall $G(1-s,p)=G(s,p)$, since $Z_{1-s}\overset{d}{=}-Z_s$ and
$N(1-s)=N(s)$. For $1<p<2$ apply
Theorem [3](#thm:R2){reference-type="ref" reference="thm:R2"} with
$(a,b)=(p,2)$, where $\|Z_s\|_2/\|Z_s\|_p=e^{-G(s,p)}$; for $p>2$ apply
it with $(a,b)=(2,p)$, where $\|Z_s\|_p/\|Z_s\|_2=e^{G(s,p)}$. By the
equality clause of Theorem [3](#thm:R2){reference-type="ref"
reference="thm:R2"}, equality in
Theorem [1](#thm:main){reference-type="ref" reference="thm:main"} holds
for $X\in\mathcal L$ if and only if $X\overset{d}{=}\lambda Z_{s_1}$ for
some $\lambda>0$ and some $s_1\in[0,1]$ that minimises $G(\cdot,p)$ on
$[0,1]$ when $1<p<2$, or maximises it when $p>2$.
Theorem [2](#thm:U){reference-type="ref" reference="thm:U"} therefore
follows from Lemma [8](#lem:S){reference-type="ref" reference="lem:S"}
below together with this symmetry.
Lemma [8](#lem:S){reference-type="ref" reference="lem:S"} shows that the
set of extremal parameters in $[0,1]$ is $\{\frac12\}$ in case (i),
$\{0,\frac12,1\}$ in case (ii), and $\{0,1\}$ in cases (iii) and (iv).
Finally, $\lambda Z_{1/2}=\frac\lambda2(E-E')$,
$\lambda Z_1=\lambda(E-1)$ and $\lambda Z_0\overset{d}{=}-\lambda(E-1)$.

#### Facts used.

Lemma [8](#lem:S){reference-type="ref" reference="lem:S"} uses only the
following statements, all established above with strict inequalities on
the stated closed sets. Here $\eta=\frac1{20}$,
$p^*_{\rm lo}=1.682457208$, $p^*_{\rm hi}=1.682457209$, and $G_0$,
$G_h$, $H=G_0-G_h$ are as in
Sections [3](#sec:reduction){reference-type="ref"
reference="sec:reduction"} and [7](#sec:bulk){reference-type="ref"
reference="sec:bulk"}.

- Lemma [9](#lem:deriv){reference-type="ref" reference="lem:deriv"}: for
  each $p>1$, $G(\cdot,p)$ is $C^2$ on $[0,\frac12]$.

- Proposition [10](#prop:endpoints){reference-type="ref"
  reference="prop:endpoints"}(i): for each $p>1$,
  $\partial_sG(0,p)=\partial_sG(\frac12,p)=0$.

- Proposition [11](#prop:local){reference-type="ref"
  reference="prop:local"}: (L1) $\partial_s^2G>\frac6{125}$ on
  $[\frac9{20},\frac12]\times(1,\frac74]$; (L2)
  $\partial_s^2G>\frac{29}{2000}$ on
  $[0,\frac1{20}]\times[\frac85,\frac{39}{20}]$; (L3)
  $\partial_s^2G<-\frac{31}{2000}$ on
  $[0,\frac1{20}]\times[\frac{41}{20},7]$.

- Proposition [12](#prop:strip){reference-type="ref"
  reference="prop:strip"}: for $p\in[\frac{39}{20},\frac{41}{20}]$ and
  $0<s\le\frac12$, $G(s,p)-G_0(p)=\int_2^p\mathcal H(s,t)\,dt$ with
  $\mathcal H<0$; hence $G(s,p)>G_0(p)$ for $\frac{39}{20}\le p<2$ and
  $G(s,p)<G_0(p)$ for $2<p\le\frac{41}{20}$.

- Proposition [13](#prop:crossing){reference-type="ref"
  reference="prop:crossing"}: $H>0$ on $[1,p^*_{\rm lo}]$, $H<0$ on
  $[p^*_{\rm hi},\frac{39}{20}]$ and $H'<0$ on
  $[p^*_{\rm lo},p^*_{\rm hi}]$. Hence $p^*_{\rm lo}<p^*<p^*_{\rm hi}$,
  $H>0$ on $[1,p^*)$, $H(p^*)=0$ and $H<0$ on $(p^*,\frac{39}{20}]$.

- Proposition [14](#prop:bulk){reference-type="ref"
  reference="prop:bulk"}, on closed rectangles: (B1) $G>G_h$ on
  $[0,\frac9{20}]\times[1,p^*_{\rm hi}-\eta]$; (B2) $G>\min\{G_0,G_h\}$
  on
  $[\frac1{20},\frac9{20}]\times[p^*_{\rm lo}-\eta,p^*_{\rm hi}+\eta]$;
  (B3) $G>G_0$ on
  $[\frac1{20},\frac12]\times[p^*_{\rm lo}+\eta,\frac{39}{20}]$; (B4)
  $G<G_0$ on $[\frac1{20},\frac12]\times[\frac{41}{20},7]$. The
  certified minimum margins are $3.2\times10^{-6}$, $4.0\times10^{-8}$,
  $1.4\times10^{-8}$ and $8.0\times10^{-8}$.

- Proposition [6](#prop:largep){reference-type="ref"
  reference="prop:largep"}: for $p\ge7$ and $0<s\le\frac12$,
  $N(s)^{p/2}F(p)-M(s,p)>0$. Since
  $G(s,p)-G_0(p)=\frac1p\log\bigl(M(s,p)/(N(s)^{p/2}F(p))\bigr)$, this
  is $G(s,p)<G_0(p)$.

::: {#lem:S .lemmaL}
**Lemma 8** (strict separation). *For $s\in[0,\frac12]$:*

- *if $1<p<p^*$, then $G(s,p)>G_h(p)$ for $s<\frac12$;*

- *if $p=p^*$, then $G(s,p^*)>G_0(p^*)=G_h(p^*)$ for $0<s<\frac12$;*

- *if $p^*<p<2$, then $G(s,p)>G_0(p)$ for $s>0$;*

- *if $p>2$, then $G(s,p)<G_0(p)$ for $s>0$.*
:::

::: proof
*Proof.* *Strict local consequences.* By (F1), (F2) and Taylor's formula
with integral remainder, (F3) gives:

- on L1, for $\frac9{20}\le s<\frac12$ and $1<p\le\frac74$:
  $G(s,p)-G_h(p)=\int_s^{1/2}(t-s)\,\partial_s^2G(t,p)\,dt\ge\frac3{125}(\frac12-s)^2>0$;

- on L2, for $0<s\le\frac1{20}$ and $\frac85\le p\le\frac{39}{20}$:
  $G(s,p)-G_0(p)=\int_0^s(s-t)\,\partial_s^2G(t,p)\,dt\ge\frac{29}{4000}s^2>0$;

- on L3, for $0<s\le\frac1{20}$ and $\frac{41}{20}\le p\le7$:
  $G(s,p)-G_0(p)\le-\frac{31}{4000}s^2<0$.

Note
$\frac85<p^*_{\rm hi}-\eta<p^*_{\rm lo}+\eta<\frac74<\frac{39}{20}$.

*(a) $1<p<p^*$.* If $p\le p^*_{\rm hi}-\eta$, B1 gives $G>G_h$ for
$s\le\frac9{20}$ and L1 for $\frac9{20}\le s<\frac12$. If
$p\in[p^*_{\rm hi}-\eta,p^*)$: at $s=0$, $G=G_0>G_h$ by (F5); for
$0<s\le\frac1{20}$, L2 gives $G>G_0>G_h$; for
$s\in[\frac1{20},\frac9{20}]$, B2 gives $G>\min\{G_0,G_h\}=G_h$; for
$s\in[\frac9{20},\frac12)$, L1 applies.

*(b) $p=p^*$.* $G_0(p^*)=G_h(p^*)$ by (F5). For $0<s\le\frac1{20}$, L2
gives $G>G_0(p^*)$; for $s\in[\frac1{20},\frac9{20}]$, B2 gives
$G>\min\{G_0,G_h\}(p^*)$; for $s\in[\frac9{20},\frac12)$, L1 gives
$G>G_h(p^*)$.

*(c) $p^*<p<2$.* If $p\in(p^*,p^*_{\rm lo}+\eta]$: for
$0<s\le\frac1{20}$, L2; for $s\in[\frac1{20},\frac9{20}]$, B2 gives
$G>\min\{G_0,G_h\}=G_0$ since $H<0$; for $s\in[\frac9{20},\frac12]$, L1
gives $G\ge G_h>G_0$. If $p\in[p^*_{\rm lo}+\eta,\frac{39}{20}]$: L2 for
$0<s\le\frac1{20}$ and B3 for $s\ge\frac1{20}$. If
$p\in[\frac{39}{20},2)$: (F4).

*(d) $p>2$.* If $p\in(2,\frac{41}{20}]$: (F4). If
$p\in[\frac{41}{20},7]$: L3 for $0<s\le\frac1{20}$ and B4 for
$s\in[\frac1{20},\frac12]$. If $p\ge7$: (F7).

The region boundaries are those of the assembly in
Section [8](#sec:assembly){reference-type="ref"
reference="sec:assembly"}. ◻
:::

::: remark
*Remark 15*. The same bounds give quadratic lower bounds in $s$ near the
extremal parameters; for example $G(s,p)-G_0(p)\ge\frac{29}{4000}s^2$ on
$[0,\frac1{20}]\times[\frac85,\frac{39}{20}]$. Combined with
Theorem [3](#thm:R2){reference-type="ref" reference="thm:R2"}, this is a
stability statement for the reduced parameter $s_0$ only. It is not
stability in law, and it is not claimed as such.
:::

# Qualitative stability in distribution {#sec:stability}

The equality classification has the following compactness consequence.
The convergence results used here are standard; no quantitative
deficit-to-distance rate is asserted.

Let $\mathcal C$ be the laws of mean-zero, variance-one log-concave
random variables. Write $$\nu_L=\operatorname{Law}((E-E')/\sqrt2),\qquad
 \nu_+=\operatorname{Law}(E-1),\qquad
 \nu_-=\operatorname{Law}(1-E),$$ and define the variance-normalised
equality set $$\mathcal A_p=
 \begin{cases}
 \{\nu_L\},&1<p<p^*,\\
 \{\nu_L,\nu_+,\nu_-\},&p=p^*,\\
 \{\nu_+,\nu_-\},&p^*<p<2\text{ or }p>2.
 \end{cases}$$ For $X$ with law $\mu\in\mathcal C$, put $$\Delta_p(\mu)=
 \begin{cases}
 \log(\|X\|_p/m_p),&1<p<2,\\
 \log(M_p/\|X\|_p),&p>2.
 \end{cases}$$ This is nonnegative by
Theorem [1](#thm:main){reference-type="ref" reference="thm:main"}. For
$1\le q<\infty$, use the conventions
$$d_{\rm TV}(\mu,\nu)=\sup_{B\text{ Borel}}|\mu(B)-\nu(B)|,
 \qquad
 W_q(\mu,\nu)=\inf_{(X,Y)}(\mathbb E|X-Y|^q)^{1/q},$$ where the infimum
is over couplings with laws $\mu,\nu$, and set
$d_q(\mu,\mathcal A)=\min_{\nu\in\mathcal A}
[d_{\rm TV}(\mu,\nu)+W_q(\mu,\nu)]$ for finite nonempty $\mathcal A$.

::: {#cor:distributional-stability .corollary}
**Corollary 16**. *Fix $p>1$, $p\ne2$, and $1\le q<\infty$. There is a
finite, nondecreasing function $\omega_{p,q}:[0,\infty)\to[0,\infty)$
with $\omega_{p,q}(t)\to0$ as $t\downarrow0$ such that
$$d_q(\mu,\mathcal A_p)\le\omega_{p,q}(\Delta_p(\mu))
 \qquad(\mu\in\mathcal C).$$ Moreover, if $K$ is a nonempty compact
subset of one of $(1,p^*)$, $(p^*,2)$ or $(2,\infty)$, a single function
$\omega_{K,q}$ has this property simultaneously for every $p\in K$.*
:::

::: proof
*Proof.* We first record a uniform tail bound. For $X$ with law in
$\mathcal C$, let $S(t)=\Pr(X\ge t)$. This function is log-concave:
apply preservation under marginalisation to
$f(t+u)\mathbf1_{\{u\ge0\}}$, where $f$ is the density of $X$; see [@SW
Theorem 3.3]. Chebyshev's inequality gives $S(2)\le1/4$ and
$S(-2)\ge3/4$. If $S(2)>0$, the decreasing secant slopes of $\log S$
imply, for $t\ge2$,
$$S(t)\le S(2)\left(\frac{S(2)}{S(-2)}\right)^{(t-2)/4}
 \le\frac14\,3^{-(t-2)/4}.$$ The bound is immediate if $S(2)=0$, and
also holds where $S(t)=0$. Applying it to $-X$ gives $$\begin{equation}
\label{eq:stability-uniform-tail}
 \sup_{\mu\in\mathcal C}\Pr_\mu(|X|\ge t)
 \le\frac12\,3^{-(t-2)/4}\qquad(t\ge2).
\end{equation}$$ In particular every fixed absolute moment is uniformly
integrable on $\mathcal C$.

The class is tight, so any sequence has a weakly convergent subsequence.
Log-concavity of probability measures is preserved under weak limits
[@SW Proposition 3.6]. By
[\[eq:stability-uniform-tail\]](#eq:stability-uniform-tail){reference-type="eqref"
reference="eq:stability-uniform-tail"}, the limit retains mean zero and
variance one. It is therefore nondegenerate and has a log-concave
density [@SW Theorem 2.7], so lies in $\mathcal C$. All fixed absolute
moments converge. Convergence also holds in $W_q$: couple by a common
uniform quantile, which gives almost sure convergence, and use the
uniform higher moments to obtain uniform integrability of the $q$th
power of the coupling difference. Finally, weak convergence to this
nondegenerate log-concave density implies convergence in total variation
by [@CS Proposition 2]. Its hypothesis that the limit has a density has
just been verified; no continuity of that density at its support
endpoints is required.

Suppose $\Delta_p(\mu_n)\to0$. Along any convergent subsequence as
above, moment continuity gives $\Delta_p(\mu)=0$.
Theorem [2](#thm:U){reference-type="ref" reference="thm:U"} then gives
$\mu\in\mathcal A_p$, while convergence in both distances gives
$d_q(\mu_n,\mathcal A_p)\to0$ along this subsequence. A subsequence
bounded away from $\mathcal A_p$ is consequently impossible. One may
therefore take $$\omega_{p,q}(t)=\sup\{d_q(\mu,\mathcal A_p):
                 \mu\in\mathcal C,\ \Delta_p(\mu)\le t\}.$$ This
function is finite by
[\[eq:stability-uniform-tail\]](#eq:stability-uniform-tail){reference-type="eqref"
reference="eq:stability-uniform-tail"}, and the preceding argument
proves its asserted limit (and $\omega_{p,q}(0)=0$).

For the uniform assertion, the set $\mathcal A_p$ is constant on $K$. If
uniformity failed, choose $p_n\in K$ and $\mu_n\in\mathcal C$ with
$\Delta_{p_n}(\mu_n)\to0$ but distances bounded away from zero. Pass to
$p_n\to p\in K$ and a weakly convergent subsequence of $\mu_n$. The tail
bound implies joint continuity of $\mathbb E_{\mu_n}|X|^{p_n}$: on
bounded intervals the integrands converge uniformly, and the tails are
controlled by any fixed moment order larger than $\max K$. The endpoint
constants are continuous in $p$. Hence $\Delta_p(\mu)=0$, and the same
contradiction applies. Taking the supremum also over $p\in K$ defines
$\omega_{K,q}$. ◻
:::

#### Limits of exponent-uniformity.

The fixed-exponent statement includes $p=p^*$, with all three
variance-normalised laws in $\mathcal A_{p^*}$. It cannot, however, be
made uniform through $p^*$ with the exact exponent-dependent equality
sets: if $p_n\uparrow p^*$ with $p_n<p^*$, then
$\Delta_{p_n}(\nu_+)\to0$, but $d_{\rm TV}(\nu_+,\nu_L)>0$ is fixed. Nor
is there such a uniform modulus as $p\to2$. For the standard Gaussian
law $\gamma$, $\Delta_p(\gamma)\to0$ as $p\to2$, whereas
$\min_{\nu\in\{\nu_+,\nu_-\}}d_{\rm TV}(\gamma,\nu)>0$. At $p=2$ every
law in $\mathcal C$ is an equality law. Thus the corollary concerns
convergence to the equality-law set, not selection of a unique branch at
the transition, and provides no explicit rate in the moment deficit.

# Supplementary material and reproducibility {#supplementary-material-and-reproducibility .unnumbered}

The bundle contains four supplementary notes with the complete closed
forms, derivative identities, enclosure lemmas and box-level details
used by the certificates: `supplement/bulk.md`
(Propositions [13](#prop:crossing){reference-type="ref"
reference="prop:crossing"}--[14](#prop:bulk){reference-type="ref"
reference="prop:bulk"}), `supplement/local.md`
(Proposition [11](#prop:local){reference-type="ref"
reference="prop:local"}), `supplement/strip.md`
(Proposition [12](#prop:strip){reference-type="ref"
reference="prop:strip"}) and `supplement/largep.md`
(Proposition [6](#prop:largep){reference-type="ref"
reference="prop:largep"}). It also contains the certifiers, their
outputs, independent floating-point checks and a one-command replay
script (`code/reproduce.sh`). All certificates use python-flint 0.9.0
(Arb) [@Johansson] with exact rational region boundaries, outward
rounding, and failure on any undecided comparison.

#### Rigour of the library functions.

Arb's elementary and special functions used here (powers, $\exp$,
$\log$, $\Gamma$, $\psi$, the upper incomplete gamma function,
${}_1F_1$, ${}_2F_2$ and Tricomi $U$) return, for ball inputs, a ball
containing $f(x)$ for every $x$ in the input balls, or a non-finite ball
when this cannot be guaranteed (for instance when an input ball meets a
pole or a branch cut). Every certifier treats a non-finite output, like
an undecided comparison, as a failure of the box. Integrals in the strip
certificate use Arb's rigorous integration [@JohanssonInt]: the
integrand callback returns enclosures on complex balls and signals
non-analyticity (logarithms are evaluated with the analytic flag, so a
ball touching the branch cut yields a non-finite value and the box
fails). When an external parameter $s$ enters as a ball, every callback
value encloses the integrand for all $s$ in that ball, so the resulting
enclosure holds uniformly in $s$. Where monotonicity in a parameter is
used for an exponent that may be negative (e.g. $\mathbb E(d+sR)^q$ with
$q=p-2<0$ in L2), the certifier takes the hull of all corner
evaluations, which is valid for either direction of monotonicity.

#### Proof-to-code map.

The domains and exact sign statements are given in the cited
propositions; the following table identifies their executable
counterparts. Paths in the last two columns are relative to `code/` and
`certificates/`. All four supplements are included in the same public
archive.

::: center
  Obligation                                                                              Certifier                    Output
  --------------------------------------------------------------------------------------- ---------------------------- ---------------------------
  Bulk B1--B4 (Prop. [14](#prop:bulk){reference-type="ref" reference="prop:bulk"})        `bulk/bulk_cert.py`          `bulk/cert_B*.json`
  Crossing (Prop. [13](#prop:crossing){reference-type="ref" reference="prop:crossing"})   `bulk/crossing_cert.py`      `bulk/cert_X.json`
  Local L1--L3 (Prop. [11](#prop:local){reference-type="ref" reference="prop:local"})     `local/local_cert.py`        `local/cert_L*.json`
  Strip S1--S2 (Prop. [12](#prop:strip){reference-type="ref" reference="prop:strip"})     `strip/strip_cert.py`        `strip/S*_cert.json`
  Large $p$ (Prop. [6](#prop:largep){reference-type="ref" reference="prop:largep"})       `largep/largep_certify.py`   `largep/largep_cert.json`
:::

The large-$p$ conclusion is a written universal argument; the executable
checks its rational constants and optional diagnostics. Numerical
samples do not replace any universal claim. In particular the strip
helper `G_direct` is a truncated diagnostic, not a certified enclosure
of the full function; production S1/S2 use the separately bounded
routines.

#### Status and limitations.

This is an unrefereed computer-assisted theorem candidate. Ball
arithmetic checks the stated regional inequalities, not the complete
written argument in a proof assistant. No unaffiliated reproduction or
journal peer review is asserted. The stability corollary is qualitative,
and the general sharp two-exponent moment problem remains open here. No
arithmetic classification of $p^*$, including irrationality, is claimed.
The public archive is identified by DOI
[10.5281/zenodo.23405173](https://doi.org/10.5281/zenodo.23405173).

::: thebibliography
99 S. Brazitikos, C. Pandis. Sharp inequalities for symmetric
polynomials, Hunter's conjecture, and moments of exponential random
variables. arXiv:2512.12254 (2025). S. Brazitikos, C. Tang, T. Tkocz.
Moments of sums of exponentials, beyond CHS. arXiv:2602.03058 (2026).
Y. Eitan. The centered convex body whose marginals have the heaviest
tails. *Studia Math.* 274 (2024), no. 3, 201--215,
doi:10.4064/sm211027-20-11; arXiv:2110.14382. A. Eskenazis, P. Nayar,
T. Tkocz. Sharp comparison of moments and the log-concave moment
problem. *Adv. Math.* 334 (2018) 389--416. F. Johansson. Arb: efficient
arbitrary-precision midpoint-radius interval arithmetic. *IEEE
Trans. Comput.* 66 (2017) 1281--1292. F. Johansson. Numerical
integration in arbitrary-precision ball arithmetic. In *Mathematical
Software -- ICMS 2018*, LNCS 10931, Springer (2018) 255--263.
J. Melbourne, M. Roysdon, C. Tang, T. Tkocz. From simplex slicing to
sharp reverse Hölder inequalities. *J. London Math. Soc.* (2) 113
(2026), no. 3, e70461; arXiv:2505.00944v2. G. Pólya, G. Szegő. *Problems
and Theorems in Analysis II*. Springer, Berlin (1976). D. Murawski.
Comparing moments of real log-concave random variables.
arXiv:2211.05210. A. Saumard, J. A. Wellner. Log-concavity and strong
log-concavity: A review. *Statist. Surv.* 8 (2014), 45--114,
doi:[10.1214/14-SS107](https://doi.org/10.1214/14-SS107). M. Cule,
R. Samworth. Theoretical properties of the log-concave maximum
likelihood estimator of a multidimensional density. *Electron. J.
Stat.* 4 (2010), 254--270,
doi:[10.1214/09-EJS505](https://doi.org/10.1214/09-EJS505).
:::
