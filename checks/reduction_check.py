# Numerical sanity check of the comparison lemma (Thm R): for centred log-concave X,
# find s0 with rho(Z_s0)=rho(X), scale, and verify E|X|^b <= E|lam Z_s0|^b; also check
# rho(Z_0) <= rho(X) <= rho(Z_1).
import numpy as np
def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)
from scipy import integrate, optimize, stats
def mom(dens, lo, hi, fn):
    g=lambda x: fn(x)*dens(x); q=lambda u,v: integrate.quad(g,u,v,limit=400,epsabs=1e-13,epsrel=1e-12)[0]
    cuts=sorted(set([lo,hi]+[t for t in (-1.0,0.0,1.0) if lo<t<hi]))
    return sum(q(cuts[i],cuts[i+1]) for i in range(len(cuts)-1))
def zs(s):
    c=1-s; d=1-2*s
    if s==0: return (lambda z: np.exp(z-1)), -np.inf, 1.0
    if s==1: return (lambda z: np.exp(-(z+1))), -1.0, np.inf
    return (lambda z: np.exp(-(z-d)/s) if z>=d else np.exp((z-d)/c)), -np.inf, np.inf
def stats_of(dens,lo,hi,a,r,b):
    A=mom(dens,lo,hi,lambda x:abs(x)**a); Rr=mom(dens,lo,hi,lambda x:abs(x)**r*np.sign(x)); B=mom(dens,lo,hi,lambda x:abs(x)**b)
    return A,Rr,B
def rho(A,Rr,a,r): return Rr/A**(r/a)
tests={
 'gamma2':(lambda x: stats.gamma(2,loc=-2).pdf(x),-2,np.inf),
 'neg_gamma3':(lambda x: stats.gamma(3,loc=-3).pdf(-x),-np.inf,3),
 'uniform':(lambda x:0.5,-1,1),
 'gauss':(lambda x: stats.norm.pdf(x),-np.inf,np.inf),
 'gumbel':(lambda x: stats.gumbel_r(loc=-np.euler_gamma).pdf(x),-np.inf,np.inf),
 'tri_skew':(lambda x: stats.triang(0.2,loc=-0.4,scale=1.2).pdf(x+0.0),-0.4,0.8),
}
# centre tri_skew numerically
m=mom(*tests['tri_skew'],lambda x:x); f0=tests['tri_skew'][0]; tests['tri_skew']=(lambda x,f0=f0,m=m: f0(x+m),-0.4-m,0.8-m)
for (a,b) in [(1.3,2),(1.68,2),(1.9,2),(2,2.5),(2,4),(2,7.5)]:
    r=(a+b)/2
    Z0=stats_of(*zs(0),a,r,b); Z1=stats_of(*zs(1),a,r,b)
    for name,(f,lo,hi) in tests.items():
        require(abs(mom(f,lo,hi,lambda x:x))<1e-8 and abs(mom(f,lo,hi,lambda x:1)-1)<1e-8, name)
        A,Rr,B=stats_of(f,lo,hi,a,r,b); rx=rho(A,Rr,a,r)
        require(rho(*Z0[:2],a,r)-1e-12<=rx<=rho(*Z1[:2],a,r)+1e-12,(name,a,b))
        g=lambda s: rho(*stats_of(*zs(s),a,r,b)[:2],a,r)-rx
        s0=optimize.brentq(g,1e-9,1-1e-9,xtol=1e-13)
        As,Rs,Bs=stats_of(*zs(s0),a,r,b); lam=(A/As)**(1/a)
        gap=lam**b*Bs-B
        print(f"a={a} b={b} {name:10s} s0={s0:.5f} E|Y|^b-E|X|^b={gap:.3e}")
        require(gap>0, (name,a,b,gap))
print("ALL OK")
