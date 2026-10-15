# Random test: nontrivial combos of {1,x,|x|^a,|x|^r sgn x,|x|^b} have <=4 sign changes on a fine grid;
# and V(c)+V(c*) <= n combinatorial lemma.
import numpy as np, itertools
def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)
rng=np.random.default_rng(1)
def V(c):
    c=[x for x in c if x!=0]; return sum(1 for i in range(len(c)-1) if c[i]*c[i+1]<0)
for n in range(1,9):
    for c in itertools.product([-1,0,1],repeat=n+1):
        if any(c): require(V(c)+V([(-1)**i*x for i,x in enumerate(c)])<=n-(1 if c[0]==0 else 0), (n,c))
x=np.concatenate([-np.logspace(-6,2,40000)[::-1],[0.0],np.logspace(-6,2,40000)])
worst=0
for t in range(20000):
    a,r,b=np.sort(rng.uniform(1.0001,8,3))
    U=np.vstack([np.ones_like(x),x,abs(x)**a,abs(x)**r*np.sign(x),abs(x)**b])
    # construct combos with many zeros: interpolate zeros at 4 random points then perturb
    z=np.sort(rng.uniform(-3,3,4)); M=np.vstack([[1,t_,abs(t_)**a,abs(t_)**r*np.sign(t_),abs(t_)**b] for t_ in z])
    c=np.linalg.svd(M)[2][-1]+rng.normal(0,1e-3,5)*(t%2)
    y=c@U; s=np.sign(y); s=s[s!=0]; k=int(np.sum(s[1:]!=s[:-1])); worst=max(worst,k)
    require(k<=4,(a,r,b,c,k))
print("combinatorial lemma OK; max sign changes found:",worst)
