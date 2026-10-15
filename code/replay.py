#!/usr/bin/env python3
"""Clean same-code replay; semantic output checks are not independent verification."""
import argparse, hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction as Q
from certificate_io import require, runtime
ROOT = Path(__file__).resolve().parent
JOBS = [(r, 'bulk/bulk_cert.py', [r]) for r in ('B1','B2','B3','B4')]
JOBS += [('X','bulk/crossing_cert.py',[])]
JOBS += [(r,'local/local_cert.py',[r,'--threshold',t]) for r,t in [('L1','6/125'),('L2','29/2000'),('L3','31/2000')]]
JOBS += [(r,'strip/strip_cert.py',['--lemma',r]) for r in ('S1','S2')]
JOBS += [('LP','largep/largep_certify.py',[])]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def semantic(name, c):
    require(c.get('status') in ('PASS','CERTIFIED'), 'certificate failed')
    if name.startswith('B'):
        from bulk.bulk_cert import REGIONS
        s0,s1,p0,p1,refs,sense = REGIONS[name]
        require(c['s_range'] == list(map(str,(s0,s1))) and c['p_range'] == list(map(str,(p0,p1))), 'bulk domain')
        require(c['region'] == name and c['failures'] == 0 and c['leaves_certified'] > 0, 'bulk failure/empty')
        require(Q(c['min_certified_margin_lower_bound']) > 0, 'bulk margin')
    elif name == 'X':
        from bulk.crossing_cert import H
        from flint import arb, fmpq, ctx
        ctx.prec = 128
        lo,hi = map(Q,c['pstar_refined_bracket'])
        require(Q('1682457208/1000000000') <= lo < hi <= Q('1682457209/1000000000'), 'collapsed/outside crossing bracket')
        require(H(arb(fmpq(str(lo)))).lower() > 0 and H(arb(fmpq(str(hi)))).upper() < 0, 'crossing endpoint signs')
        for part in ('X1','X2','X3'):
            require(c[part]['failures'] == 0, 'crossing failures')
        require(c['X1']['interval'] == ['1',c['pstar_lo']] and c['X2']['interval'] == [c['pstar_hi'],'39/20'] and c['X3']['interval'] == [c['pstar_lo'],c['pstar_hi']], 'crossing domains')
        require(all(Q(c[x]['min_certified_margin']) > 0 for x in ('X1','X2')), 'crossing margins')
        require(Q(c['X3']["max_certified_upper_bound_of_H'"]) < 0, 'derivative sign')
    elif name in ('L1','L2','L3'):
        from local.local_cert import REGIONS
        r = c['result']; spec = REGIONS[name]; threshold = Q(dict(L1='6/125',L2='29/2000',L3='31/2000')[name])
        require(r['s_range'] == list(map(str,spec[:2])) and r['p_range'] == list(map(str,spec[2:4])), 'local domain')
        require(r['failures'] == 0 and r['coverage_area_ok'] is True and r['leaves'] > 0, 'local coverage/failures')
        require(Q(r['threshold']) == threshold and Q(r['min_certified_abs_Gss']) > threshold, 'local threshold')
        lo,hi = map(Q,r['worst_leaf'][4])
        require(lo <= hi and (lo > threshold if name != 'L3' else hi < -threshold), 'local worst bound')
    elif name in ('S1','S2'):
        from strip.strip_cert import verify_cover
        s = [Q('1/20'),Q('1/2')] if name == 'S1' else [Q(0),Q('1/20')]
        q = [Q('-1/20'),Q('1/20')]
        require(c['lemma'] == name and list(map(Q,c['s_range'])) == s and list(map(Q,c['p_range'])) == [x+2 for x in q], 'strip domain')
        verify_cover(c['leaves'],s,q)
        require(c['n_leaves'] == len(c['leaves']), 'strip count')
        for leaf in c['leaves']:
            require(Q(leaf['lower']) <= Q(leaf['upper']) < 0, 'strip sign/order')
        require(Q(c['max_upper_bound']) == max(Q(x['upper']) for x in c['leaves']), 'strip maximum')
    elif name == 'LP':
        a = c['part_A_constants']
        require(Q(a['e_upper']) == Q('87/32') and Q(a['kappa_lower(p>=7)']) == Q('1/16'), 'large-p constants')
        require(Q(a['K_lower_at_p=7']) == 7*(45-Q('87/64'))-Q('2001/8192'), 'large-p bound')
        require(len(c['part_B_audit']['rows']) == 9 and all(r['n_tau'] == 64 for r in c['part_B_audit']['rows']), 'audit grid')
        require(all(c['part_C_negative_control'][p].startswith('failed as expected:') for p in ('6','34/5')), 'negative controls')
    else:
        raise ValueError('unknown certificate')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,help='new directory; must not already exist')
    ap.add_argument('--optimized',action='store_true')
    args = ap.parse_args(); versions = runtime(); out = Path(args.out).resolve()
    out.mkdir(parents=False,exist_ok=False)
    manifest = {str(p.relative_to(ROOT.parent)):sha(p) for d in (ROOT,ROOT.parent/'checks') for p in d.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    record = dict(environment=versions,python=sys.version,executable=sys.executable,platform=platform.platform(),optimized=args.optimized,input_sha256=manifest,started=datetime.now(timezone.utc).isoformat(),runs=[])
    (out/'environment.json').write_text(json.dumps(record,indent=2)+'\n')
    all_ok = True
    for name,script,argv in JOBS:
        dest = out/(name+'.json'); command = [sys.executable,'-B']+(['-O'] if args.optimized else [])+[str(ROOT/script),*argv,'--out',str(dest)]
        start = time.monotonic(); stamp = datetime.now(timezone.utc).isoformat()
        with (out/(name+'.stdout')).open('x') as stdout, (out/(name+'.stderr')).open('x') as stderr:
            proc = subprocess.run(command,stdout=stdout,stderr=stderr,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        error = None
        try:
            require(proc.returncode == 0,'nonzero process exit')
            cert = json.loads(dest.read_text())
            require(cert['script_sha256'] == sha(ROOT/script),'script hash mismatch')
            semantic(name,cert)
        except Exception as exc:
            error = str(exc); all_ok = False
        record['runs'].append(dict(name=name,command=command,started=stamp,elapsed_seconds=time.monotonic()-start,returncode=proc.returncode,semantic_error=error,certificate_sha256=sha(dest) if dest.exists() else None,stdout_sha256=sha(out/(name+'.stdout')),stderr_sha256=sha(out/(name+'.stderr'))))
        print(name,'PASS' if error is None else 'FAIL: '+error,flush=True)
    record.update(finished=datetime.now(timezone.utc).isoformat(),status='PASS' if all_ok else 'FAIL')
    (out/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    return 0 if all_ok else 1

if __name__ == '__main__':
    raise SystemExit(main())
