#!/usr/bin/env python3
"""Fail-closed regression tests. Same behavior under Python -O; no bare asserts."""
import copy
import json
from pathlib import Path
import sys
import unittest
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from flint import arb, fmpq, ctx
from certificate_io import interval, point, positive_int, rational
from replay import semantic
from strip.strip_cert import adaptive, verify_cover
from bulk.crossing_cert import sign_cert, deriv_cert
from bulk.bulk_cert import PInt, box_enclosures
from local.local_cert import certify, parts_box
from largep.largep_certify import part_b


class CertificateTests(unittest.TestCase):
    def test_exact_roundtrip_and_containment(self):
        ctx.prec = 128
        for x in [arb(0),arb(1),arb('-1 +/- 0.01'),arb('1e-100 +/- 1e-110'),arb('1e100 +/- 1e90'),arb('1.682457208455160698 +/- 1e-30')]:
            lo,hi = interval(x)
            self.assertLessEqual(Q(lo),Q(hi))
            self.assertTrue(arb(fmpq(lo)) <= x.lower())
            self.assertTrue(arb(fmpq(hi)) >= x.upper())
            self.assertEqual(point(x.lower()),lo)
        x=arb('1.682457208455160698 +/- 1e-30')
        self.assertEqual(x.lower().str(20,radius=False),x.upper().str(20,radius=False))
        self.assertLess(*map(Q,interval(x)))

    def test_nonfinite_rejected(self):
        for text in ('nan','inf','-inf'):
            with self.assertRaises(ValueError): interval(arb(text))
        with self.assertRaises(ValueError): point(arb('1 +/- 0.1'))

    def test_invalid_domains(self):
        for value in (True,False,0,-1,1.5,float('nan')):
            with self.assertRaises(ValueError): positive_int(value,'count')
            with self.assertRaises(ValueError): part_b([Q(7)],value)
        for value in (True,0.1,float('nan'),float('inf')):
            with self.assertRaises(ValueError): rational(value)
        for a,b in [(2,1),(1,1),(0,1),(1,3)]:
            with self.assertRaises(Exception): sign_cert(a,b,1)
            with self.assertRaises(Exception): deriv_cert(a,b)
        with self.assertRaises(Exception): sign_cert(Q(1),Q('3/2'),True)
        with self.assertRaises(ValueError): certify('L1',max_boxes=True)
        with self.assertRaises(ValueError): certify('L1',threshold=True)
        with self.assertRaises(Exception): PInt(fmpq(0),fmpq(1))
        with self.assertRaises(Exception): PInt(True,2)
        with self.assertRaises(Exception): box_enclosures(fmpq(0),fmpq(1),fmpq(1),fmpq(2),('0',))
        with self.assertRaises(Exception): parts_box(0,Q('1/20'),1,2,'small_s')

    def test_strip_cover_failclosed(self):
        s=[Q(0),Q('1/20')]; q=[Q('-1/20'),Q('1/20')]
        leaf=dict(s=list(map(str,s)),q=list(map(str,q)),lower='-2',upper='-1')
        self.assertTrue(verify_cover([leaf],s,q))
        for leaves in ([],[leaf,leaf]):
            with self.assertRaises(Exception): verify_cover(leaves,s,q)
        with self.assertRaises(Exception): verify_cover([], [Q(0),Q(0)],q)
        for count in (True,0,-1):
            with self.assertRaises(ValueError): adaptive('S2',None,s,q,count,1)
        with self.assertRaises(Exception): adaptive('S2',None,s[::-1],q,1,1)
        with self.assertRaises(Exception): adaptive('S2',lambda *a:(arb('nan'),{}),s,q,1,1)
        with self.assertRaises(Exception): adaptive('S2',lambda *a:(arb(1),{}),s,q,1,1)

    def test_semantic_mutations(self):
        c=dict(status='CERTIFIED',lemma='S2',s_range=['0','1/20'],p_range=['39/20','41/20'],n_leaves=1,
               max_upper_bound='-1',leaves=[dict(s=['0','1/20'],q=['-1/20','1/20'],lower='-2',upper='-1')])
        semantic('S2',c)
        mutations=[lambda d:d.update(status='FAILED'),lambda d:d.update(leaves=[]),lambda d:d.update(n_leaves=0),
                   lambda d:d['leaves'][0].update(upper='0'),lambda d:d['leaves'][0].update(lower='1'),
                   lambda d:d['leaves'].append(copy.deepcopy(d['leaves'][0])),lambda d:d.update(max_upper_bound='-2'),
                   lambda d:d.update(p_range=['39/20','2'])]
        for mutation in mutations:
            d=copy.deepcopy(c);mutation(d)
            with self.assertRaises(Exception):semantic('S2',d)
        with self.assertRaises(Exception): semantic('X',dict(status='CERTIFIED',pstar_refined_bracket=['1.6824572084551606980']*2))
        with self.assertRaises(Exception): semantic('unknown',dict(status='PASS'))

    def test_actual_certificates_if_supplied(self):
        # Optional fresh replay directory; failure to find any certificate is an error.
        if not REPLAY: return
        from replay import JOBS
        for name,_,_ in JOBS:
            c=json.loads((Path(REPLAY)/(name+'.json')).read_text());semantic(name,c)
            broken=copy.deepcopy(c);broken['status']='FAIL'
            with self.assertRaises(Exception):semantic(name,broken)
            broken=copy.deepcopy(c)
            if name.startswith('B'): broken['min_certified_margin_lower_bound']='0'
            elif name=='X': broken['pstar_refined_bracket']=[c['pstar_refined_bracket'][0]]*2
            elif name in ('L1','L2','L3'): broken['result']['threshold']='100'
            elif name in ('S1','S2'): broken['leaves'][0]['upper']='0'
            else: broken['part_A_constants']['K_lower_at_p=7']='1000000'
            with self.assertRaises(Exception):semantic(name,broken)


REPLAY = None
if __name__ == '__main__':
    if len(sys.argv)>1 and not sys.argv[1].startswith('-'):
        REPLAY=sys.argv.pop(1)
    unittest.main(verbosity=2)
