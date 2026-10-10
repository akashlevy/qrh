"""Generate explicit Lean polynomial identities; Lean checks every generated identity.

No proof depends on Python arithmetic, the JSON signs, or native_decide.
The leaf rectangles are merely proposed data; exact ring identities and
positivity in Lean produce the complete proof terms.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

def add(*polys):
    out={}
    for p in polys:
        for k,v in p.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(a,p):return {k:a*v for k,v in p.items() if a*v}
def mul(p,q):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():r[i+k,j+l]=r.get((i+k,j+l),F(0))+a*b
    return {k:v for k,v in r.items() if v}
def const(v):return {(0,0):F(v)}
def rat(v):
    v=F(v)
    if v.denominator==1:return str(v.numerator)
    return f'({v.numerator} / {v.denominator} : ℝ)'

d={(1,0):F(1)};x={(0,1):F(1)};q=mul(d,x)
D=add(const(3),scale(-F(17,9),x))
P=mul(add(const(2),scale(-F(8,9),x)),add(const(1),scale(-1,x)))
J=add(mul(add(const(F(5,6)),scale(-1,d)),D),mul(d,P))
AP=mul(mul(add(const(F(5,6)),scale(-1,d)),d),P)
RJ=add(mul(add(const(1),scale(-1,d)),J),scale(F(1,2),AP))
ly=F(919,2000);h=F(1587,2000);ell=F(167,1000);theta=F(10499,12000)
rest=add(const((1-ly)/2-ell/2-h/6-theta),scale((1-ly)/2+h/2,d),scale(ell,q))
old=add(mul(rest,J),scale(h,RJ))
corr=add(scale(h/4,add(J,scale(-1,RJ))),scale((ly-h)/2,mul(d,J)),scale(-(2*ly-h)/9,mul(q,J)))
new=add(old,corr)
polys=[scale(-1,add(old,scale(F(1,10000),J))),scale(-1,add(new,scale(F(1,10000),J)))]

def decomposition(poly,box):
    da,db,xa,xb=box
    n=max(i for i,j in poly);m=max(j for i,j in poly)
    power={}
    for (i,j),a in poly.items():
        for k in range(i+1):
            for l in range(j+1):
                v=a*comb(i,k)*da**(i-k)*(db-da)**k*comb(j,l)*xa**(j-l)*(xb-xa)**l
                power[k,l]=power.get((k,l),F(0))+v
    terms=[]
    for i in range(n+1):
        for j in range(m+1):
            v=sum(power.get((k,l),F(0))*F(comb(i,k),comb(n,k))*F(comb(j,l),comb(m,l)) for k in range(i+1) for l in range(j+1))
            assert v>=0
            v*=F(comb(n,i)*comb(m,j),(db-da)**n*(xb-xa)**m)
            fs=[rat(v)]
            for factor,pow_ in [(f'(δ - {rat(da)})',i),(f'({rat(db)} - δ)',n-i),(f'(x - {rat(xa)})',j),(f'({rat(xb)} - x)',m-j)]:
                if pow_==1:fs.append(factor)
                elif pow_>1:fs.append(f'{factor} ^ {pow_}')
            terms.append(' * '.join(fs))
    return '\n      + '.join(terms)

base=Path(__file__).resolve().parent
cert=json.loads((base.parent/'weighted_numerator_certificate.json').read_text())
s='''import WeightedQRH.Parameters

namespace WeightedQRH

noncomputable section

set_option maxRecDepth 4096
set_option maxHeartbeats 4000000

/-- D_x in the manuscript. -/
def countD (x : ℝ) : ℝ := 3 - (17 / 9 : ℝ) * x
/-- P_x in the manuscript. -/
def countP (x : ℝ) : ℝ := (2 - (8 / 9 : ℝ) * x) * (1 - x)
/-- The common strictly positive denominator J. -/
def countDen (δ x : ℝ) : ℝ := (5 / 6 - δ) * countD x + δ * countP x
/-- The numerator R J of the selected row count. -/
def countNum (δ x : ℝ) : ℝ := (1 - δ) * countDen δ x +
  (5 / 6 - δ) * δ * countP x / 2

def rowCount (δ x : ℝ) : ℝ := countNum δ x / countDen δ x

def oldRest (δ x : ℝ) : ℝ :=
  (1 + δ) * (1 - ly) / 2 - ell / 2 + δ * x * ell + h * δ / 2 - h / 6 - theta

def oldCert (δ x : ℝ) : ℝ :=
  -((oldRest δ x + endpointMargin) * countDen δ x + h * countNum δ x)

def weightedCert (δ x : ℝ) : ℝ := oldCert δ x -
  (h / 4 * (countDen δ x - countNum δ x) +
   δ * (ly - h) / 2 * countDen δ x -
   δ * x * (2 * ly - h) / (12 * kappa) * countDen δ x)

theorem countDen_pos {δ x : ℝ}
    (hδ0 : 0 ≤ δ) (hδ1 : δ ≤ 3 / 4) (hx0 : 0 ≤ x) (hx1 : x ≤ 1 / 2) :
    0 < countDen δ x := by
  have ha : 0 < (5 / 6 : ℝ) - δ := by linarith
  have hb : 0 < countD x := by unfold countD; linarith
  have hc : 0 ≤ countP x := by
    have hc1 : 0 ≤ (2 : ℝ) - 8 / 9 * x := by linarith
    have hc2 : 0 ≤ (1 : ℝ) - x := by linarith
    exact mul_nonneg hc1 hc2
  exact add_pos_of_pos_of_nonneg (mul_pos ha hb) (mul_nonneg hδ0 hc)

'''
for typ in range(2):
    name='old' if typ==0 else 'weighted'
    for idx,leaf in enumerate(cert['branches'][typ]['leaves']):
        box=list(map(F,leaf['box']));da,db,xa,xb=box
        expr=decomposition(polys[typ],box)
        s+=f'''/-- Exact unnormalised Bernstein identity on certificate leaf {idx}. -/
private theorem {name}_identity_{idx} (δ x : ℝ) :
    {name}Cert δ x =
      {expr} := by
  unfold {name}Cert'''
        if typ:s+=' oldCert'
        s+=' oldRest countNum countDen countD countP\n'
        s+='  norm_num [theta, ell, h, ly, kappa, endpointMargin]\n  ring\n\n'
        s+=f'''theorem {name}_leaf_{idx} {{δ x : ℝ}}
    (hδ0 : {rat(da)} ≤ δ) (hδ1 : δ ≤ {rat(db)})
    (hx0 : {rat(xa)} ≤ x) (hx1 : x ≤ {rat(xb)}) :
    0 ≤ {name}Cert δ x := by
  have ha : 0 ≤ δ - {rat(da)} := by linarith
  have hb : 0 ≤ {rat(db)} - δ := by linarith
  have hc : 0 ≤ x - {rat(xa)} := by linarith
  have hd : 0 ≤ {rat(xb)} - x := by linarith
  rw [{name}_identity_{idx}]
  positivity

'''
s+='''theorem oldCert_nonneg {δ x : ℝ}
    (hδ0 : 0 ≤ δ) (hδ1 : δ ≤ 1 / 3) (hx0 : 0 ≤ x) (hx1 : x ≤ 1 / 2) :
    0 ≤ oldCert δ x := old_leaf_0 hδ0 hδ1 hx0 hx1

theorem weightedCert_nonneg {δ x : ℝ}
    (hδ0 : 1 / 3 ≤ δ) (hδ1 : δ ≤ 3 / 4) (hx0 : 0 ≤ x) (hx1 : x ≤ 1 / 2) :
    0 ≤ weightedCert δ x := by
  by_cases hxa : x ≤ 1 / 4
  · apply weighted_leaf_0 <;> linarith
  · by_cases hda : δ ≤ 13 / 24
    · by_cases hxb : x ≤ 3 / 8
      · apply weighted_leaf_1 <;> linarith
      · by_cases hdb : δ ≤ 7 / 16
        · apply weighted_leaf_2 <;> linarith
        · by_cases hxc : x ≤ 7 / 16
          · apply weighted_leaf_3 <;> linarith
          · by_cases hdc : δ ≤ 47 / 96
            · by_cases hxd : x ≤ 15 / 32
              · apply weighted_leaf_4 <;> linarith
              · by_cases hdd : δ ≤ 89 / 192
                · apply weighted_leaf_5 <;> linarith
                · apply weighted_leaf_6 <;> linarith
            · apply weighted_leaf_7 <;> linarith
    · apply weighted_leaf_8 <;> linarith

end
end WeightedQRH
'''
(base/'WeightedQRH/CertificateBernstein.lean').write_text(s)
print('Generated',len(s),'characters; polynomial degrees',[(max(i for i,j in p),max(j for i,j in p)) for p in polys])
