import numpy as np, sys
sys.path.insert(0,'.')
from gi import *
dt=np.dtype([('u','<i8'),('v','<i8'),('N','<i8'),('re','<f8'),('im','<f8'),('kl','<i4'),('ki','<i4'),('om','<i4')])
A=np.fromfile('sqf2_2e6.bin',dtype=dt)
A=A[A['N']<=3000]
P=sorted([tuple(map(int,l.split()[:2])) for l in open('primes_2e4.txt')],key=gnorm)
bad=0;n=0
for r in A:
    c=(int(r['u']),int(r['v']))
    if c==(1,0): continue
    fs,_=factor(c,P)
    kl=qsym((1,1),fs); ki=qsym((0,1),fs)
    n+=1
    if kl!=r['kl'] or ki!=r['ki'] or len(fs)!=r['om']: bad+=1; print(c,kl,ki,r)
print("checked",n,"bad",bad)
