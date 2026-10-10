from gm2 import best
import numpy as np
# saturated marked (cM=1) + plain cP scan; + alpha scan
for cP in [0.25,0.3,0.33,0.36,0.38,0.4]:
    (w,arg),b=best(cM=1,cap82=False,cP=cP); print('mk sat, cP',cP,f"{w:+.5f}",arg,b)
for al in [0.8,0.78,0.76,0.75,0.74]:
    (w,arg),b=best(cM=1,cap82=False,alpha=al); print('mk sat, alpha',al,f"{w:+.5f}",arg,b)
# minimal mark capacity with cP=0.4
for cM in [1.25,1.3,1.35]:
    (w,arg),b=best(cM=cM,cap82=False,cP=0.4); print('cP=.4 cM',cM,f"{w:+.5f}",arg,b)
for wM in [0.12,0.14,0.16,0.18]:
    (w,arg),b=best(wM=wM,cap82=False,cP=0.4); print('cP=.4 wM',wM,f"{w:+.5f}",arg,b)
