"""Exploratory exponent model, NOT an analytic proof or exact certificate.

Independent candidate: use a fourth moment for the numerator on Re(w)=1/2.
Only the original kappa=3/4 is used. Contour choice is by zero bin delta,
outside the supremum in the prime-amplitude variable x.
"""
import numpy as np

delta=np.linspace(1/50, 3/4, 2401)[:,None]
x=np.linspace(0,1/2,201)[None,:]
c=2/9
alpha=5/6
D=3-(1+4*c)*x
P=(2-4*c*x)*(1-x)
J=(alpha-delta)*D+delta*P
R=1-delta+(alpha-delta)*delta*P/(2*J)

def endpoint(ell,b,details=False,weighted=False):
    theta=11/12-ell/4
    ly=(1-ell+b)/2
    h=(1+3*ell+b)/2
    old=(1+delta)*(1-ly)/2-ell/2+x*delta*ell+h*(R+delta/2-1/6)-theta
    central=old+h*((1-R)/4-delta/2)+ly*delta/2
    capacity=max(0,(2*ly/h-1)/(6*(3/4)))
    if weighted:
        central=central-h*x*delta*capacity/2
    old_worst=old.max(axis=1)
    central_worst=central.max(axis=1)
    combined=np.minimum(old_worst, central_worst)
    j=int(combined.argmax())
    if details:
        return dict(theta=theta,ell=ell,b=b,h=h,ly=ly,
                    worst=float(combined[j]),critical_delta=float(delta[j,0]),
                    old_worst=float(old_worst[j]),central_worst=float(central_worst[j]),
                    old_x=float(x[0,old[j].argmax()]),central_x=float(x[0,central[j].argmax()]),
                    central_bins=int((central_worst<old_worst).sum()))
    return float(combined[j])

if __name__=='__main__':
    import sys
    weighted='--weighted' in sys.argv
    for ell in [.16684,.167,.168,.169,.17,.175,.18,.19,.2]:
        best=min((endpoint(ell,float(b),weighted=weighted),float(b)) for b in np.linspace(.005,.2,131))
        print(endpoint(ell,best[1],True,weighted=weighted),flush=True)
