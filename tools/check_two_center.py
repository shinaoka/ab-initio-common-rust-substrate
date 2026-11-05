# SPDX-License-Identifier: MIT
# Numerical check of the two-center overlap/kinetic formulas in
# src/modules/two-center-overlap-kinetic.md, using Gaussian-type radial
# functions R(r) = r^l exp(-alpha r^2) and a direct 3D grid integration.
# Run: uv run --with scipy --with sympy --with numpy python tools/check_two_center.py

import numpy as np
from scipy.special import spherical_jn, sph_harm_y
from scipy.integrate import quad
from sympy.physics.wigner import gaunt
# complex Y_lm(theta,phi) with Condon-Shortley; scipy sph_harm(m,l,azimuth,polar)
def Y(l,m,th,ph): return sph_harm_y(l,m,th,ph)
def G(l1,m1,l2,m2,L,M):
    # ∫ conj(Y_l1m1) Y_l2m2 Y_LM ; conj(Y_lm)=(-1)^m Y_l,-m
    return float((-1)**m1*gaunt(l1,l2,L,-m1,m2,M))
a,b=0.8,1.3
def Rt(l,al,k): return quad(lambda r: r*r*spherical_jn(l,k*r)*r**l*np.exp(-al*r*r),0,12,limit=200)[0]
def S_formula(l,m,lp,mp,tau):
    t=np.linalg.norm(tau); th=np.arccos(tau[2]/t); ph=np.arctan2(tau[1],tau[0])
    s=0
    for L in range(abs(l-lp),l+lp+1):
        M=m-mp
        if abs(M)>L: continue
        g=G(l,m,lp,mp,L,M)
        if g==0: continue
        I=quad(lambda k: k*k*spherical_jn(L,k*t)*Rt(l,a,k)*Rt(lp,b,k),0,25,limit=200)[0]
        s+=8*(1j)**(l-lp-L)*g*np.conj(Y(L,M,th,ph))*I
    return s
def S_direct(l,m,lp,mp,tau,n=70,Lbox=7.0):
    x=np.linspace(-Lbox,Lbox,n); h=x[1]-x[0]
    X,Yg,Z=np.meshgrid(x,x,x,indexing='ij')
    def phi(l,m,al,X,Yg,Z):
        r=np.sqrt(X*X+Yg*Yg+Z*Z)+1e-300; th=np.arccos(np.clip(Z/r,-1,1)); ph=np.arctan2(Yg,X)
        return r**l*np.exp(-al*r*r)*Y(l,m,th,ph)
    return np.sum(np.conj(phi(l,m,a,X,Yg,Z))*phi(lp,mp,b,X-tau[0],Yg-tau[1],Z-tau[2]))*h**3
tau=np.array([0.4,-0.3,0.9])
for (l,m,lp,mp) in [(0,0,0,0),(1,0,0,0),(0,0,1,1),(1,1,1,-1),(2,1,1,0),(1,-1,2,-2)]:
    print((l,m,lp,mp), np.round(S_formula(l,m,lp,mp,tau),6), np.round(S_direct(l,m,lp,mp,tau),6))
t=np.linalg.norm(tau); mu=a*b/(a+b)
I2=quad(lambda k: k**4*spherical_jn(0,k*t)*Rt(0,a,k)*Rt(0,b,k),0,25,limit=200)[0]
T=4*G(0,0,0,0,0,0)*np.conj(Y(0,0,0,0))*I2
S=S_formula(0,0,0,0,tau)
print("T formula",T.real,"analytic",(mu*(3-2*mu*t*t)*S).real)
