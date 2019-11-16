#!/usr/bin/python

from mpmath import *

r=.15/2 # radius of cylinder (m)
l=400.e-9/2 # half-height of cylinder (m)
h=.029 # height of scan above plane (m)
Deltax=.2 # range of scan (m)
dx=.001 # step size for scan (m)

def squiggle_p(z):
    return z+l

def squiggle_m(z):
    return z-l

def alpha_p(rho,z):
    return 1./sqrt(squiggle_p(z)**2+(rho+r)**2)

def alpha_m(rho,z):
    return 1./sqrt(squiggle_m(z)**2+(rho+r)**2)

def beta_p(rho,z):
    return squiggle_p(z)*alpha_p(rho,z)

def beta_m(rho,z):
    return squiggle_m(z)*alpha_m(rho,z)

def gamma(rho):
    return (rho-r)/(rho+r)

def k_p_2(rho,z):
    return (squiggle_p(z)**2+(rho-r)**2)/(squiggle_p(z)**2+(rho+r)**2)

def k_m_2(rho,z):
    return (squiggle_m(z)**2+(rho-r)**2)/(squiggle_m(z)**2+(rho+r)**2)

def k_p(rho,z):
    return sqrt(k_p_2(rho,z))

def k_m(rho,z):
    return sqrt(k_m_2(rho,z))

def curlyk(k):
    #return ellipk(sqrt(1-k**2))
    #See difference in definition between Eq. (6) of paper and sympy documentation
    #https://docs.sympy.org/0.7.1/modules/mpmath/functions/elliptic.html#ellipk
    return ellipk(1-k**2)

def curlye(k):
    #return ellipe(sqrt(1-k**2))
    return ellipe(1-k**2)

def curlyp(gamma,k):
    #return ellippi(1-gamma**2,sqrt(1-k**2))
    return ellippi(1-gamma**2,1-k**2)

def p1(k):
    return curlyk(k)-2*(curlyk(k)-curlye(k))/(1-k**2)

def p3(k,rho):
    return (curlyk(k)-curlye(k))/(1-k**2)-(gamma(rho)**2/(1-gamma(rho)**2))*(curlyp(gamma(rho),k)-curlyk(k))

def p4(k,rho):
    return (gamma(rho)/(1-gamma(rho)**2))*(curlyp(gamma(rho),k)-curlyk(k))+(gamma(rho)/(1-gamma(rho)**2))*(gamma(rho)**2*curlyp(gamma(rho),k)-curlyk(k))-p1(k)



def phi_t(rho,phi,z):
    return (r*cos(phi)/pi)*(beta_p(rho,z)*p3(k_p(rho,z),rho)-beta_m(rho,z)*p3(k_m(rho,z),rho))

def hrho_t(rho,phi,z):
    return (r*cos(phi)/(2*pi*rho))*(beta_p(rho,z)*p4(k_p(rho,z),rho)-beta_m(rho,z)*p4(k_m(rho,z),rho))

def hphi_t(rho,phi,z):
    return (r*sin(phi)/(pi*rho))*(beta_p(rho,z)*p3(k_p(rho,z),rho)-beta_m(rho,z)*p3(k_m(rho,z),rho))

def hz_t(rho,phi,z):
    return (r*cos(phi)/(pi))*(alpha_p(rho,z)*p1(k_p(rho,z))-alpha_m(rho,z)*p1(k_m(rho,z)))

def hx_cart_t(x,y,z):
    rho=sqrt(x**2+y**2)
    phi=atan2(y,x)
    return hrho_t(rho,phi,z)*cos(phi)-hphi_t(rho,phi,z)*sin(phi)

def hy_cart_t(x,y,z):
    rho=sqrt(x**2+y**2)
    phi=atan2(y,x)
    return hrho_t(rho,phi,z)*sin(phi)+hphi_t(rho,phi,z)*cos(phi)

def hz_cart_t(x,y,z):
    rho=sqrt(x**2+y**2)
    phi=atan2(y,x)
    return hz_t(rho,phi,z)

import matplotlib.pyplot as plt 
import numpy as np 
  
x=np.arange(-Deltax,Deltax,dx)
y=0
z=h

xnonzero=[xp for xp in x if abs(xp)>1.e-10]
hxes=[hx_cart_t(xp,y,z) for xp in xnonzero]
hyes=[hy_cart_t(xp,y,z) for xp in xnonzero]
hzes=[hz_cart_t(xp,y,z) for xp in xnonzero]

plt.plot(xnonzero,hxes,label='Hx/Mx')
plt.plot(xnonzero,hyes,label='Hy/Mx')
plt.plot(xnonzero,hzes,label='Hz/Mx')

plt.xlabel('$x$ (m)')
plt.ylabel('$H_{measured}/M_{disk}$ (dimensionless)')
plt.legend()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,-5))
plt.show()


x=0
y=np.arange(-Deltax,Deltax,dx)
z=h

ynonzero=[yp for yp in y if abs(yp)>1.e-10]
hxes=[hx_cart_t(x,yp,z) for yp in ynonzero]
hyes=[hy_cart_t(x,yp,z) for yp in ynonzero]
hzes=[hz_cart_t(x,yp,z) for yp in ynonzero]

plt.plot(ynonzero,hxes,label='Hx/Mx')
plt.plot(ynonzero,hyes,label='Hy/Mx')
#plt.plot(ynonzero,hzes,label='Hz/Mx')

plt.xlabel('$y$ (m)')
plt.ylabel('$H_{measured}/M_{disk}$ (dimensionless)')
plt.legend()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,-5))
plt.show()
