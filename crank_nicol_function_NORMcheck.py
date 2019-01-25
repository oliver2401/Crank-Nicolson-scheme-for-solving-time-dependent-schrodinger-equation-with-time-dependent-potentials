#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 15:39:00 2019

@author: Oliver
"""

import scipy as sp
import numpy as np
from scipy import integrate, sparse, linalg,special
import scipy.sparse.linalg
import pylab as pl
from multiprocessing import Pool
import multiprocessing

def crank_nicol(tf):
    ax=-10. #punto inicial
    bx=10. #punto final
    #ti=0. #tiempo inicial
    #tf=628 #tiempo final
    nx= 500#numero de puntos en el x-grid
    Lx=bx-ax #longitud total
    #T=tf-ti #tiempo total
    delx= Lx/nx #x-step size=Lx/nx
    dt= 1.e-3#T/mt # time step delta t
    mt=int((2*np.pi*tf*0.333)/dt)#mt=int(tf/dt)
    print (' MT = ', mt)
    nux= 1j*dt/(4.0*(delx**2)) #nu q son iguales a mis c_i=b_i 
    
    gridx = sp.zeros(nx) #space grid
    igridx = sp.array(range(nx)) # time grid
    psi = sp.zeros(nx)
    pot = sp.zeros(nx)
    gridx = delx*(igridx - nx/2)
    #x=gridx

###############Parte de dubinko del potencial cuartico#############
    sigma=0.08
    alpha= 2.614
    D=1#75
    #R_ave=2.75 #time dependent part
    R=2.95#R_ave-R_ave*0.1*np.cos()#2.95 # this is the parameter that either is a constant or varies in time, is DAD
    r_e=.9699
    R_e=R-2*r_e #static one
    r_0=(1./alpha)*np.arccosh(0.5*np.exp(alpha*R_e/2.0))#checar esta parte que debe de estar en el loop?? Static
    x=alpha*gridx
    X_e=alpha*R_e
    E_barrier=D*( 1-2*np.exp(-alpha*R_e/2.) )**2
    ##
    # Compute the initial energy and print it 
    
    #xmin=(1/alpha)*np.arccosh(0.5*np.exp(alpha*0.812/2.0)) time dependent
####################################################
    #psi=(2*sp.pi*sigma**2)**(-1/4.) *np.exp( -(gridx+r_0)**2 /(4*sigma**2) )
    #psi /= sp.sqrt(sp.integrate.simps(psi*np.conjugate(psi), x))
    
    psi=(2*sp.pi*sigma**2)**(-0.25) *np.exp( -(x+alpha*r_0)**2 /(4*sigma**2) )
    psi /= sp.sqrt(sp.integrate.simps(psi*np.conjugate(psi), x))



    print(nux)
#######Creacion de la matriz A+######################
    Adiag=sp.empty(nx,dtype=complex) #puede q el dtype estple mal asi q se puede quitar, esta va a ser mi diag 
    Asup=sp.empty(nx,dtype=complex) # igual q arriba con el dtype, esta va a ser mi diagonal superior
    Asub=sp.empty(nx,dtype=complex) #igual q arriba con el dtype, esta va a ser mi diagonal inferior
    Adiag.fill(1+2.0*nux) #llena la diag principal
    Asup.fill(-nux)
    Asub.fill(-nux)
    aplus=sp.sparse.spdiags([Adiag,Asup,Asub],[0,1,-1],nx,nx) #crea la matriz A+
#######Creacion de la matrix A-########
    adiag=sp.empty(nx,dtype=complex)
    asup=sp.empty(nx,dtype=complex)
    asub=sp.empty(nx,dtype=complex)
    adiag.fill(1-2.0*nux)
    asup.fill(nux)
    asub.fill(nux)
    aminus=sp.sparse.spdiags([adiag,asup,asub],[0,1,-1],nx,nx) #crea mi matriz A-
#############################################################

###############parametros de time-dependent eigenfreq con nx=1000,mt=100000
   
    pot=(D*( (1-np.exp(-(X_e/2 + x)))**2 + (1-np.exp(-(X_e/2 -x)))**2   - (1 - 2*np.exp(-X_e)) ))    #pot= (D*( (1-np.exp(-alpha*(R_e/2 + gridx)))**2 + (1-np.exp(-alpha*(R_e/2 -gridx)))**2   - (1 - 2*np.exp(-alpha*R_e)) ))
    #omega=0.333#2*np.pi#time dependent part
    time=[]
    ik = 0
    xave = []
    norm = []
    ts = []
    t1 = []
    t2 = []
    iw1 = 0
    iw2 = 0
    ns = 10
    naxv = int(float(mt)/ns)+1
    t1 = np.zeros(naxv)
    t2 = np.zeros(naxv)
    xav_old = -r_0
    for t in range(mt) :
        #R=R_ave+R_ave*0.1*np.sin(t*dt)#time dependent, si pones omega alenta por 400 veces el movimiento del potencial
        #R_e=R-2*r_e#time dependent
        #pot=pot=(75*( (1-np.exp(-alpha*(R_e/2 + x)))**2 + (1-np.exp(-alpha*(R_e/2 -x)))**2   - (1 - 2*np.exp(-alpha*R_e)) ))
        #pot=(D*( (1-np.exp(-(X_e/2 + x)))**2 + (1-np.exp(-(X_e/2 -x)))**2   - (1 - 2*np.exp(-X_e)) ))
        psi= sp.exp(-1j*dt*pot)*psi
        psi = sp.sparse.linalg.bicg(aplus, aminus*psi)[0]
        #psi/=sp.sqrt(sp.integrate.simps(psi*np.conjugate(psi),x))#,dx=delx))
        psi_area=psi[np.where(x>0)]
        #areaofinterest=sp.integrate.simps(psi_area*np.conjugate(psi_area),x[np.where(x>0)])
        #if areaofinterest >0.5:
        #    time.append(t)
        # Sampling
        if np.mod(t,ns) == 0:           
            ts.append(t*dt)
            xav = sp.integrate.simps(psi*np.conjugate(psi)*x,x)
            if xav < 0 :
                t1[iw1] += dt*ns
            if xav > 0 :
                t2[iw2] += dt*ns
           
            if xav > 0 and xav_old < 0:  # transition 1 --> 2
                iw2  += 1
            if xav < 0 and xav_old > 0:  # transition 2 --> 1
                iw1  += 1  
             
                
            xave.append(xav)
            nn = sp.sqrt(sp.integrate.simps(psi*np.conjugate(psi),x))
            norm.append(nn)    #,dx=delx))  
            xav_old = xav
            
            
        
    first=np.gradient(psi,x)
    AvJcurr=sp.integrate.simps(np.imag(np.conjugate(psi)*first)*psi*np.conjugate(psi),x )#np.sum( np.imag(psi*first)*psi*np.conjugate(psi)*delx)
    #xave=sp.integrate.simps(psi*np.conjugate(psi)*x,x)#np.sum(psi*np.conjugate(psi)*x*delx)
    time=np.array(time)
    norm=np.array(norm)
    t1=np.array(t1)
    t2=np.array(t2)
    return x, ts, psi*np.conjugate(psi),tf,pot,AvJcurr,np.array(xave),time, norm, t1, t2, iw1, iw2

x, ts, prob,tf,pot,AvJcurr,xave,time, norm, t1, t2, iw1, iw2 = crank_nicol(2.)

pl.subplot(2,2,1)
pl.plot(ts,norm,'r-')

pl.subplot(2,2,2)
pl.plot(ts,xave,'g-')
pl.plot(ts,np.zeros(len(ts)),'b--')

pl.subplot(2,2,3)
pl.plot(t1[0:iw1],'r-')

pl.subplot(2,2,4)
pl.plot(t2[0:iw2],'g-')


pl.show()

