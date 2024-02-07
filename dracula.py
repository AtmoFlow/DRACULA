#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 06:43:56 2024

@author: ygaillard
"""

import numpy as np
from scipy.integrate import solve_bvp
from dracula_coeff import a_lmn,A_lmn,B_lmn,C_lmn,D_lmn,E_lmn,F_lmn,G_lmn,a_n,b_n,c_n
from scipy.special import legendre
import JsonSim as JS

class draculaSolver:
    
    def __init__(self,order,eta,mu,Re,nPoints=100):
        self.order=order
        self.eta=eta
        self.mu=mu
        self.Re=Re
        
        self.Nt=2*self.order-1
        
        print('Nt:'+str(self.Nt))
        
        self.r = np.linspace(self.eta, 1.0, nPoints )
        
   
 
    def load(self):
        
        jsonSim=JS.JsonSim('dracula_output_order'+str(self.order)+'_eta'+str(self.eta)+'_mu'+str(self.mu)+'_Re'+str(self.Re)+'.json',
                   'order'+str(self.order)+'_eta'+str(self.eta)+'_mu'+str(self.mu)+'_Re'+str(self.Re),
                   access='r')
        if jsonSim.fileExisted:        
            self.f=jsonSim.siminformation['f']
            self.fp=jsonSim.siminformation['fp']
            self.g=jsonSim.siminformation['g']
            self.gp=jsonSim.siminformation['gp']
            self.gpp=jsonSim.siminformation['gpp']
            self.gppp=jsonSim.siminformation['gppp']
            
            del jsonSim
            return True
        else:
            print('No file available')
            return False
                
    
    def save(self):
        
        jsonSim=JS.JsonSim('dracula_output_order'+str(self.order)+'_eta'+str(self.eta)+'_mu'+str(self.mu)+'_Re'+str(self.Re)+'.json',
                   'order'+str(self.order)+'_eta'+str(self.eta)+'_mu'+str(self.mu)+'_Re'+str(self.Re),
                   access='rw')
        
        jsonSim.write('f', self.f)
        jsonSim.write('fp', self.fp)
        jsonSim.write('g', self.g)
        jsonSim.write('gp', self.gp)
        jsonSim.write('gpp', self.gpp)
        jsonSim.write('gppp', self.gppp)
        
        del jsonSim
            
        
    def solve(self):
       
        y = np.zeros( ( (int((self.Nt-1)/2)+1)*6, len( self.r ) ) )
        y[0] = np.linspace( self.eta**2/self.mu, 1, len( self.r ) )        # initial guess to fit boundary conditions

       
        result = solve_bvp( self.ode, self.bc, self.r, y, verbose=2 )
        rplot = np.linspace( self.eta, 1.0, 100 )

        self.f, self.fp, self.g, self.gp, self.gpp, self.gppp=np.zeros((6,self.Nt+1,np.shape(y)[1]))

        for n in range(0,int((self.Nt-1)/2)+1):
            self.f[2*n]=result.sol( rplot )[6*n+0]
            self.fp[2*n]=result.sol( rplot )[6*n+1]
            self.g[2*n+1]=result.sol( rplot )[6*n+2]
            self.gp[2*n+1]=result.sol( rplot )[6*n+3]
            self.gpp[2*n+1]=result.sol( rplot )[6*n+4]
            self.gppp[2*n+1]=result.sol( rplot )[6*n+5]
            
    def getFlowResult(self,theta):
        Psi=0
        Omega=0

        for n in range(0,self.Nt+1):
            if (n%2==0):  
                Omega+=np.outer(np.sin(theta)**2*np.polyval(legendre(n), np.cos(theta)),self.f[n])
            else:
                Psi+=np.outer(np.sin(theta)**2*np.polyval(legendre(n), np.cos(theta)),self.g[n])
                
        return (Psi,Omega)

        
        
    
    
    def ode(self, r, y ):  
            # y is [ f0, f0', g1, g1', g1'', g1''' ]
        f, fp, g, gp, gpp, gppp=np.zeros((6,self.Nt+1,np.shape(y)[1]))
        
        
        for n in range(0,int((self.Nt-1)/2)+1):
            f[2*n]=y[6*n]
            fp[2*n]=y[6*n+1]
            g[2*n+1]=y[6*n+2]
            gp[2*n+1]=y[6*n+3]
            gpp[2*n+1]=y[6*n+4]
            gppp[2*n+1]=y[6*n+5]
    
                
        EquList=[]
        
        for n in range(0,self.Nt+1):
            if (n%2==0):    
                subSum=0
                for l in range(0,self.Nt+1):
                    if l>=n+2:
                        subSum+=2*(2*n+1)/r**2*f[l]
                    for m in range(0,self.Nt+1):
                        if (m+l)>n-2:    
                            subSum+=(2*n+1)/2*self.Re/r**2*(a_lmn(l, m, n)*(g[l]*fp[m]-gp[m]*f[l]))
                RHSf = (n+1)*(n+2)/r**2*f[n]  + subSum
                
                EquList.extend((fp[n],RHSf))
            else:
                subSum=0
                for l in range(0,self.Nt+1):
                        if l>=n+2: 
                            subSum+=4*(2*n+1)/r**2*(
                                gpp[l]-2/r*gp[l]+(2-l*(l+1))/r**2*g[l]
                                )
                        for m in range(0,self.Nt+1):
                            if (m+l)>n-2:    
                                subSum+=(2*n+1)/2*self.Re/r**3*(
                                     A_lmn(l, m, n)*r*  fp[l]*f[m]
                                    +B_lmn(l, m, n)*    f[l]*f[m]
                                    +C_lmn(l, m, n)*r*  gppp[l]*g[m]
                                    +D_lmn(l, m, n)*r*  gpp[l]*gp[m]
                                    +E_lmn(l, m, n)*    gpp[l]*g[m]
                                    +F_lmn(l, m, n)/r*  gp[l]*g[m]
                                    +G_lmn(l, m, n)/r**2*g[l]*g[m]
                                    )
                              
                RHSg = -a_n(n)/r**2*gpp[n] -b_n(n)/r**3*gp[n] -c_n(n)/r**4*g[n]+subSum
                
                EquList.extend((gp[n],gpp[n],gppp[n],RHSg))
               
        return np.vstack( EquList )
    
    
    def bc(self,ya, yb ):
        bcList=[]
        for n in range(0,int((self.Nt-1)/2)+1):
            if n==0:
                bcList.extend((ya[0]-self.eta**2/self.mu, ya[2], ya[3],      yb[0] - 1 , yb[2], yb[3]))
            else:
                bcList.extend((ya[6*n+0], ya[6*n+2], ya[6*n+3],      yb[6*n+0]  , yb[6*n+2], yb[6*n+3]))   
        
        return np.array( bcList )
    

        