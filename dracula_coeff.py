#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 06:50:35 2024

@author: ygaillard
"""
from sympy import Function, integrate, legendre
from sympy.abc import x

class N_lmn(Function):
    @classmethod
    def eval(cls, l, m, n):
        # if (l+m>=n) and (l>=m) and (l+m+n)%2==0:
            return (integrate( 
                legendre(l,x)*
                legendre(m,x)*
                legendre(n,x)
                ,(x,-1,1)))
        # else:
        #     return 0
    
def d_lmn(l, m, n):
    if (l+m<=n-2) or (l+m+n)%2==0:
        return 0
    else:
        return float(
            1/(2*n+1)*(
                (n+1)*N_lmn(l,m,n+1)+
                n*N_lmn(l-1,m,n)
                )
            )
    
def e_lmn(l, m, n):
    if (l+m<=n-2) or (l+m+n)%2==0:
        return 0
    else:
        return float(
            l*(l+1)/(2*l+1)*(
                N_lmn(l+1,m,n)
                -N_lmn(l-1,m,n)) 
            )
    
def a_lmn(l, m, n):
    if (l+m<=n-2) or (l+m+n)%2==0:
        return 0
    else:
        return 2*d_lmn(l,m,n)+e_lmn(l,m,n)

def k_lmn(l, m, n):
    if (l+m<=n-2) or (l+m+n)%2==0:
        return 0
    else:
        return float(
            k_lmn(l-2, m, n)
            -e_lmn(l, m, n)
            -(l-1)*(l-2)/(2*l-3)*N_lmn(l-3,m,n)
            +((l-1)*(l-2)/(2*l-3)-(2*l-1))*N_lmn(l-1,m,n)
                )

def s_lmn(l, m, n):
    if (l+m<=n-2) or (l+m+n)%2==0:
        return 0
    else:
        return float(
            m*(
                k_lmn(l, m, n)
               -k_lmn(l-1, m-1, n)
               -e_lmn(l-1, m-1, n)
               +l*N_lmn(l,m-1,n))
            )

def a_n(n):
    return -2*(n+1)*(n+2)

def b_n(n):
    return 4*(n+1)*(n+2) 

def c_n(n):
    return (n-1)*(n+2)*(n+1)*(n+4)

def A_lmn(l,m,n):
    return 2*d_lmn(n, m, l) ## !!! Care l and n switched not verified yet

def B_lmn(l,m,n):
    return -2*a_lmn(l, n, m) ## !!! Care n and m switched not verified yet

def C_lmn(l,m,n):
    return a_lmn(m, l, n)

def D_lmn(l,m,n):
    return -e_lmn(l, m, n)

def E_lmn(l,m,n):
    return -2*a_lmn(m, l, n)

def F_lmn(l,m,n):
    return (m*(m+1)-l*(l+1)+2)*a_lmn(m, l, n)+2*s_lmn(m, l, n)-4*k_lmn(n, l, m)-8

def G_lmn(l,m,n):
    return 4*(m*(m+1)+2)*a_lmn(l, m, n)-8*s_lmn(l, m, n)-16*k_lmn(l, m, n)


