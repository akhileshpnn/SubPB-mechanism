# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:50:09 2020

@author: nandan
"""

import numpy as np

class Model:
    def reaction(self, u, v, pher):
        return
    
    
    
class WavePinning:
    
    ctot = 2.268
    k0 = 0.067 # s-1
    g = 1 # s-1 # 
    K = 1 # uM
    d = 1 # s-1
    du = 0.1 # um2s-1
    dv = 10 # um2s-1

    def reaction(self, t, y, s):
        u, v = y 
        fu = v*(self.k0 + (self.g*u*u)/(self.K**2 + u**2)) - self.d*u + s*v
        fv = -fu
        return [fu, fv]


class SubPB:
    
    ctot=2.155 # ctot value used in the one-dimensional projection model is slightly different (ctot=2.21).
    # The difference in the Snpb is arising because in the one-dimensional projection model space is not
    # explicitly considered. One-dimensional model is only an approximation of the limiting case of full PDE model.
    
    k0 = 0.067 # s-1
    g = 1 # s-1 # 
    K = 1 # uM
    d = 1 # s-1
    du = 0.1 # um2s-1
    dv = 10 # um2s-1

    def reaction(self, t, y, s):
        u, v = y 
        fu = v*(self.k0 + (self.g*u*u)/(self.K**2 + u**2)) - self.d*u + s*v
        fv = -fu
        return [fu, fv]
    
class Otsuji:
    
    ctot = 2
    a1=2.5;a2=0.7;eps=0.01    
    du = 0.1 # um2s-1
    dv = 10 # um2s-1

    def reaction(self, t, y, s):
        u, v = y
        fu = self.a1*(v - (u+v)/(self.a2*(u+v)+1)**2) +s*v
        fv = -fu
        return [fu, fv]
    
class Legi:
    
    ctot = 1
    
    k1a=k1b=2
    k2a=k2b=1
    k3a=k3b=1

    du = 0.1 # um2s-1
    dv = 10 # um2s-1
    dw = 0.1 # um2s-1

    def reaction(self, t, y, s):
        u, v, w = y 
        fu = self.k1a*s-self.k1b*u
        fv = self.k2a*s-self.k2b*v
        fw = self.k3a*u*(self.ctot-w)-self.k3b*v*w
        return [fu, fv, fw]
    
