# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 08:21:15 2023

@author: Akhilesh Nandan
"""
import numpy as np
from stimulus_repo import *

class wavepinning_1d:
    
    def reaction_terms(self,t, y,input_params):
            
        (total, sLmax, sRmax)=input_params
    
        k0 = 0.067 # s-1
        g = 1 # s-1 
        K = 1 # uM
        delta = 1 # s-1
        d = 0.01 # um2s-1
        d = 0.01 # um2s-1
        
        uL, uR = y 
        v=0.5*(2*total-(uL+uR))
        
        sL,sR=self.stimulus_type.add_stimulus(t,input_params)
    
        duL = v*(k0 + (g*uL*uL)/(K**2 + uL**2)) - delta*uL  -d*(uL-uR) + sL*v
        duR = v*(k0 + (g*uR*uR)/(K**2 + uR**2)) - delta*uR  -d*(uR-uL) + sR*v
        
        return np.array([duL, duR])
    
class legi_1d:
    
    def reaction_terms(self,t, y, input_params):
            
        (wt, sLmax, sRmax)=input_params
        
        k1a=2;k1b=2
        k2a=1;k2b=1;k3a=1;k3b=1
        du=0.5;dw=0.5
        
        wL,wR = y 
        
        sL,sR=self.stimulus_type.add_stimulus(t,input_params)
    
        vLs=(k2a/k2b)*(sL+sR)/2
        vRs=(k2a/k2b)*(sL+sR)/2
        
        uLs=(k1a/(2*du+k1b))*(sL+du*(sL+sR)/k1b)
        uRs=((k1a/k1b)*(sR+sL))-uLs
        
        dwL = k3a*uLs*(wt-wL)-k3b*vLs*wL-dw*(wL-wR)
        dwR = k3a*uRs*(wt-wR)-k3b*vRs*wR-dw*(wR-wL)
        
        return np.array([dwL, dwR])