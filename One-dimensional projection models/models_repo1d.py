# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 08:21:15 2023

@author: Akhilesh Nandan
"""
import numpy as np


'''
see 'Description of the different cell polarity models' in the
'Materials and methods' in the paper for more details.

'''

class wavepinning_1d:
    
    def reaction_terms(self,t, y,input_params):
            
        (ctot, sLmax, sRmax)=input_params
    
        k0 = 0.067 # s-1
        g = 1 # s-1 
        K = 1 # uM
        delta = 1 # s-1
        d = 0.01 # um2s-1
        d = 0.01 # um2s-1
        
        uL, uR = y 
        v=0.5*(2*ctot-(uL+uR)) # mass conservation. 
        # assumes that due to large diffusivity, the value of compunent v is uniform
        # in both compartments.
        
        sL,sR=self.stimulus_type.add_stimulus(t,input_params)
    
        duL = v*(k0 + (g*uL*uL)/(K**2 + uL**2)) - delta*uL  -d*(uL-uR) + sL*v
        duR = v*(k0 + (g*uR*uR)/(K**2 + uR**2)) - delta*uR  -d*(uR-uL) + sR*v
        
        return np.array([duL, duR])
    
class legi_1d:
    
    def reaction_terms(self,t, y, input_params):
            
        (ctot, sLmax, sRmax)=input_params
        
        k1a=2;k1b=2
        k2a=1;k2b=1;k3a=1;k3b=1
        du=0.5;
        
        uL,uR = y 
        
        sL,sR=self.stimulus_type.add_stimulus(t,input_params)
    
        vLs=(k2a/k2b)*(sL+sR)/2
        vRs=(k2a/k2b)*(sL+sR)/2
        
        wLs=(k1a/(2*du+k1b))*(sL+du*(sL+sR)/k1b)
        wRs=((k1a/k1b)*(sR+sL))-wLs
        
        duL = k3a*wLs*(ctot-uL)-k3b*vLs*uL-du*(uL-uR)
        duR = k3a*wRs*(ctot-uR)-k3b*vRs*uR-du*(uR-uL)
        
        return np.array([duL, duR])
