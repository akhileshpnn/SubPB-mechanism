# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:57:58 2020

@author: nandan
"""

import numpy as np
from scipy.signal import gaussian

class single_gradient:
    
    t_beg1 = 0;
    
    beg_range=2; end_range= 2
    
    def add_stimulus(self, t):
        t_end1 = self.tF
        # t_end1=70
        extend1=((self.end_range-self.beg_range)*(t-self.t_beg1)**2/(t_end1-self.t_beg1)**2)+self.beg_range


        if self.t_beg1<=t<=t_end1:
            grad = self.stimu_strength*gaussian(self.N,extend1,sym=False)
        else:
            grad = np.zeros(self.N)
        return grad    

class single_gradient_transient:
    
    t_beg1 = 10;
    
    beg_range=2; end_range= 2
    
    def add_stimulus(self, t):
        t_end1=70
        extend1=((self.end_range-self.beg_range)*(t-self.t_beg1)**2/(t_end1-self.t_beg1)**2)+self.beg_range
        if self.t_beg1<=t<=t_end1:
            grad = self.stimu_strength*gaussian(self.N,extend1,sym=False)
        else:
            grad = np.zeros(self.N)
        return grad  


class gradient_reversal:
    
    t_beg1 = 0; t_end1 = 300;
    t_beg2 = 300
    
    beg_range=2; end_range= 2
    
    def add_stimulus(self, t):
        t_end2=self.tF
        extend1=((self.end_range-self.beg_range)*(t-self.t_beg1)**2/(self.t_end1-self.t_beg1)**2)+self.beg_range
        extend2=((self.end_range-self.beg_range)*(t-self.t_beg2)**2/(t_end2-self.t_beg2)**2)+self.beg_range

        
        if self.t_beg1<=t<=self.t_end1:  
            grad = self.threshold_stimu*gaussian(self.N,extend1,sym=False)     
        elif self.t_beg2<=t<=t_end2:  
            grad = self.stimu_strength*gaussian(self.N,extend2,sym=False)    
            grad=np.roll(grad,10)
        else:
            grad = np.zeros(self.N)
        return grad

    
class simultaneous_gradients:
    t_beg1=0;  
    beg_range1=2; end_range1=2 #static gradient
    beg_range2=2; end_range2=2 #static gradient
 
    
    def add_stimulus(self, t):
        t_end1=self.tF
        
        extend1=((self.end_range1-self.beg_range1)*(t-self.t_beg1)**2/(t_end1-self.t_beg1)**2)+self.beg_range1
        extend2=((self.end_range2-self.beg_range2)*(t-self.t_beg1)**2/(t_end1-self.t_beg1)**2)+self.beg_range2

        if self.t_beg1<=t<=t_end1:
            grad1 = self.threshold_stimu*gaussian(self.N,extend1,sym=False)
            grad2 = self.stimu_strength*gaussian(self.N,extend2,sym=False)
            # ligand2 = 4*egft*gaussian(self.N,1.75,sym=False)
            grad2=np.roll(grad2,10)  
            grad = grad1+grad2
            grad=np.round(grad,3)
            
            # grad=grad1
        else:
            grad = np.zeros(self.N)

        return grad

#######################################################################################
















class seq_gradient:
    t_beg1=10;t_end1=50
    t_beg2=90;t_end2=110;
    t_beg3=150;t_end3=190; 
 
    beg_range=2; end_range=5 #dynamic gradient
    beg_range2=2; end_range2=5 #dynamic gradient

    def add_stimulus(self, t):
        L = self.cellmem[-1]
        extend1=((self.end_range-self.beg_range)*(t-self.t_beg1)**2/(self.t_end1-self.t_beg1)**2)+self.beg_range
        extend2=((self.end_range2-self.beg_range2)*(t-self.t_beg2)**2/(self.t_end2-self.t_beg2)**2)+self.beg_range2
        extend3=((self.end_range-self.beg_range)*(t-self.t_beg3)**2/(self.t_end3-self.t_beg3)**2)+self.beg_range

        if self.t_beg1<=t<=self.t_end1:
            grad = self.stimu_strength*gaussian(self.N,extend1,sym=False)
            # grad = self.threshold_stimu*(1-self.cellmem/L)*np.ones(self.N)      
        elif self.t_beg2<=t<=self.t_end2:
            grad = self.stimu_strength*gaussian(self.N,extend2,sym=False)
        #     grad=np.roll(grad,10)
        elif self.t_beg3<=t<=self.t_end3:
            grad = self.stimu_strength*gaussian(self.N,extend3,sym=False)
            grad=np.roll(grad,10)
        else:
            grad = np.zeros(self.N)
    #    ligand=egft*np.ones(N)
    #    ligand=np.roll(ligand,5)
        return grad