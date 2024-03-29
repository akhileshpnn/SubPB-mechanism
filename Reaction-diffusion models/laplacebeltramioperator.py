# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 14:58:50 2020

@author: nandan
"""

import numpy as np
from scipy.sparse import spdiags

class LaplaceBeltrami:
    def matrix(self, N):
        return


class Periodic:
    
    def set_matrix(self, N):
        e = np.ones(N);
        LB = spdiags([e,-2*e,e],[-1,0,1],N,N).toarray(); 
        LB[0,N-1] = 1; 
        LB[N-1,0] = 1;        
        return LB

class Dirichlet:
    
    def set_matrix(self, N):
        e = np.ones(N);
        LB = spdiags([e,-2*e,e],[-1,0,1],N,N).toarray(); 
        return LB

class Neumann:
    
    def set_matrix(self, N):
        e = np.ones(N);
        LB = spdiags([e,-2*e,e],[-1,0,1],N,N).toarray(); 
        LB[0,0] = -1; 
        LB[N-1,N-1] = -1;
        return LB	
