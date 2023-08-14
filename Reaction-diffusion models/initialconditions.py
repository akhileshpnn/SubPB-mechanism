# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:08:26 2020

@author: nandan
"""

import numpy as np
from scipy.integrate import solve_ivp

class InitialCondition:

    def set_initial_condition(self, N):
        return
    

class around_steadystate_2vars:
    
    zeta_per=0.01 # amplitude of random perturbation around the homogenous steady state
    
    def set_initial_condition(self, model,N):
        '''
        Parameters
        ----------
        model : name of the reaction diffusion model class
        N : Number of bins used in the simulation

        Returns
        -------
        W0 : numpy array of initial values. Shape: (1 , 2 * N)

        '''
        
        f=model
        uin=f.ctot*0.5;vin=f.ctot*0.5
        s=0
        sol_timeseries = solve_ivp(f.reaction, [0, 500], [uin, vin], args=(s,), dense_output=True)        
        [uss,vss] = sol_timeseries.y[:,-1] # are the homogeneous steady state of the system
        
        # generate initial condition for the RD system
        e = np.ones(N);
        We = (np.concatenate([uss*e, vss*e])).transpose()
        random_num=np.random.uniform(0,1,N)
        Wper = (np.concatenate([self.zeta_per*(1-random_num), self.zeta_per*(random_num-1)])).transpose()
        W0 = We + Wper;
        if len(np.argwhere(W0<0))!=0:
            '''
            sometimes the Wper results in negative number which doesnot reflect physiological conditions.
            decrease the zeta_per parameter in that case
            '''
            W0=np.nan
            print('Negative numbers in initial values. Decrease zeta_per')
        else:
            return W0

class around_steadystate_legi:
    
    '''
    please follow notations are explanations as in the class 'around_steadystate_2vars'.
    Everything is the same except that this class if for generating initial conditions
    for Legi model and the output W0 has a shape (1 , 3 * N)
    '''
    
    zeta_per=0.005
    
    def set_initial_condition(self, model, N):
        f=model
        uin=0.0;vin=0.0;win=0.45
        s=0
        sol_timeseries = solve_ivp(f.reaction, [0, 5000], [uin, vin, win], args=(s,), dense_output=True)        
        [uss,vss,wss] = sol_timeseries.y[:,-1]

        
        e = np.ones(N);
        We = (np.concatenate([uss*e, vss*e, wss*e])).transpose()
        Wper = (np.concatenate([self.zeta_per*(1-np.random.rand(N)), self.zeta_per*(1-np.random.rand(N)), self.zeta_per*(1-np.random.rand(N))])).transpose()
        W0 = We + Wper;
        if len(np.argwhere(W0<0))!=0:
            W0=np.nan
            print('Negative numbers in initial values. Decrease zeta_per')
        else:
            return W0        
    
