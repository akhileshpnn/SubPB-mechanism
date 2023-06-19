# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:08:26 2020

@author: nandan
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class InitialCondition:

    def set_initial_condition(self, N):
        return
    

class around_steadystate_2vars:
    
    dsigma=0.01
    def set_initial_condition(self, model,N):
        f=model
        uin=f.total*0.5;vin=f.total*0.5
        s=0
        sol_timeseries = solve_ivp(f.reaction, [0, 500], [uin, vin], args=(s,), dense_output=True)
#        
        [uss,vss] = sol_timeseries.y[:,-1]
        
        # u=sol_timeseries.y[0]
        # v=sol_timeseries.y[1]
        
        # plt.figure()
        # plt.title('Finding uniforn SS')
        # plt.plot(sol_timeseries.t,u,'r-',lw=2.0,label='u')
        # plt.plot(sol_timeseries.t,v,'k-',lw=2.0,label='v')
        # plt.legend()
        # plt.show()
        
        e = np.ones(N);
        # ess,pss=0.2,0.2
        We = (np.concatenate([uss*e, vss*e])).transpose()
        #We = np.array([ce*e; ge*e]).transpose()
        # random_num=np.random.rand(N)
        random_num=np.random.uniform(0,1,N)
        Wper = (np.concatenate([self.dsigma*(1-random_num), self.dsigma*(random_num-1)])).transpose()
        W0 = We + Wper;
        # W0 = We - Wper
        if len(np.argwhere(W0<0))!=0:
            W0=np.nan
            print('Negative numbers in initial values')
        else:
            return W0

class around_steadystate_legi:
    
    dsigma=0.005
    def set_initial_condition(self, model, N):
        f=model
        uin=0.0;vin=0.0;win=0.45
        s=0
        sol_timeseries = solve_ivp(f.reaction, [0, 5000], [uin, vin, win], args=(s,), dense_output=True)
#        
        [uss,vss,wss] = sol_timeseries.y[:,-1]
        # 
        # wss=0.45
        
        u=sol_timeseries.y[0]
        v=sol_timeseries.y[1]
        w=sol_timeseries.y[2]
        
        # plt.figure()
        # plt.title('Finding uniforn SS')
        # plt.plot(sol_timeseries.t,u,'r-',lw=2.0,label='u')
        # plt.plot(sol_timeseries.t,v,'k-',lw=2.0,label='v')
        # plt.plot(sol_timeseries.t,w,'b-',lw=2.0,label='w')
        # plt.legend()
        # plt.show()
        
        e = np.ones(N);
        # ess,pss=0.2,0.2
        We = (np.concatenate([uss*e, vss*e, wss*e])).transpose()
        #We = np.array([ce*e; ge*e]).transpose()
#        random_num=np.random.rand(N)
        Wper = (np.concatenate([self.dsigma*(1-np.random.rand(N)), self.dsigma*(1-np.random.rand(N)), self.dsigma*(1-np.random.rand(N))])).transpose()
        W0 = We + Wper;
        if len(np.argwhere(W0<0))!=0:
            W0=np.nan
            print('Negative numbers in initial values')
        else:
            return W0        
    
