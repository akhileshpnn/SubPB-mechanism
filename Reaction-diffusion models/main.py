# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 01:44:15 2020

@author: nandan
"""

import numpy as np
from scipy.integrate import solve_ivp
import sdeint
from scipy.linalg import block_diag
import matplotlib.pyplot as plt

from laplacebeltramioperator import *
from initialconditions import *
from stimulus_repo_rd import *
from model_repo_rd import *
from colormaps import parula_map


class ReactionDiffusion1D:
    
    R = 2; A = np.pi*R**2; L = 2*np.pi*R; # length of the cell membrane contour
    N = 20 # number of nodes
    tF = 251; # total time of intergration
    dt=0.01 # integration time step
    t_eval = np.arange(0,tF,dt)
   

    def __init__(self, model, initial_condition, lbo,stimulus):
        
        self.model = model # model with reaction terms and parameters
        self.lbo = lbo # specifies boundary condition
        self.initial_condition = initial_condition 
        self.stimulus=stimulus
    
    def initialize_system(self):
        
        
        self.stimulus.N=self.N
        self.cellmem = np.linspace(0,self.L,self.N) # define the cell membrane contour
        self.dsigma=self.cellmem[1]-self.cellmem[0] # spatial grid size
        
        self.du=self.model.du # setting diffusion constant of u as a class instance
        self.dv=self.model.dv # setting diffusion constant of v as a class instance
        if type(model).__name__=='Legi':
            self.dw=self.model.dw
            self.Z = np.zeros(3*self.N) # array for storing the simulated data for a single time point. 
            # 3*N because LEGI has three variables and each variable is simulated in N spatial bins
        else:
            self.Z = np.zeros(2*self.N) # wavepinning, Turing and SubPB have only two variables
        
        self.Z = self.initial_condition.set_initial_condition(self.model,self.N) # generates the intial conditions in the entire grid
        
        self.Stimulus=np.zeros((self.tF,self.N)) # stimulus values are initiated to be zero
                
    def get_input_profile(self):
        """
        generate stimulus profile
        
        self.stimulus class has the information about the type of stimulus. Please
        check 'stimulus_repo_rd.py' for more details.
        
        """
        for t in range(self.tF):
            self.Stimulus[t]=np.round(self.stimulus.add_stimulus(t),3)       
    

    def F_det(self,t,W):
        
        """
        for deterministic simulation of the system using solve_ivp when add_noise=None.
        
        The PDEs are converted to ODEs using method of lines and this fucntion does the
        implementation of this method. Please see the subsection 'Model implementation
        ' in the section 'Matrials and methods' of the main article for more information.
        
        Some changes in notation:
        dtheta is the step size. here it is dsigma.
        
        
        
        """
        
        stimulus_input=np.round(self.stimulus.add_stimulus(t),3) # import the stimulus values from the stimulus class
            
        LB = self.lbo.set_matrix(self.N) # the Laplace-Beltrami operator matrix for diffusion and specified boundary condition 
        LB = (1/self.dsigma**2)*LB # dsigma is the spatial step size and this captures the diffusion term in 
        
        if type(model).__name__=='Legi':
            A = block_diag(self.du*LB, self.dv*LB, self.dw*LB)  
            y = [W[:self.N],W[self.N:2*self.N],W[2*self.N:]]
            fu, fv,fw = self.model.reaction(t, y, stimulus_input)
            return np.matmul(A,W) + np.concatenate([fu, fv,fw]).transpose() # First part: diffusion. Second part: reaction
        else:
            A = block_diag(self.model.du*LB, self.model.dv*LB);
            y = [W[:self.N],W[self.N:]]
            fu, fv = self.model.reaction(t, y, stimulus_input)
        
            return np.matmul(A,W) + np.concatenate([fu, fv]).transpose()# First part: diffusion. Second part: reaction
    
    def F_stocha(self,W,t):
        
        """
        for stochastic integration of the system using sdeint when add_noise=True.
        This function is to flip the position of dependent variable(W) and independent variable (t).
        so that the reaction terms are compatible with Sde int package in python. 
        
        """
        return self.F_det(t,W)
        
    
    def G_stocha(self,W,t):
        
        """
        adding noise term to the system.
        returns block diagonal matrix with diaginal elements specifying the
        noise intensity
        
        """
        
        
        if type(model).__name__=='Legi':
            noise_ampl=0.0005
            Gu = Gw=np.diag(noise_ampl*np.ones(self.N)) # noise added to u and w
            Gv = np.zeros((self.N,self.N)) # no noise to v
            G = block_diag(Gu, Gv, Gw)
        else: 
            noise_ampl=0.005
            Gu = np.diag(noise_ampl*np.ones(self.N)) # noise is added only to u variable
            Gv = np.zeros((self.N,self.N)) # no noise added to v
            G = block_diag(Gu, Gv)
        return G
    
    def simulate(self):
        
        '''
        outputs:
            sol_det: deterministic solution. If add_noise=True, sol_det=None
            sol_stocha: stochastic solution. If add_noise=None, sol_stocha=None
        '''
        
        self.initialize_system()
        
        if self.add_noise==True:
            sol_stocha = sdeint.itoint(self.F_stocha, self.G_stocha, self.Z, tspan=self.t_eval)
            sol_det = None
        else:
            sol_det = solve_ivp(self.F_det, [0,self.tF], self.Z, t_eval=self.t_eval);
            sol_stocha=None
        
        self.get_input_profile()
        
        return sol_det,sol_stocha

    
    def plot_kymo(self,u,tt):
        
        '''
        plots the kymograph
        
        inputs:
            u: value of the activity array (2D)
            tt: array of timepoints that specifies the stimulus postions namely the onset and end.
        '''
        
        u=np.roll(u,shift=0,axis=0)
        fig,ax = plt.subplots()
        im = ax.imshow(u.T, extent=[0,self.N-1,len(u.T),0],cmap=parula_map, aspect = 'auto',vmin=np.min(u),vmax=np.max(u))
        ax.figure.colorbar(im)
        ax.set_ylabel('time(sec)',fontsize=20)
        ax.set_xlabel(r'plama membrane contour$(\theta)$',fontsize=20)
        ax.set_xticks([0,5,10,15,19])
        ax.set_xticklabels(['0',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2\pi$'],fontsize=20)
        ax.set_yticks(np.arange(0,len(u.T),50))
        ax.set_yticklabels(np.arange(0,len(u.T),50),fontsize=20)
        if tt is not None:
            for t in tt:
                ax.axhline(y=t,color='k',lw=5)        
        plt.show()
    
    def plot_profile(self,stimulus,tt):
            
        s_t= np.zeros((len(tt),len(stimulus))) 
        for i,t in enumerate(tt):
            s_t[i]=stimulus[:,t]

        plt.figure()
        plt.figure(figsize=(6,2))
        for i,t in enumerate(tt):        
            plt.plot(s_t[i],'g-',lw=2.0)
        plt.ylabel(r'$S(\theta)$',fontsize=20)
        # plt.ylim(0,0.02)
        plt.xlabel(r'plama membrane contour$(\theta)$',fontsize=20)
        ax=plt.gca()
        ax.set_xticks([0,5,10,15,19])
        ax.set_xticklabels(['0',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2\pi$'],fontsize=20)
        plt.show()    
    
    def plot_timeseries(self,u,bins,tt):
        
        [front,back]=bins
        response_front=u[front]
        response_back=u[back]
        time=np.arange(0,np.shape(u)[1])
        
        plt.figure()
        plt.figure(figsize=(9,3))
        plt.plot(time,response_front,'k-',lw=3.0)
        plt.plot(time,response_back,'k--',lw=3.0)
        if tt is not None:
            for t in tt:
                plt.axvline(x=t,color='k',lw=5.0)
        plt.ylabel(r'$u$',fontsize=20)
        plt.xlabel(r'$time(sec)$',fontsize=20)
        plt.xlim(0,time[-1])
        
        plt.show()       
    
        
if __name__ == '__main__':
     
    lbo = Periodic() # periodic boundary condition
    
    ## select the model to simulate and corresponding threshold stimulus strength (from Figure 3A)
    # model = WavePinning();threshold_stimu=0.003
    model = SubPB();threshold_stimu=0.012
    # model = Otsuji();threshold_stimu=0.001
    # model = Legi();threshold_stimu=0.005
       
    
    initial_condition = around_steadystate_2vars() # initial condition  
    # initial_condition = around_steadystate_legi()
    
    ## select type of stimulus
    # stimulus_type=single_gradient()
    stimulus_type=single_gradient_transient()
    # stimulus_type=gradient_reversal()
    # stimulus_type = simultaneous_gradients()
    
    # stimulus_type.stimu_strength=0.1 # set stimulus strength.
    stimulus_type.stimu_strength=2*threshold_stimu
    
    ## create class instance for the given model
    rd = ReactionDiffusion1D(model, initial_condition, lbo, stimulus_type)
    stimulus_type.tF=rd.tF
    
    add_noise=None
    rd.add_noise=add_noise # assigning the condition as a 'rd' class 
    
    ## integrate the system and find solution
    ## If add_noise=True, stochastic simulation. If add_noise=None, deterministic simulation.
    sol_det, sol_stocha=rd.simulate()
    
    if type(model).__name__=='Legi':
        
        if sol_det is None:
            out=sol_stocha[:,2*rd.N:].T             
        else:
            out=sol_det.y[2*rd.N:]
    else:
        if sol_det is None:
            out=sol_stocha[:,:rd.N].T             
        else:
            out=sol_det.y[:rd.N]
    
    tout=rd.t_eval[::int(1/rd.dt)]
    out = out[:,::int(1/rd.dt)]
    
    tt=[stimulus_type.t_beg1,stimulus_type.tF-1]
    rd.plot_profile(rd.Stimulus.T,tt)
    rd.plot_kymo(out,tt)     
    rd.plot_timeseries(out,[10,0],tt)
    
    
