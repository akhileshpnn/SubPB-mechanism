# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 12:20:21 2023

@author: Akhilesh Nandan
"""
import numpy as np
from scipy.optimize import fsolve
import sympy as sp
import itertools
from scipy.integrate import solve_ivp
from scipy import ndimage
from skimage.segmentation import find_boundaries


class StabilityAnalysis:
    
    '''
    given a 2 variable ode model, this class numerically finding the fixed points, evaluates the 
    jacobian around those points, infer their stability.
    '''
    
    
    # defining grid in the state space for estimating possible steady states
    uL_range=np.linspace(0,1.2,101)
    uR_range=np.linspace(0,1.2,101)
    grid_ss = np.meshgrid(uL_range, uR_range)
  
    def __init__(self, timepoint,model,input_params):
        
        self.model = model
        self.timepoint=timepoint
        self.input_params=input_params
        
    def find_root(self,root_guess):
        
        '''
        given the fixed point guesses, this function finds the state space points where the
        given ode has a root using 'fsolve' function 
        '''
          
        func = lambda y : self.model.reaction_terms(self.timepoint,y,self.input_params) # t=10 because of sL          
        # roots=scipy.optimize.fixed_point(func, root_guess,xtol=1e-8,maxiter=100,method='iteration')
        roots = fsolve(func, x0=root_guess,xtol=1e-10, maxfev=0,factor=10,epsfcn=1e-15)
        roots = np.round(roots,2)
        return roots
    
    def partial_derivatives(self):
        
        '''
        defines the symbolic partial derivatives in the jacobian. Here for 2 variable system 4 of the 
        jacobian terms are calculated.
        
        u1, u2 are the general variables
        G1,G2 are the symbolic reaction terems.
        G1_u1 is the partial derivative of G1 qith respect to u1.
        similarly for others

        '''
        
        u1, u2= sp.symbols('u1 u2')
    
        func = lambda u1,u2 : self.model.reaction_terms(self.timepoint,[u1,u2],self.input_params)
        
        G1 = lambda u1, u2 : func(u1,u2)[0]
        G2 = lambda u1, u2 : func(u1,u2)[1]
            
        G1_u1 = sp.lambdify((u1, u2),sp.diff(G1(u1, u2),u1))
        G1_u2 = sp.lambdify((u1, u2),sp.diff(G1(u1, u2),u2))
        
        G2_u1 = sp.lambdify((u1, u2),sp.diff(G2(u1, u2),u1))
        G2_u2 = sp.lambdify((u1, u2),sp.diff(G2(u1, u2),u2))
        
        return [G1_u1,G1_u2,G2_u1,G2_u2]
        
    def check_stability(self,uLs,uRs,round_float):
        
        [G1_u1,G1_u2,G2_u1,G2_u2]=self.partial_derivatives()
        
        a11 = G1_u1(uLs,uRs)
        a12 = G1_u2(uLs,uRs)
        a21 = G2_u1(uLs,uRs)
        a22 = G2_u2(uLs,uRs)
        J=[[a11,a12],[a21,a22]]
    
        eig_value = np.linalg.eig(J)[0]
        # max_eigevalue=np.max(eig_value).round(2)
        max_eigevalue=float(format(np.max(eig_value), round_float))  ## 0.1f was not identifying the unstable steady state for IHSS regime
        if max_eigevalue==0:
            ss_type = 'unknown'
        else:
            ss_eig_sign = np.sign(max_eigevalue) 
            
            if ss_eig_sign>=0:
                ss_type = 'unstable'
            elif ss_eig_sign<0:
                ss_type='stable'
    
        return ss_type
    
    
    def estimate_fixed_points(self,round_float):
        
        '''
        given the ode system, this fuction numerically search for possible fixed pooints 
        in the given grid and infer its stability.
        
        input:
            round_float: the precision for rounding off the eigenvalue
        return:
            separate list of stable and usnatbl;e fixed points
        
        '''
        
        uLgs=self.uL_range[::50]
        uRgs=self.uR_range[::50]
        
        root_guesses=np.zeros((len(uLgs),2))
        root_guesses[:,0]=uLgs
        root_guesses[:,1]=uRgs
        
        
        fps_dummy=[]
        
        for root_guess in itertools.product(uLgs, uRgs):
            # root_guess=root_guesses[i]
            [uLs,uRs]=self.find_root(root_guess=root_guess)      
            fps_dummy.append([uLs,uRs])
            
        fps_dummy=np.array(fps_dummy).reshape((len(fps_dummy),2))
        fps_dummy=np.unique(fps_dummy,axis=0)
        
        fps=[]
        stab=[]
        
        for i in range(len(fps_dummy)): 
            stability = self.check_stability(fps_dummy[i][0],fps_dummy[i][1],round_float)
            if stability!='unknown':
                stab.append(stability)
                fps.append(fps_dummy[i])
        
        fps=np.array(fps)
        
        unstable_fps=[]
        stable_fps=[]
        
        for j in range(len(stab)):
            if stab[j]=='unstable':  
                unstable_fps.append(fps[j])
            else:
                stable_fps.append(fps[j])
     
        return np.array(stable_fps),np.array(unstable_fps)
        
    def finding_separatrix(self,usfp,eps,uniform_pert=False):
        # finding seperatrix
    
        def new_reaction_term(t, y):
            uL,uR = y
            if self.uL_range[0] < uL < self.uL_range[-1] and self.uR_range[0] < uR < self.uR_range[-1]:
                return -self.model.reaction_terms(t, y,self.input_params)
            else:
                return np.array([0, 0])
        
        separatrices=[]
        
        tF_temp=100
        t_eval_temp=np.arange(0,tF_temp,0.001)
        
        uLin,uRin=usfp+[eps,0]
        if uniform_pert==True:
            uLin,uRin=usfp+[eps,eps]   
        sol_timeseries = solve_ivp(new_reaction_term, [0, tF_temp], [uLin, uRin], t_eval=t_eval_temp, dense_output=True)
        uLs1=sol_timeseries.y[0][::1000]
        uRs1=sol_timeseries.y[1][::1000]
        
        uLin,uRin=usfp-[5*eps,0]
        if uniform_pert==True:
            uLin,uRin=usfp-[eps,eps]
        sol_timeseries = solve_ivp(new_reaction_term, [0, tF_temp], [uLin, uRin], t_eval=t_eval_temp, dense_output=True)
        uLs2=sol_timeseries.y[0][::1000]
        uRs2=sol_timeseries.y[1][::1000]
        
        separatrices.append([uLs1,uRs1])
        separatrices.append([uLs2,uRs2])
        
        return separatrices
    
class BasinBoundary:
    
    """
    This class is used to find the quasi-potential value (Qbound) at the boundary 
    of basin of attraction of steady state using gaussin curvature method. Given the 
    quasi-potential values for the parameter setting, the objects in this class estimates
    gaussin curvature and finds the boundary.
    
    """
    
    def __init__(self,Q,grid_pot):
        
        self.Q = Q
        self.grid_pot=grid_pot
        
    def gaussian_curvature(self):
        
        self.Xpot=self.grid_pot[0];self.Ypot=self.grid_pot[1]

        Z=self.Q.copy()
        xmin=self.Xpot[0,0];xmax=self.Xpot[-1,0]
        ymin=self.Ypot[0,0];ymax=self.Ypot[0,-1]
        self.h=self.Xpot[1,0]-self.Xpot[0,0]
        
        Zy, Zx = np.gradient(Z,self.h)   # axis zero is y direction, axis=1 is x direction                                                  
        Zxy, Zxx = np.gradient(Zx,self.h)                                                  
        Zyy, _ = np.gradient(Zy,self.h)                                                    
        K = (Zxx * Zyy - (Zxy ** 2)) /  (1 + (Zx ** 2) + (Zy **2)) ** 2             
        return K
    
    def find_boundary(self):
        
        Z=self.Q.copy()
        K=self.gaussian_curvature()
    
        ##################  identifying blobs from curvature
        
        """
        threshold is a free parameter in the method that determines the Qbound.
        Thereshold curvature value is used to find high curvature areas that maybe
        a steady state. For a new model one needs to optimize this thereshold value by trial and error.
        
        Here threshold is set to np.mean(K)+0.1*np.std(K) (for wavepinning model) and 
        np.nanmean(K)/1000 (for legi model).
        
        """
        
        if type(self.model).__name__=='wavepinning_1d':     
            blobs = K > np.mean(K)+0.1*np.std(K) 
        elif type(self.model).__name__=='legi_1d':  
            blobs = K>np.nanmean(K)/1000
        else:
            blobs = K > np.mean(K)+0.1*np.std(K)
        
        # label connected regions that satisfy this condition
        labels, nlabels = ndimage.label(blobs)
        xc, yc = np.vstack(ndimage.center_of_mass(K, labels, np.arange(nlabels) + 1)).T
    
        # fig, ax = plt.subplots()
        # plt.rcParams.update({'font.size': 20})
        # im=ax.imshow(K,cmap='hot',origin='lower',extent=[xmin,xmax,ymin,ymax],vmin=-10,vmax=10)
        # plt.colorbar(im)
        # ax.set_xlim([0,1])
        # ax.set_ylim([0,1])
        # ax.set(xlabel='$u_{left}$', ylabel='$u_{right}$')
        # ax.set_xticks(ticks=[0,0.5,1])
        # ax.set_yticks(ticks=[0,0.5,1])
        # ax.set_aspect('auto')
        # plt.show()
    
    
        ### first and second order partial derivatives 
        Zy,Zx=np.gradient(Z,self.h) # axis 0 is y; axis 1 is x
        Zyy,Zyx=np.gradient(Zy,self.h)
        Zxy,Zxx=np.gradient(Zx,self.h)
    
        Contour_max=[]
        
        '''
        go through each of the labelled regions and apply the condition of
        slopes. For stable fixed point slopes are distributed aroud zero.
        see the main article for more details.
        '''
        
        for idx in range(len(xc)):
               
            Zy_b=Zy[np.where(labels==idx+1)] # the slope values in the y direction in the lablled region
            Zx_b=Zx[np.where(labels==idx+1)] # the slope values in the x direction in the lablled region
            
            ## visualize these slope distributions if necessary
            # plt.figure()
            # plt.hist(Zx_b,bins=np.arange(-2,2,0.1))
            # plt.title('Slope along uleft')
            # plt.show()
            
            # plt.figure()
            # plt.hist(Zy_b,bins=np.arange(-2,2,0.1))
            # plt.title('Slope along uright')
            # plt.show()
            
            '''
            checking the slope distributions.
            '''
            if (np.min(Zx_b)<0 and  np.min(Zy_b)<0) and (np.max(Zx_b)>0 and  np.max(Zy_b)>0):
                
                '''
                this if loop ensures that only the labelled regions with slope histogram 
                that lie across zero value is considered. implicitly this is a stable well in the potential landscape.
                '''
                
                labelled_image_new=np.zeros(np.shape(labels),dtype='int32')
                labelled_image_new[np.where(labels==idx+1)]=1
                
                bndry=find_boundaries(labelled_image_new, mode='outer', background=0)*1
                Q_bndry=self.Q[np.where(bndry==1)]
                
                Qm=np.mean(Q_bndry)
                Contour_max.append(Qm)
    
        Qmin=np.min(self.Q)
        Qbound=np.max(Contour_max)
        
        return Qmin,Qbound