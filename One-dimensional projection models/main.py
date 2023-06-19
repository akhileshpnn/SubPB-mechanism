# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:42:45 2023

@author: Akhilesh Nandan
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
from math import isclose
import os
import matplotlib.pylab as pylab
import warnings
warnings.filterwarnings("ignore")

params = {'legend.fontsize': 15,
          'axes.labelsize': 20,
          'axes.labelpad' : 15,
          'axes.titlesize':20,
          'xtick.labelsize':20,
          'ytick.labelsize':20}
pylab.rcParams.update(params)

from models_repo import *
from stimulus_repo import *
from utils import *
from quasi_potential_landscape import *

class Model1D:
    
    tF = 251;
    t_eval = np.arange(0,tF,0.01)

    def __init__(self, model,stimulus_type):
        
        self.model = model
        self.stimulus_type=stimulus_type
    
    # def initialize_system(self):
        
    
    def add_stimulus(self):
        self.stimulus_left=np.zeros(len(self.t_eval))
        self.stimulus_right=np.zeros(len(self.t_eval))
        for i, t in enumerate(self.t_eval):
            self.stimulus_left[i]=self.stimulus_type.add_stimulus(t,self.input_params)[0]
            self.stimulus_right[i]=self.stimulus_type.add_stimulus(t,self.input_params)[1]
        
        
    def solve_timeseries(self):
        
        self.add_stimulus()
        uLin,uRin=initial_condition
        
        sol_timeseries = solve_ivp(self.model.reaction_terms, [0, self.tF], [uLin, uRin], t_eval=self.t_eval,args=(self.input_params,), dense_output=True)
              
        [uLss,uRss] = sol_timeseries.y[:,-1]
        
        uL=sol_timeseries.y[0]
        uR=sol_timeseries.y[1]
        
        return uL,uR
    
    
    
    def plot_timeseries(self,x,y1,y2):
        
        self.vertical_lines = np.linspace(self.stimulus_type.stimulus_beg,self.stimulus_type.stimulus_end,6)

        plt.figure()
        for i in range(len(self.vertical_lines)):
            plt.axvspan(self.vertical_lines[i], self.vertical_lines[-1-i], alpha=0.2, color='green')
        plt.xlabel('time(sec)')
        plt.plot(x,y1,'k-',lw=2.0,label='$u_{left}$')
        plt.plot(x,y2,'k--',lw=2.0,label='$u_{right}$')
        plt.xlim(0,self.tF)
        plt.ylabel('response')
        plt.ylim(0,1)
        plt.legend()
        plt.show()
        
    def plot_statespace(self,time_point,traj=None,Qlimits=None):
        
        sa=StabilityAnalysis(time_point,self.model,self.input_params)
        
        stable_fps,unstable_fps=sa.estimate_fixed_points('.1f')
     
        ### plotting state space with separatrix and trajectories

        X_ss=sa.grid_ss[0];Y_ss=sa.grid_ss[1]

        func = lambda u1,u2 :self.model.reaction_terms(time_point,[u1,u2],self.input_params)

        u,v=func(X_ss,Y_ss)
        
        fig, ax = plt.subplots()
        ax.streamplot(X_ss,Y_ss,u,v,density=0.5,color=[0.1,0.1,0.1,0.1])
        if not traj is None:
            uL,uR=traj
            if time_point>40:
                ax.plot(uL[time_point-40:time_point],uR[time_point-40:time_point],lw=3.0,c='blue')
                ax.plot(uL[time_point-40:time_point][-1],uR[time_point-40:time_point][-1],marker='^',ms=10,c='magenta')     
            else:
                ax.plot(uL[0:time_point],uR[0:time_point],lw=3.0,c='blue')
                ax.plot(uL[0:time_point][-1],uR[0:time_point][-1],marker='^',ms=10,c='magenta')
        
        ## plot stable and unstable fixed points
        for usfp in unstable_fps:
            sep=sa.finding_separatrix(usfp,eps=1e-2,uniform_pert=True)
            sep_uL=np.concatenate((np.flip(sep[0][0]),sep[1][0]))
            sep_uR=np.concatenate((np.flip(sep[0][1]),sep[1][1]))
            # ax.plot(sep_uL,sep_uR,ls='--',color='gray',lw=3)

        for kk in range(len(stable_fps)):
            ls,rs=stable_fps[kk]
            if isclose(ls,rs):        
                ax.scatter(ls,rs,marker='o',s=100,color='black',edgecolors='black',alpha=1)
            else:
                ax.scatter(ls,rs,marker='s',s=100,color='black',edgecolors='black',alpha=1)
        if len(unstable_fps)!=0:
            ax.scatter(unstable_fps[:,0],unstable_fps[:,1],marker='o',s=100,color='gray',edgecolors='black',alpha=1)
        if Qlimits is not None:
            Qmin,Qbound=Qlimits    
            cset = ax.contour(self.grid_pot[0],self.grid_pot[1], self.Q,levels=np.linspace(Qmin,Qbound,18), zdir='z', offset=0, cmap='coolwarm',alpha=0.25)
    

        ax.set_xlim(0,1.2)
        ax.set_ylim(0,1.2)

        ax.set_xticks(ticks=[0,0.5,1])
        ax.set_yticks(ticks=[0,0.5,1])
        ax.set_aspect('equal')  
        ax.set_xlabel('$u_{left}$')
        ax.set_ylabel('$u_{right}$')
        # ax.set_title('sL='+str(stimulus_left[time_point])+', time='+str(time_point))
        ax.set_title('t= '+str(time_point)+'(sec)',pad=10)
        plt.show()
    
    def plot_3D(self,grid_pot,Q,zlims,Qlimits): 
        
        Xpot,Ypot=grid_pot
        Qmin,Qbound=Qlimits   
        
        Q_contour=Q.copy()
        Q[Q>zlims[1]]=np.nan
        Q=Q-(np.nanmin(Q)-0.5)
    
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.rcParams.update({'font.size': 10})
        surf = ax.plot_surface(Xpot,Ypot, Q,  rstride=1, cstride=1, cmap='hot',alpha=0.3,
                    linewidth=0.25,antialiased=True,edgecolor='k',vmin=np.nanmin(Q),vmax=np.nanmin(Q)+2)
        cset = ax.contour(Xpot,Ypot, Q_contour,levels=np.linspace(Qmin,Qbound,18), zdir='z', offset=0, cmap='coolwarm',alpha=0.25)

        # uL,uR=traj
        # if time_point>40:
        #     ax.plot3D(uL[time_point-40:time_point],uR[time_point-40:time_point],np.zeros(len(uR[time_point-40:time_point])),lw=3.0,c='blue')
        #     ax.plot3D(uL[time_point-40:time_point][-1],uR[time_point-40:time_point][-1],0,marker='o',ms=10,c='blue')
        # else:
        #     ax.plot3D(uL[0:time_point],uR[0:time_point],np.zeros(len(uR[0:time_point])),lw=3.0,c='blue')
        #     ax.plot3D(uL[0:time_point][-1],uR[0:time_point][-1],0,marker='o',ms=10,c='blue')
           
    
        # ax.scatter(X_bndry,Y_bndry,Us_bndry,'o',color='k')
        ax.set_zlim([0,np.nanmin(Q)+1])
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        # ax.set(xlabel='$u_{left}$', ylabel='$u_{right}$')
        ax.set_xticks(ticks=[0,0.5,1])
        ax.set_yticks(ticks=[0,0.5,1])
        ax.set_zticks(ticks=[0,1])
        ax.set_yticklabels(labels=[])
        ax.set_xticklabels(labels=[])
        ax.set_zticklabels(labels=[])
        ax.view_init(elev=26, azim=-85)
        ax.set_aspect('auto')
        plt.show()             
        
        
if __name__ == '__main__':
    
    ## select model equation and type of stimulus from model_repo.py and stimulus_repo.py respectively
    
    
    """
    Wave-pinning model parameters
    
    """
    model = wavepinning_1d()
    stimulus_type=stimulus_ramp()
    # stimulus_type=stimulus_single()
    model.stimulus_type=stimulus_type
    
    ## simulation conditions
    initial_condition = [0.1,0.05]
    
    # model parameters 
    sLmax=0.02 # maximum stimulus strength on the left side of the cell (sleft)
    sRmax=0.0 # maximum stimulus strength on the right side of the cell (sright)

    # total=2.15 # region I
    total=2.21 # region II, criticality  
    # total=2.26 # region III
    # total=2.32 # region IV

    input_params=[total,sLmax,sRmax]
    
    """
    Legi model parameters
    
    """
    
    # model = legi_1d()
    # stimulus_type=stimulus_ramp()
    # # stimulus_type=stimulus_single()
    # model.stimulus_type=stimulus_type
    
    
    # sL=10
    # sR=0.1

    # # system parameters
    
    # ## simulation conditions
    # initial_condition = [0.2,0.1]
    
    # # model parameters 
    # sLmax=10 # maximum stimulus strength on the left side of the cell (sleft)
    # sRmax=0.1 # maximum stimulus strength on the right side of the cell (sright)

    # total=1   

    # input_params=[total,sLmax,sRmax]
    
    
    ####################################### run the model
    m1d=Model1D(model,stimulus_type) # generating one-dimensional projection model instance
    m1d.input_params=input_params # contains bifurcation parameter and stimulus amplitudes
             
    uL,uR=m1d.solve_timeseries()
    
    ## cropping the time series
    uL=uL[::100]
    uR=uR[::100]
    t_eval=m1d.t_eval[::100]
    
    m1d.stimulus_left=m1d.stimulus_left[::100]
    m1d.stimulus_right=m1d.stimulus_right[::100]
       
    m1d.plot_timeseries(t_eval,uL,uR) # plot solution time series
    
    # time points to plot state space snapshots as in Figure1G
    time_points_plot=np.pad(m1d.vertical_lines, (1, 2), 'constant', constant_values=(0,0))
    time_points_plot[0]=9
    time_points_plot[-3]=100
    time_points_plot[-2]=125
    time_points_plot[-1]=250

    plot_idx=7 # change this index to change the snapshot time
    time_point=int(time_points_plot[plot_idx])########################################################################

    
    """
    estimating quasi-potential from steady state probability distribution. For more details refer,
    
    Wang, J., Xu, L., Wang, E., and Huang, S. (2010). The potential landscape of genetic circuits imposes the arrow of time 
    in stem cell differentiation. Biophysical journal, 99:29–39.
    
    """

    try:
        folder_load=os.path.abspath(os.getcwd())+'\\qpl output\\'
        Pt=np.load(os.path.join(folder_load,'Probability_'+str(total)+'_'+str(m1d.stimulus_left[time_point])+'.npy'))
        grid_pot=np.load(os.path.join(folder_load,'Grid_'+str(total)+'_'+str(m1d.stimulus_left[time_point])+'.npy'))      
    except:
        print('estimating quasi-potential landscape. This might take sometime.')
        qpl=QuasiPotentialLandscape(time_point,model,input_params)
        grid_pot,Pt=qpl.find_potential()  
         
        # folder_save=os.path.abspath(os.getcwd())+'\\qpl output\\'    
        # np.save(os.path.join(folder_save,'Probability_'+str(self.input_params[0])+'_'+str(s)+'.npy'),Q)
        # np.save(os.path.join(folder_save,'Grid_'+str(self.input_params[1])+'_'+str(s)+'.npy'),sim.grid)
    
    Q= -np.log(Pt) # quasi-potential value
    zlims=[np.min(Q)-0.5,np.min(Q)+1]
    
    bb=BasinBoundary(Q,grid_pot) #estimating asymptotically attracting state space region from Q.
    bb.model=model
    Qmin,Qbound=bb.find_boundary()
    m1d.grid_pot=grid_pot
    m1d.Q=Q
    
    """
    plotting state space trajectory with landscape projections
    
    """
    
    traj=[uL,uR]
    # traj=None
    Qlimts=[Qmin,Qbound]
    
    m1d.plot_statespace(time_point,traj,Qlimts)
    m1d.plot_3D(grid_pot,Q,zlims,Qlimts)
    
    


    

           