import numpy as np
from scipy import constants
from fp_solver import fokker_planck, boundary

from models_repo import *

class QuasiPotentialLandscape:
    
    """
    estimation of quasi potential landscape of a given dynamical system (currently for
    optimized for two variable systems). codes adapted from....
    
    """
    
    ####### solver settings and parameters
    D=0.02
    drag = 1
    temperature=D*drag/constants.k
    
    ###### spatial grid parameters
    L=1.2 # length of the spatial domain
    h=0.02  # step size
    N=int(np.ceil(L/h)) # number of grid points
    
    def __init__(self, time_point,model,input_params):
        
        self.time_point=time_point
        self.model = model
        self.input_params=input_params
    
    def F(self,u1, u2):        
        y=[u1,u2]       
        return self.model.reaction_terms(self.time_point, y,self.input_params)

    def random_pdf(self):
    
        def pdf(*args):
            values = np.ones_like(args[0])
    
            for i, arg in enumerate(args):
                values *= np.random.uniform(0,1,(self.N,self.N))
            return values
    
        return pdf
    
    def find_potential(self):

        sim = fokker_planck(temperature=self.temperature, drag=self.drag, extent=[self.L, self.L],
                    resolution=self.h, boundary=boundary.reflecting, force=self.F)
    
        ### time-evolved solution
        pdf = self.random_pdf()
        p0 = pdf(*sim.grid)
    
        Tmax=100
        dt=0.01
        Nsteps = int(Tmax/dt)
        Pt = sim.propagate_interval(pdf, Tmax, Nsteps=Nsteps)
        
        return sim.grid,Pt
        
        
