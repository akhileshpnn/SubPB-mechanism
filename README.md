# SubPB-mechanism


## <div align="center"><b><a href="README.md">SubPB</a></b></div>

<div align="center">

🚩 [**Paper**](#-Paper) **|** ⚡ **Quick Guide** **|** 


</div>

---


## 🚩 Paper

This repository is for the models described in the paper:

"Non-asymptotic transients away from steady states determine cellular responsiveness to dynamic spatial-temporal signals
" [[biorRXiv]]([https://www.biorxiv.org/content/10.1101/2023.06.01.543361v1](https://www.biorxiv.org/content/10.1101/2023.02.03.526969v1)) 

## 🔧 Dependencies 
  - Python >= 3.6 
  - Numpy == 1.22.3
  - Scipy == 1.7.3

## ⚡ Quick Guide

Contains two folders.

1. One-dimensional projection models: 
Recalcuates evolution of state-space tracjectories with quasi-potential landscapes and their projections in Figs 1 and 2 (including corresponding supplementary figures) of the paper.

```python
# file 'main.py' has class Model1D that has all attributes for integrating the system, estimating quasi-potential landscape and plotting them.
# example with one-dimensional projection of Wave-pinning model. Variable notations same as in text. Any new model can be added to 'models_repo1d.py'. Currently works only for models with two variables.

 model = wavepinning_1d()
stimulus_type=stimulus_ramp()
model.stimulus_type=stimulus_type

## simulation conditions
initial_condition = [0.1,0.05]

# model parameters 
sLmax=0.02 # maximum stimulus strength on the left side of the cell (sleft)
sRmax=0.0 # maximum stimulus strength on the right side of the cell (sright)

# values are taken from Fig 1D
ctot=2.21 # region II, criticality  

# integrating the ODEs

m1d=Model1D(model,stimulus_type) # generating one-dimensional projection model instance
m1d.input_params=input_params # contains bifurcation parameter and stimulus amplitudes
uL,uR=m1d.solve_timeseries()

# Estimating quasi-potential
qpl=QuasiPotentialLandscape(time_point,model,input_params)
grid_pot,Pt=qpl.find_potential() 
Q=-np.log(Pt)
```
![one_dimensional_projection](https://github.com/akhileshpnn/SubPB-mechanism/assets/41164857/2c3f1670-9632-4227-a01f-8fbc44aa028e)

2. Reaction-diffusion models: 
For simulating different mechanisms for cellular polarization discussed in the article (SubPB, LEGI, Wave-pinning and Turing) upon stimulation with spatial-temporal signal
gradients. 

 ```python

    # file 'main.py' has class ReactionDiffusion1D that has all attributes for numerically solving the (two variable or three variable) partial differential equation.
    # in 1D space. Any new model can be added to 'models_repo_rd.py'.

    '''
    This piece of example code generates the Kymograph of SubPB upon transient gradient stimulation (See the image below).
    Copy and paste it below 'if __name__ == '__main__':' in the 'main.py' code. For more details of functions please see the 'main.py' file. 
    '''

    lbo = Periodic() # periodic boundary condition
    model = SubPB();threshold_stimu=0.013 # model selction
    initial_condition = around_steadystate_2vars() # intial value selection    
    stimulus_type=single_gradient() # stimulus selction
    stimulus_type.stimu_strength=2*threshold_stimu 
    stimulus_type.t_end=70
    
    ## create class instance for the given model using initial, boundary and stimulus conditions.
    rd = ReactionDiffusion1D(model, initial_condition, lbo, stimulus_type)

    ## deterministic or stochastic
    add_noise=None
    rd.add_noise=add_noise # assigning the condition as a 'rd' class instance
    
    ## integrate the system and find solution
    ## If add_noise=True, stochastic simulation. If add_noise=None, deterministic simulation.
    sol_det, sol_stocha=rd.simulate()
    
 
    if sol_det is None:
        out=sol_stocha[:,:rd.N].T             
    else:
        out=sol_det.y[:rd.N]
    
    ## slicing the solution array
    tout=rd.t_eval[::int(1/rd.dt)]
    out = out[:,::int(1/rd.dt)] # output u activity from each model
    
    '''
    plotting functions
    '''
    tt=[stimulus_type.t_beg,stimulus_type.t_end] # specify thr stimuls time points
    rd.plot_profile(rd.Stimulus.T,tt)
    rd.plot_kymo(out,tt)     
    rd.plot_timeseries(out,bins=[10,0],tt=tt)
 ```
![reaction_diffusion](https://github.com/akhileshpnn/SubPB-mechanism/assets/41164857/1b94e60e-460f-45d8-82cb-c2cca18bd02a)





  
