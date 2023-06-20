# SubPB-mechanism


## <div align="center"><b><a href="README.md">SubPB</a></b></div>

<div align="center">

🚩 [**Paper**](#-Paper) **|** ⚡ **Quick Guide** **|** 


</div>

---


## 🚩 Paper

This repository is for the models described in the paper

"Non-asymptotic transients away from steady states determine cellular responsiveness to dynamic spatial-temporal signals
" [[biorRXiv]]([https://www.biorxiv.org/content/10.1101/2023.06.01.543361v1](https://www.biorxiv.org/content/10.1101/2023.02.03.526969v1)) 

## 🔧 Dependencies 
  - Python >= 3.6 
  - Numpy == 1.22.3
  - Scipy == 1.7.3

## ⚡ Quick Guide

The repository has two folders,
1. One-dimensional projection models: 
Recaluates evolution of state-space tracjectories with quasi-potential landscapes and their projections in Figure1 and Figure2 (including corresponding supplementary figures) of the paper.

```python
# file 'main.py' has class Model1D that has all attributes for integrating the system, estimating quasi-potential landscape and plotting them.
# example with one-dimensional projection of Wave-pinning model. Variable notations same as in text. Any new model can be added to 'models_repo.py'. Currently works only for models with two variables.

model = wavepinning_1d()
stimulus_type=stimulus_ramp()
# stimulus_type=stimulus_single()
model.stimulus_type=stimulus_type

## simulation conditions
initial_condition = [0.1,0.05]

# model parameters 
sLmax=0.02 # maximum stimulus strength on the left side of the cell (sleft)
sRmax=0.0 # maximum stimulus strength on the right side of the cell (sright)
total=2.21 # region II, criticality  

input_params=[total,sLmax,sRmax]

# integrating the ODEs

m1d=Model1D(model,stimulus_type) # generating one-dimensional projection model instance
m1d.input_params=input_params # contains bifurcation parameter and stimulus amplitudes
uL,uR=m1d.solve_timeseries()

# Estimating quasi-potential
qpl=QuasiPotentialLandscape(time_point,model,input_params)
grid_pot,Pt=qpl.find_potential() 
Q=-np.log(Pt)
```

2. Reaction-diffusion models: 
Reaction-diffusion simulation for recalculating all features compared between polarization mechanisms outlined in Figure 3 and its supplementary.
Example script to obtain Kymograph in Figure1C
 ```python

    # file 'main.py' has class ReactionDiffusion1D that has all attributes for numerically solving the partial differential equation.
    # Any new model can be added to 'models_repo.py'. Currently works only for models with two and three variables.

    lbo = Periodic() # periodic boundary condition
    
    ## select the model to simulate and corresponding threshold stimulus strength (from Figure 3A)
    model = SubPB();threshold_stimu=0.012       
    
    initial_condition = around_steadystate_2vars() # initial condition around homogeneous steady state  
    
    ## select type of stimulus
    stimulus_type=single_gradient_transient()
    
    stimulus_type.stimu_strength=0.02 # set stimulus strength.
    
    ## create class instance for the given model
    rd = ReactionDiffusion1D(model, initial_condition, lbo, stimulus_type)
    stimulus_type.tF=rd.tF
    
    add_noise=None # for performing stochastic simulation, set add_noise=True
    rd.add_noise=add_noise
    
    ## integrate the system and find solution
    ## If add_noise=True, stochastic simulation. If add_noise=None, deterministic simulation.
    sol_det, sol_stocha=rd.simulate()
    
    if sol_det is None:
        out=sol_stocha[:,:rd.N].T             
    else:
        out=sol_det.y[:rd.N]

    tout=rd.t_eval[::int(1/rd.dt)]
    out = out[:,::int(1/rd.dt)]
    
    tt=[10,70]
    rd.plot_profile(rd.Stimulus.T,tt)
    rd.plot_kymo(out,tt)     
    rd.plot_timeseries(out,[10,0],tt)

---
## &#x1F308; Demo

  
