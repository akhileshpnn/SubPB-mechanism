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
Recaluates evolution of state-space tracjectories with quasi-potential landscapes and their projections in Figure 1 and Figure 2 (including corresponding supplementary figures)

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

2. Reaction diffusion models: 
Recalculates all feature comparison done with polarization mechanims outlined in Figure 3 and supplementary.

 

---
## &#x1F308; Demo

  
