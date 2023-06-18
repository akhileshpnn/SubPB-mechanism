# SubPB-mechanism


## <div align="center"><b><a href="README.md">SubPB</a></b></div>

<div align="center">

🚩 [**Paper**](#-Paper) **|** 🏰 [**Models and excecution**](#-Model-Zoo)**|** [**Demo**](#-Demo)


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
1. One-dimensional projection models 
Recaluates evolution of state-space tracjectories with quasi-potential landscapes and their projections in Figure 1 and Figure 2 (Including their supplementary figures)

 ```bash
    # Simulate one dimensional projection model
     file main.py has class Model1D that has all attributes for integrating the system, estimating quasi-potential landscape and plotting them.
  ```

2. Reaction diffusion models 
Recalculates all feature comparison done with polarization mechanims outlined in Figure 3 and supplementary.

 

---
## &#x1F308; Demo

  ### 1. STORM Denoising
  <p align="center">
  <img src="assets/storm_vis.png" width='800'>
  </p>

  ### 2. Localization and reconstruction of STORM
  <p align="center">
  <img src="assets/storm_rec.png" width='800'>
  </p>

  ### 3. Calcium Denoising at 30Hz
  <p align="center">
  <img src="assets/cad_30hz.png" width='800'>
  </p>

  ### 4. Calcium Denoising at 0.3Hz
  <p align="center">
  <img src="assets/cad_0.3hz.png" width='800'>
  </p>

  ### 4. Videos
  1. Three-dimensional visualization of large neuronal populations in a 510 × 510 × 548 μm volume (100 planes, 0.3-Hz volume rate). Left, low-SNR raw volume. Right, the same volume denoised with SRDTrans.
  <p align="center">
    <img src="assets/Supplementary Video_5.gif" width='800'>
    </p>
  2. Comparison of denoising performance of DeepCAD, CNN, and SRDTrans at 0.3Hz calcium imaging data.
  <p align="center">
    <img src="assets/Supplementary Video_4.gif" width='800'>
    </p>
  3. Validation experiments on synthetic MNIST datasets of different moving speed [0.5,0.5,10] (pixel/s). Faster movement of digits in MNIST means a lager discrepancy gap between adjacent frames. 

  <p align="center">
    <img src="assets/Supplementary Video_3.gif" width='800'>
    </p>

[def]: #-demos-videos
