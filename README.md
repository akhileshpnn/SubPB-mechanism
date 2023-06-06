# SubPB-mechanism

<p align="center">
  <img src="assets/SRDTrans.gif" width='600'>
</p> 

## <div align="center"><b><a href="README.md">SRDTrans</a></b></div>

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




## ⚡ Quick Inference
1. Pretrained model

    Download the pretrained model

2. Data preparation 

    Please delete the "_\_init__.py" file used for occupancy. Then, you can download the demo data(.tif file) and put the clean data into datasets/clean/.

3. Test

  ```bash
    # Simulated Calcium imaging dataat 0.3hz
    python -u test.py --denoise_model cad_03hz --patch_x 128 --patch_t 128 --GPU 0 --ckp_idx [test_idx] --datasets_folder noisy --test_datasize 1000 --datasets_path datasets/ --clean_path datasets/clean/clean.tif
  ```

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
