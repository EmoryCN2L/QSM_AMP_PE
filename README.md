# Robust Quantitative Susceptibility Mapping via Approximate Message Passing with Parameter Estimation
We propose a probabilistic Bayesian approach for QSM with built-in parameter estimation, and incorporate the nonlinear formulation of the dipole inversion to achieve a robust recovery of the susceptibility maps. 
In addition, we propose a morphology mask for the image wavelet coefficients to incorporate anatomical structural information into the reconstruction.

* If you use this code and find it helpful, please cite the above paper. Thanks :smile:
```
@ARTICLE{QSM_AMP_PE:2022,
    author    = {Shuai Huang and James J. Lah and Jason W. Allen and Deqiang Qiu},
    title     = {Robust Quantitative Susceptibility Mapping via Approximate Message Passing with Parameter Estimation},
    journal   = {arXiv preprint},
    volume    = {arXiv:2207.14709},
    year      = {2022},
    url       = {https://arxiv.org/abs/2207.14709},
}
```

## Summary
```
    ./amp_pe_qsm	-- This folder contains MATLAB files to perform nonlinear dipole inversion using the AMP-PE approach.
    ./functions 	-- This folder contains MATLAB files to perform pre-processing steps in the QSM pipeline.
```

Detailed comments and instructions are provided in the file `qsm_rec_src.m`. 

Please update the locations of the phase and magnitude images accordingly, then run it in `MATLAB`.

The proposed approach first perform a preliminary reconstruction of the susceptibility map using a single Gaussian to mdoel the noise. Based on the preliminary reconstruction, the 


The QSM pipeline described in `qsm_rec_src.m` mainly contains:

* 1) Set up the parameters
* 2) Generate a region-of-interest (ROI) mask using the bet tool
* 3) Use 3D best-path phase unwrapping to unwrap the phase images
* 4) Use the projection onto dipole fields (PDF) to remove the background field from unwrapped phase images
* 5) Perform the preliminary reconstruction of the susceptibility map via AMP-PE using a single Gaussian to mdoel the noise.
* 6) Perform the final reconstruction via AMP-PE using a two-component Gaussian mixture to model the noise. The mixture weights are determined based on the preliminary reconstruction and fixed to avoid over-estimation.