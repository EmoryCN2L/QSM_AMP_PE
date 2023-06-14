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

This package requires the Wavelet Toolbox of MATLAB. Please install it before you run the code.

Detailed comments and instructions are provided in the files `qsm_single_echo.m`, `qsm_multi_echo_combined.m` and `qsm_multi_echo.m`. 

* For a single-echo acquisition, we use `qsm_single_echo.m` to perform QSM.

* For a multi-echo acquisition, there are two ways to perform QSM:
1) `qsm_multi_echo_combined.m`: The local fields (in Hz) are computed for each echo time and then averaged. The averaged local field is converted to a simulated phase image (in radian) by multiplying it with a chosen TE. The susceptibility map is then recovered from the simulated phase image.
2) `qsm_multi_echo.m`: AMP-PE naturally supports recovering the susceptibility map from multi-echo phase images. In this case, the raw phase images go through phase unwrapping and background field removal to produce the processed multi-echo phase images. The susceptibility map is then recovered using from the processed multi-echo phase image.

Please try the first way first. The first way is faster since it uses data from only one simulated echo (with the averaged local field). Usually the first way would be good enough. The second way is slower since it uses data from multiple echoes. The second way performs a bit better than the first way.

## Parameter setting

**Wavelet basis**: Choosing an apropriate wavelet basis is important. We usually use the `db1` or `db2` basis for QSM. The best choice typically depends on whether it is a straight or oblique scan.  Generally speaking,

 1) When the B0 direction is [0 0 1] or close to [0 0 1], the db1 wavelet basis would perform as well as or better than the db2 basis.
 2) When the B0 direction is tilted and far from [0 0 1], sometimes the db1 basis leads to pixelation artifacts due to the approximation of dipole kernel.

It is thus recommend to run experiments with both the db1 and db2 bases, and see which one performs better for your specific dataset. 

In addition to updating the locations of phase and magnitude images in "PhaseImagLoc" and "MagImgLoc", please also make sure to update the parameters in the section "set up the parameters". At the very least, the following parameters must be checked and updated accordingly. You are welcome to experiment with the rest parameters.

```
>> LN_num
>> output_dir
>> PhaseImagLoc
>> MagImgLoc
>> voxel_size
>> TE
>> simulated_TE
>> B0
```

Other important parameters worth checking:

```
>> wave_pec % the percentage threshold for generating the morphology mask to preserve anatomical structures, usually around 80%~85% would suffice
>> wave_idx % try both db1 and db2 wavelet bases by setting wave_idx to 1 and 2 respectively
>> nlevel % the wavelet basis level
>> damp_rate_sig % the damping/learning rate to estimate the signal
>> damp_rate_par % the dampling/learning rate to estimate the parameters
```

## Overview of the pipeline
The proposed approach first perform a preliminary reconstruction of the susceptibility map using a single Gaussian to mdoel the noise. Based on the preliminary reconstructin, we can estimate the amount of noise outliers. We then perform the final reconstruction using a two-component Gaussian mixture to model the noise, where the second component is for modelling the noise outliers. Detailed discussions are given in the paper.


The QSM pipeline described in `qsm_single_echo.m`, `qsm_multi_echo_combined.m` and `qsm_multi_echo.m` mainly contains:

1) Set up the parameters
2) Generate a region-of-interest (ROI) mask using the bet tool
3) Use 3D best-path phase unwrapping to unwrap the phase images
4) Use the projection onto dipole fields (PDF) to remove the background field from unwrapped phase images
5) Perform the preliminary reconstruction of the susceptibility map via AMP-PE using a single Gaussian to mdoel the noise.
6) Perform the final reconstruction via AMP-PE using a two-component Gaussian mixture to model the noise. The mixture weights are determined based on the preliminary reconstruction and fixed to avoid over-estimation.
