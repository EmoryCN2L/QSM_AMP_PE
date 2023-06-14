% add the path
addpath(genpath('./amp_pe_qsm/'))
addpath(genpath('./functions/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ***Note***: The 3D best-path phase unwrapping algorithm creates temporary files in the current working directory. To avoid overwriting the temporary files when multiple QSM reconstructions are run in parallel, please *create a different working directory* for each QSM reconstsruction, and run the experiments in each individual directory accordingly.

LN_num = 10;    % the maximum number computing threads
LN = maxNumCompThreads( LN_num );   % set the number of computing threads
output_dir = './results/';    % the directory where the reconstruction is saved
system(['mkdir -p ', output_dir])   % create the directory incase it was not created

% the locations of the multi-echo phase images and magnitude images
% the images should be saved in nifti format, make sure that the header of the nifti file is correct
% please normalize the phase image if necessary so that it spans a complete cycle: [-pi, pi] or [0,2*pi]
PhaseImagLoc = 'Please_put_the_phase_image_location_here';
MagImagLoc = 'Please_put_the_magnitude_image_location_here';
rotMat = niftiQFormRot(PhaseImagLoc);
foo = niftiLoadNii(PhaseImagLoc);

sx = size(foo.img,1);
sy = size(foo.img,2);
sz = size(foo.img,3);
Ne = size(foo.img,4);

mat_sz = [sx sy sz];

voxel_size = [0.6875 0.6875 0.7];   % Ideally, the voxel size should be isotropic
TE = [7.32 16 24.7 33.39]*1e-3; % echo time for each GRE image, in second
B0_dir = rotMat(3,:);       % main magnetic field direction, [x,y,z]. It is recommended to perform straight acquisiton where B0_dir is [0 0 1]. Oblique acquisition often brings additional artifacts in the recovered susceptiblity map due to the finite approximation of the dipole kernel...
B0 = 3;                 % magnetic field strength, in Tesla
gyro_ratio = 42.58;     % gyromagnetic ratio
pdf_method = 1;         % option 1: use our own PDF implementation; option 2: use the PDF method from the MEDI toolbox. Option 1 performs zero-padding of the image before Fourier transform and produces better results, though a bit slower. Option 2 sometimes produce strange banding artifacts...
pdf_ite = 100;          % the number of PDF iterations, which should be be at least 100 when option 1 is used and 200 when option 2 is used. It can be increased to improve the results a little bit
pdf_tol = 1e-3;         % the convergence threshold for PDF
erosion_width_1 = 2;    % the first erosion width to remove low-SNR voxels
erosion_width_2 = 2;    % the second erosion width to the voxels at the boundary of ROI after PDF



nlevel = 3;     % the number of levels in the wavelet tranform, ususally 3-4 is enough
wave_idx = 1;   % the order of the Daubechies wavelets, the db2 wavelet basis is recommended for QSM. You can also try the db1 wavelet basis by setting wave_idx to 1. Higher-order wavelet bases could capture more high-frequency info from the streaking artifacts and are not recommended.
wave_pec = 0.85;    % the percentage threshold (with respect to the l1-norm of wavelet coefficients) used to generate the morphology mask for wavelet coefficients, so that anatomical information could be incorporated into reconstruction. the higherer the wave_pec, the more high-frequency information is retained.
max_linearization_ite = 25;  % the maximum number of linearization steps to linearize the complex exponential measurement mdoel

l2beta_reg = 2e-2;      % regularization parameter for the l2-norm minimization solution that is used to initialize the distribution parameters only

damp_rate_sig = 0.01;   % damp rate for the signal update. If the algorithm does not converge, try decreasing the damp_rate
damp_rate_par = 0.1;    % damp rate for the parameter estimation. If the algorithm does not converge, try decreasing the damp_rate
gamp_cvg_thd = 1e-6;    % convergence rate of the AMP updates
max_pe_spar_ite = 5;    % maximum number of iterations for the linear AMP updates
max_pe_est_ite = 5;     % maximum number of iterations for the parameter estimation updates

%%%%%%%%%%%%%%%%%%%%
%% Pre-processing %%
%%%%%%%%%%%%%%%%%%%%

foo = niftiLoadNii(MagImagLoc);
magnitude_image = double(foo.img);

% make sure the phase image spans a complete cycle: [-pi, pi] or [0,2*pi]
foo = niftiLoadNii(PhaseImagLoc);
phase_image = double(foo.img);

iField = magnitude_image.*exp(1i*phase_image);

% Estimate the frequency offset in each of the voxel using a complex fitting
[iFreq_raw N_std relres initial_phase] = Fit_ppm_complex_TE(iField,TE);

% Computed a *weighted* magnitude image, it will be used to generate the morphology mask for wavelet coefficients
iMagWtd = sqrt(sum(abs(iField).^2,4));
% The default threshold of BET is 0.5, if it does not produce a good brain mask, try changing the threshold 
mask = BET(iMagWtd,mat_sz,voxel_size);
%mask = BET(iMagWtd,mat_sz,voxel_size, 0.5);

% Apply a small erosion to remove low SNR voxels
mask = MaskErode(mask, mat_sz, voxel_size, erosion_width_1);

% Remove the initial phase before phase unwrapping
iField_phase = zeros(size(iField));
for (i=1:size(iField,4))
    iField_phase(:,:,:,i) = angle(iField(:,:,:,i)) - initial_phase;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform phase unwrapping %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iField_phase_unwrap = zeros(size(iField_phase));
for (i=1:size(iField_phase,4))
    iField_phase_unwrap(:,:,:,i) = UnwrapPhase_3DBestPath_voxelsize(iField_phase(:,:,:,i), mask, mat_sz, voxel_size);
end
iFreq = iField_phase_unwrap;


%%%%%%%%%%%%%%%%%
%% perform PDF %%
%%%%%%%%%%%%%%%%%

% option 1: use PDF from our implementation
% option 2: use PDF from the MEDI toolbox
% the number of iterations should be at least 100 when option 1 is used and 200 when option 2 is used.

if (pdf_method == 1)

    % create new nifti file name
    foo.hdr.dime.dim(1)=3;
    foo.hdr.dime.dim(5)=1;
    niftiSaveNii(foo.hdr, mask, strcat(output_dir,'mask_initial.nii'));
    opt.writeSF = 1;
    opt.B0 = B0;
    opt.num_iter = pdf_ite;
    opt.plane = rotMat;

    PDF_echo = zeros(size(iField_phase_unwrap));
    for (i=1:size(iField_phase_unwrap,4))
        niftiSaveNii(foo.hdr, iField_phase_unwrap(:,:,:,i), strcat(output_dir, 'unwrapped_phase_e',num2str(i),'.nii'));
        opt.TE = TE(i);
        Dipole_fitting(strcat(output_dir, 'unwrapped_phase_e',num2str(i),'.nii'), strcat(output_dir,'mask_initial.nii'), strcat(output_dir, 'BG_removed_phase_e',num2str(i),'.img'), opt);
        goo = niftiLoadNii(strcat(output_dir, 'BG_removed_phase_e', num2str(i), '.img'));
        PDF_echo(:,:,:,i) = goo.img;
    end

else

    PDF_echo = zeros(size(iField_phase_unwrap));
    for (i=1:size(iField_phase_unwrap,4))
        PDF_echo(:,:,:,i) = PDF(iField_phase_unwrap(:,:,:,i), N_std, mask, mat_sz, voxel_size, B0_dir, pdf_tol, pdf_ite);
    end

end

mask_ori = mask;
% Apply a second erosion to remove inaccurate voxels close to the boundary
mask = MaskErode(mask, mat_sz, voxel_size, erosion_width_2);

save(strcat(output_dir, 'mask.mat'), 'mask')

%%%%%%%%%%%%%%%%%%%%
%% QSM via AMP-PE %%
%%%%%%%%%%%%%%%%%%%%

% Read the measurements within the ROI defined by the mask
phase_image_new = zeros(sum(mask,'all'),length(TE));
for (i=1:size(phase_image_new,2))
    phase_image_tmp = PDF_echo(:,:,:,i);
    phase_image_new(:,i) = phase_image_tmp(mask==1);
end
phase_image_ori = phase_image;
phase_image = phase_image_new;

% multi-echo magnitude image is used to weight the complex expnential measurements
% every column of magnitude_image only contains voxel values within the mask and corresponds to an echo time
magnitude_image_new = zeros(sum(mask,'all'),length(TE));
for (i=1:size(magnitude_image_new,2))
    magnitude_image_tmp = magnitude_image(:,:,:,i);
    magnitude_image_new(:,i) = magnitude_image_tmp(mask==1);
end
magnitude_image_ori = magnitude_image;
magnitude_image = magnitude_image_new;

% Generate the dipole kernel
dipoleFT = GenerateDipoleFT3Drot(mat_sz, voxel_size, rotMat);
% Generate the measurement operator
A = SI_operator_withmask(dipoleFT, mat_sz, mask);

% use wavelet transform to reveal the sparse representation of an image 
X0 = zeros(sx, sy, sz);

% construct the wavelet transform
dwtmode('per');
switch wave_idx
case 1
    C1=wavedec3(X0,nlevel,'db1'); 
    ncoef1=length(C1.dec);
case 2
    C2=wavedec3(X0,nlevel,'db2'); 
    ncoef2=length(C2.dec);
case 3
    C3=wavedec3(X0,nlevel,'db3'); 
    ncoef3=length(C3.dec);
case 4
    C4=wavedec3(X0,nlevel,'db4'); 
    ncoef4=length(C4.dec);
case 5
    C5=wavedec3(X0,nlevel,'db5'); 
    ncoef5=length(C5.dec);
case 6
    C6=wavedec3(X0,nlevel,'db6'); 
    ncoef6=length(C6.dec);
case 7
    C7=wavedec3(X0,nlevel,'db7'); 
    ncoef7=length(C7.dec);
otherwise
    C8=wavedec3(X0,nlevel,'db8'); 
    ncoef8=length(C8.dec);
end


switch wave_idx
case 1
    Psi = @(x) [wavedec3(x,nlevel,'db1')']; 
    Psit = @(x) (waverec3(x));
    wav_vessel = C1; 
case 2
    Psi = @(x) [wavedec3(x,nlevel,'db2')']; 
    Psit = @(x) (waverec3(x));
    wav_vessel = C2; 
case 3
    Psi = @(x) [wavedec3(x,nlevel,'db3')']; 
    Psit = @(x) (waverec3(x));
    wav_vessel = C3; 
case 4
    Psi = @(x) [wavedec3(x,nlevel,'db4')']; 
    Psit = @(x) (waverec3(x));
    wav_vessel = C4; 
case 5
    Psi = @(x) [wavedec3(x,nlevel,'db5')']; 
    Psit = @(x) (waverec3(x));
    wav_vessel = C5; 
case 6
    Psi = @(x) [wavedec3(x,nlevel,'db6')']; 
    Psit = @(x) (waverec3(x));
    wav_vessel = C6; 
case 7
    Psi = @(x) [wavedec3(x,nlevel,'db7')'];
    Psit = @(x) (waverec3(x));
    wav_vessel = C7;
otherwise
    Psi = @(x) [wavedec3(x,nlevel,'db8')'];
    Psit = @(x) (waverec3(x));
    wav_vessel = C8;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use the truncated threshold method to compute a solution for the initialization of distribution parameters %%
%% Note that this is NOT the initialization of the susceptiblity map, which is acutally an all-zero vector    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use the averaged local field to compute the initialization for parameters only
f_image = zeros(size(phase_image(:,1)));
for (i=1:size(phase_image,2))
    f_image = f_image + phase_image(:,i)/(gyro_ratio*TE(i)*B0*2*pi);
end
f_image = f_image/size(phase_image,2);

NN=size(mask);
spatial_res = voxel_size;
kernel = dipole_kernel_angulated( NN, spatial_res, B0_dir );
phs_tissue = zeros(size(mask));
phs_tissue(mask==1) = f_image;
l2beta = l2beta_reg;    % regularization parameter
chi_L2 = chiL2( phs_tissue, mask, kernel, l2beta, NN );

% this initialization is only used to estimate the parameters, not to initialize the susceptibility map
X_init_par = real(chi_L2);
X_init_par(mask==0) = 0;

% construct 3D wavelet operator to extract coefficients and create MATLAB structures
E_coef = @(in) extract_3d_wav_coef(in);
C_struct = @(in) construct_3d_wav_struct(in, wav_vessel);

% Generate a morphology mask for the wavelet coefficients to incorporate anatomical structural information
% iMagWtd is the weighted magnitude image
iMagWtd_with_mask = iMagWtd;
iMagWtd_with_mask(mask==0) = 0;
iMagWtd_with_mask_wavelet = E_coef(Psi(iMagWtd_with_mask));
iMagWtd_with_mask_wavelet_abs_sort = sort(abs(iMagWtd_with_mask_wavelet), 'descend');
iMagWtd_with_mask_wavelet_abs_sort_cumsum = cumsum(iMagWtd_with_mask_wavelet_abs_sort)/sum(iMagWtd_with_mask_wavelet_abs_sort);

iMagWtd_with_mask_wavelet_abs_sort_thd = iMagWtd_with_mask_wavelet_abs_sort(length(iMagWtd_with_mask_wavelet_abs_sort_cumsum(iMagWtd_with_mask_wavelet_abs_sort_cumsum<=wave_pec)));
wave_mask = zeros(length(iMagWtd_with_mask_wavelet_abs_sort),1);
wave_mask(abs(iMagWtd_with_mask_wavelet)>iMagWtd_with_mask_wavelet_abs_sort_thd) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preliminary reconstruction via AMP-PE with one Gaussian component %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the gamp parameters
gamp_par.damp_rate = damp_rate_sig;  % the damp rate for the susceptibility map updates
gamp_par.max_pe_spar_ite = max_pe_spar_ite;   % the number of AMP updates in every linearized model step
gamp_par.max_pe_est_ite = max_pe_est_ite;    % the number of parameter estimation iterations

gamp_par.cvg_thd = gamp_cvg_thd;    % convergence threshold
gamp_par.kappa = damp_rate_par;   % the damp rate for the parameter estimation updates
gamp_par.sx = sx;
gamp_par.sy = sy;
gamp_par.sz = sz;

% add the ROI mask and the morphology mask
gamp_par.mask = mask;
gamp_par.wave_mask = wave_mask;


weight_vect = magnitude_image;  % the multi-echo magnitude images are used as the weights

% gamp recovery with magnitude info 

X_init_par_psi_struct = Psi(X_init_par);
X_init_par_psi = extract_3d_wav_coef(X_init_par_psi_struct);

wav_coef_len = length(X_init_par_psi);

M = mat_sz(1)*mat_sz(2)*mat_sz(3);
N = wav_coef_len;

A_wav_single_3d = A_wav_single_3d_LinTrans(M,N,mat_sz,wav_coef_len,Psit,Psi,C_struct,E_coef);

X_init = zeros(size(X_init_par));   % the susceptiblity map is initialized with all zeros

for (iter = 1:max_linearization_ite)

    % In each inner iteration, we need to update the linearized measurement model

    if (iter==1)
        tau_w_1 = 1e-12; % the initialization of noise variance
        output_par.tau_w_1 = tau_w_1;   % set the output channel parameter
    else
        tau_w_1 = output_par_new.tau_w_1;   % initialize the noise variance use the previous estimation
        output_par.tau_w_1 = tau_w_1;   % set the output channel parameter
    end


    if (iter==1)
        gamp_par.x_hat_meas = X_init;  % initialize the susceptibility map with all zeros
        gamp_par.tau_x_meas = var(X_init_par(:));   % the variance should not be initizalized with zeros
        gamp_par.s_hat_meas_1 = zeros(size(phase_image));

        gamp_par.x_hat_psi = zeros(size(X_init_par_psi));   % initialize the wavelet coefficients with all zeros
        gamp_par.p_hat_psi = zeros(size(X_init_par));
        gamp_par.tau_x_hat_psi = var(X_init_par_psi(:));    % the variance should not be initizalized with zeros
        gamp_par.tau_p_psi = A_wav_single_3d.multSq(gamp_par.tau_x_hat_psi);    % the variance should not be initizalized with zeros

        input_par.lambda_x_hat_psi = 1/sqrt(var(abs(X_init_par_psi))/2);    % set the input channel parameter

    else
        gamp_par.x_hat_meas = res.x_hat_meas;
        gamp_par.tau_x_meas = res.tau_x_meas;
        gamp_par.s_hat_meas_1 = res.s_hat_meas_1;

        gamp_par.x_hat_psi = res.x_hat_psi;
        gamp_par.p_hat_psi = res.p_hat_psi;
        gamp_par.tau_x_hat_psi = res.tau_x_hat_psi;
        gamp_par.tau_p_psi = res.tau_p_psi;

        input_par.lambda_x_hat_psi = input_par_new.lambda_x_hat_psi;    % set the input channel parameter

    end


    % Calculated the linearized measurement model
    M_row = size(phase_image,1);
    M_col = size(phase_image,2);
    M=M_row*M_col;
    N=mat_sz(1)*mat_sz(2)*mat_sz(3);

    mut_cst = gyro_ratio*B0*2*pi*TE;

    A_X_init = zeros(size(phase_image));
    for (i=1:size(A_X_init,2))
        A_X_init(:,i) = real(A.times(X_init)) * mut_cst(i);
    end

    % the first-order derivative
    der_1st = 1i*exp(1i*A_X_init);

    % generate the measurement operator in the linearized model
    A_qsm_weighted_nw_combine_real = A_qsm_weighted_nw_combine_real_LinTrans(M, N, mat_sz, der_1st.*weight_vect, A, mut_cst, M_row, M_col);

    % update the measurements in the linearized model
    phase_image_updated = weight_vect .* (der_1st.*A_X_init+exp(1i*phase_image)-exp(1i*A_X_init));

    % use AMP-PE to recover the susceptibility under the linearized model
    [res, input_par_new, output_par_new] = amp_pe_mri_qsm_awgn(A_qsm_weighted_nw_combine_real, A_wav_single_3d, phase_image_updated, gamp_par, input_par, output_par);

    % the susceptibility initialization for the next iteration
    X_init_pre = X_init;
    X_init = res.x_hat_meas;

    nonlinear_cvg_rate = norm(X_init(:)-X_init_pre(:), 'fro') / norm(X_init(:), 'fro');
    fprintf('Nonlinear convergence rate: %d\n', nonlinear_cvg_rate)
    A_X_init = zeros(size(phase_image));
    for (i=1:size(A_X_init,2))
        A_X_init(:,i) = real(A.times(X_init)) * mut_cst(i);
    end
    fprintf('Residue: %d\n', norm(weight_vect.*(exp(1i*A_X_init)-exp(1i*phase_image)), 'fro'))

end

% the recovered susceptibility map in the preliminary step (step1) with a single Gaussian to model the noise
QSM_step1 = res.x_hat_meas;
QSM_step1(mask==0) = 0;

save(strcat(output_dir, 'QSM_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step1.mat'), 'QSM_step1', '-v7.3')
save(strcat(output_dir, 'res_wave_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step1.mat'), 'res', '-v7.3')
save(strcat(output_dir, 'input_par_wave_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step1.mat'), 'input_par_new', '-v7.3')
save(strcat(output_dir, 'output_par_wave_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step1.mat'), 'output_par_new', '-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the mixture weights based on the preliminary reconstruction %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_init = res.x_hat_meas;    % the initialization for the second step of AMP-PE
A_X_init = zeros(size(phase_image));
for (i=1:size(A_X_init,2))
    A_X_init(:,i) = real(A.times(X_init)) * mut_cst(i);
end

residue_mat = zeros([sx sy sz length(TE)]);
for (i=1:size(residue_mat,4))
    residue_mat_tmp = zeros(mat_sz);
    residue_mat_tmp(mask==1) = weight_vect(:,i).*(exp(1i*A_X_init(:,i))-exp(1i*phase_image(:,i)));
    residue_mat(:,:,:,i) = residue_mat_tmp;
end

residue_mat_new = zeros(size(phase_image));
for (i=1:size(residue_mat,4))
    residue_mat_tmp = residue_mat(:,:,:,i);
    residue_mat_new(:,i) = residue_mat_tmp(mask==1);
end
residue_mat_seq = residue_mat_new(:);

% assume mean is zero
residue_mat_seq_abs = abs(residue_mat_seq);
residue_std_abs = sqrt(mean((residue_mat_seq_abs).^2));
% the estimated mixture weight gamma_est (for the noise outliers) will be fixed in the final reconstruction to avoid over-estimation
gamma_est = length(residue_mat_seq(residue_mat_seq_abs>3*residue_std_abs))/length(residue_mat_seq);
% the estimated mixture variance psi_est for the noise outliers
psi_est = var(residue_mat_seq(residue_mat_seq_abs>3*residue_std_abs));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final reconstruction via AMP-PE with 2 Gaussian mixtures %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the output channel parameters need to be changed completely

for (iter = 1:max_linearization_ite)


    if (iter==1)
        theta_output = 0;  % the means of the zero-mean Gaussian mixtures
        phi_output = output_par_new.tau_w_1;    % the variances of the Gaussian mixtures, this will be updated
        omega_output = 1;  % the weights of the Gaussian mixtures that are not reserved for noise outliers, since there is only one Gaussian component left to model the non-outliers, so it is 1 here
        num_c_output = 1;  % the number of Gaussian mixtures that are not reserved for noise outliers
        gamma_output = gamma_est;  % the weight of the outlier distribution  (a zero-mean Gaussian), this is fixed 
        psi_output = psi_est;    % the variance of the outlier distribution  (a zero-mean Gaussian), this will be updated

    else
        theta_output = output_par_new.theta_output;  % the means of the zero-mean Gaussian mixtures
        phi_output = output_par_new.phi_output;    % the variances of the Gaussian mixtures, this will be updated
        omega_output = output_par_new.omega_output;  % the weights of the Gaussian mixtures
        num_c_output = output_par_new.num_c_output;  % the number of Gaussian mixtures
        gamma_output = output_par_new.gamma_output;  % the weight of the outlier distribution (a zero-mean Gaussian), this is fixed
        psi_output = output_par_new.psi_output;    % the variance of the outlier distribution (a zero-mean Gaussian), this will be updated

    end

    output_par.theta_output     = theta_output;  % the means of the Gaussian mixtures
    output_par.phi_output       = phi_output;    % the variances of the Gaussian mixtures
    output_par.omega_output     = omega_output;  % the weights of the Gaussian mixtures
    output_par.num_c_output     = num_c_output;  % the number of Gaussian mixtures
    output_par.gamma_output     = gamma_output;  % the weight of the outlier distribution (a   zero-mean Gaussian)
    output_par.psi_output       = psi_output;    % the variance of the outlier distribution (a zero-mean Gaussian)

    input_par.lambda_x_hat_psi = input_par_new.lambda_x_hat_psi;    % input channel parameter


    gamp_par.x_hat_meas = res.x_hat_meas;
    gamp_par.tau_x_meas = res.tau_x_meas;
    gamp_par.s_hat_meas_1 = res.s_hat_meas_1;

    gamp_par.x_hat_psi = res.x_hat_psi;
    gamp_par.p_hat_psi = res.p_hat_psi;
    gamp_par.tau_x_hat_psi = res.tau_x_hat_psi;
    gamp_par.tau_p_psi = res.tau_p_psi;


    % Calculated the linearized measurement model
    M_row = size(phase_image,1);
    M_col = size(phase_image,2);
    M=M_row*M_col;
    N=mat_sz(1)*mat_sz(2)*mat_sz(3);

    mut_cst = gyro_ratio*B0*2*pi*TE;

    A_X_init = zeros(size(phase_image));
    for (i=1:size(A_X_init,2))
        A_X_init(:,i) = real(A.times(X_init)) * mut_cst(i);
    end

    % the first-order derivative
    der_1st = 1i*exp(1i*A_X_init);

    % generate the measurement operator in the linearized model
    A_qsm_weighted_nw_combine_real = A_qsm_weighted_nw_combine_real_LinTrans(M,N,mat_sz,der_1st.*weight_vect,A,mut_cst, M_row,M_col);

    % update the measurements in the linearized model
    phase_image_updated = weight_vect .* (der_1st.*A_X_init+exp(1i*phase_image)-exp(1i*A_X_init));

    % use AMP-PE to recover the susceptibility under the linearized model
    [res, input_par_new, output_par_new] = amp_pe_mri_qsm_awgn_mix(A_qsm_weighted_nw_combine_real, A_wav_single_3d, phase_image_updated, gamp_par, input_par, output_par);

    % the susceptibility initialization for the next iteration
    X_init_pre = X_init;
    X_init = res.x_hat_meas;

    nonlinear_cvg_rate = norm(X_init(:)-X_init_pre(:), 'fro') / norm(X_init(:), 'fro');
    fprintf('Nonlinear convergence rate: %d\n', nonlinear_cvg_rate)
    A_X_init = zeros(size(phase_image));
    for (i=1:size(A_X_init,2))
        A_X_init(:,i) = real(A.times(X_init)) * mut_cst(i);
    end
    fprintf('Residue: %d\n', norm(weight_vect.*(exp(1i*A_X_init)-exp(1i*phase_image)), 'fro'))

end

% the recovered susceptibility map in the final step (step2) with a two-component Gaussian mixture to model the noise
QSM_step2 = res.x_hat_meas;
QSM_step2(mask==0) = 0;

save(strcat(output_dir, 'QSM_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step2.mat'), 'QSM_step2', '-v7.3')
save(strcat(output_dir, 'res_wave_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step2.mat'), 'res', '-v7.3')
save(strcat(output_dir, 'input_par_wave_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step2.mat'), 'input_par_new', '-v7.3')
save(strcat(output_dir, 'output_par_wave_db',num2str(wave_idx),'_l',num2str(nlevel),'_',num2str(wave_pec*100),'_step2.mat'), 'output_par_new', '-v7.3')

figure; imshow3D(QSM_step2,[-0.1,0.1])

