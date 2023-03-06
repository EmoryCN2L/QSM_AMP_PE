function [res, input_par, output_par] = amp_pe_mri_qsm_awgn_mix(A_qsm_weighted_nw, A_psi_mask_nw, y, gamp_par, input_par, output_par)

	% set GAMP parameters
    max_pe_spar_ite = gamp_par.max_pe_spar_ite;
    max_pe_est_ite = gamp_par.max_pe_est_ite;
	cvg_thd = gamp_par.cvg_thd;
	kappa = gamp_par.kappa;
    damp_rate = gamp_par.damp_rate;

    mask = gamp_par.mask;
    wave_mask = gamp_par.wave_mask;

    tau_x_meas = gamp_par.tau_x_meas;
    x_hat_meas = gamp_par.x_hat_meas;
    s_hat_meas_1 = gamp_par.s_hat_meas_1;

    tau_p_psi = gamp_par.tau_p_psi;
    p_hat_psi = gamp_par.p_hat_psi;
    x_hat_psi = gamp_par.x_hat_psi;

    % set input distribution parameters
    lambda_x_hat_psi = input_par.lambda_x_hat_psi;

    % set output distribution parameters
    theta_output = output_par.theta_output;  % the means of the Gaussian mixtures
    phi_output = output_par.phi_output;    % the variances of the Gaussian mixtures
    omega_output = output_par.omega_output;  % the weights of the Gaussian mixtures
    num_c_output = output_par.num_c_output;  % the number of Gaussian mixtures
    gamma_output = output_par.gamma_output;  % the weight of the outlier distribution (a zero-mean Gaussian)
    psi_output = output_par.psi_output;    % the variance of the outlier distribution (a zero-mean Gaussian)



    sx = gamp_par.sx;
    sy = gamp_par.sy;
    sz = gamp_par.sz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize tau_w_exp with tau_x_mag_meas from the warm up step %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% first finish the sparse reconstruction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for (ite_pe = 1:max_pe_spar_ite)

            tau_p_meas_1 = A_qsm_weighted_nw.multSq(tau_x_meas);
            p_hat_meas_1 = A_qsm_weighted_nw.mult(x_hat_meas) - tau_p_meas_1 * s_hat_meas_1;

            % parameter estimation
            for (ite_pe_est = 1:max_pe_est_ite)
                [omega_output, theta_output, phi_output, gamma_output, psi_output] = output_parameter_est(y(:)-p_hat_meas_1(:), tau_p_meas_1, omega_output, theta_output, phi_output, gamma_output, psi_output, kappa);
            end

            [noise_update, tau_noise_update] = output_function(y(:)-p_hat_meas_1(:), tau_p_meas_1, omega_output, theta_output, phi_output, gamma_output, psi_output);

            z_hat = y-reshape(noise_update,size(y));
            tau_z = tau_noise_update;

            tau_s_meas_1 = 1/tau_p_meas_1*(1-tau_z/tau_p_meas_1);
            s_hat_meas_1 = 1/tau_p_meas_1*(z_hat-p_hat_meas_1);

            %%%%% the first sparse info from multi echo images
            tau_r_meas_1 = 1 / A_qsm_weighted_nw.multSqTr(tau_s_meas_1);
            r_hat_meas_1 = x_hat_meas + tau_r_meas_1 * A_qsm_weighted_nw.multTr(s_hat_meas_1);

            r_hat_meas_1 = real(r_hat_meas_1);

            %% iteration in the tv block
            tau_s_psi = 1 / (tau_r_meas_1 + tau_p_psi);
            s_hat_psi = (r_hat_meas_1 - p_hat_psi) * tau_s_psi;

            tau_r_psi = 1./ A_psi_mask_nw.multSqTr(tau_s_psi);
            r_hat_psi = x_hat_psi  + tau_r_psi * A_psi_mask_nw.multTr(s_hat_psi);

            for (ite_pe_est = 1:max_pe_est_ite)
                lambda_x_hat_psi = input_parameter_est(abs(r_hat_psi), tau_r_psi, lambda_x_hat_psi, kappa);
            end

            x_hat_psi_pre = x_hat_psi;
            [x_hat_psi, tau_x_hat_psi] = input_function(r_hat_psi, tau_r_psi, lambda_x_hat_psi, wave_mask);
            %x_hat_psi = x_hat_psi_pre + damp_rate*(x_hat_psi-x_hat_psi_pre);

            tau_p_psi = A_psi_mask_nw.multSq(tau_x_hat_psi);
            p_hat_psi = A_psi_mask_nw.mult(x_hat_psi) - tau_p_psi * s_hat_psi;

            tau_x_meas_pre = tau_x_meas;
            tau_x_meas = (tau_p_psi*tau_r_meas_1) / (tau_p_psi+tau_r_meas_1);
            tau_x_meas = tau_x_meas_pre + damp_rate * (tau_x_meas-tau_x_meas_pre);

            x_hat_meas_pre = x_hat_meas;
            x_hat_meas = (tau_p_psi*r_hat_meas_1 + tau_r_meas_1*p_hat_psi) / (tau_p_psi+tau_r_meas_1);
            x_hat_meas = x_hat_meas_pre + damp_rate * (x_hat_meas-x_hat_meas_pre);

            % enforce the roi mask on the susceptibility map?
            %x_hat_meas(mask==0) = 0;

            fprintf('Para: %d\t%d\n', lambda_x_hat_psi, gamma_output)
            fprintf('Var: %d\t%d\n', phi_output(1), psi_output)

            cvg_gamp_x_hat_meas = norm(x_hat_meas(:)-x_hat_meas_pre(:), 'fro')/norm(x_hat_meas(:), 'fro');

            cvg_gamp = max([ cvg_gamp_x_hat_meas]);
            if (cvg_gamp<cvg_thd) 
                break;
            end


        fprintf('Ite %d CVG PE: %d\n', ite_pe, cvg_gamp)

        cvg_pe = max([ cvg_gamp_x_hat_meas]);
        if ((cvg_pe<cvg_thd)&&(ite_pe>2))
            break;
        end

    end

    res.x_hat_psi = x_hat_psi;
    res.tau_x_hat_psi = tau_x_hat_psi;
    res.s_hat_psi = s_hat_psi;
    res.tau_s_psi = tau_s_psi;

    res.p_hat_psi = p_hat_psi;
    res.tau_p_psi = tau_p_psi;

    res.x_hat_meas = x_hat_meas;
    res.tau_x_meas = tau_x_meas;
    res.s_hat_meas_1 = s_hat_meas_1;

    input_par.lambda_x_hat_psi = lambda_x_hat_psi;

    output_par.theta_output     = theta_output;  % the means of the Gaussian mixtures
    output_par.phi_output       = phi_output;    % the variances of the Gaussian mixtures
    output_par.omega_output     = omega_output;  % the weights of the Gaussian mixtures
    output_par.num_c_output     = num_c_output;  % the number of Gaussian mixtures
    output_par.gamma_output     = gamma_output;  % the weight of the outlier distribution (a zero-mean Gaussian)
    output_par.psi_output       = psi_output;    % the variance of the outlier distribution (a zero-mean Gaussian)


end

function [x0_hat, tau_x0] = input_function(r_hat, tau_r, lambda, wave_mask)

    r_hat_ori = r_hat;

    thresh = lambda * tau_r;
    x0_hat = max(0, abs(r_hat)-thresh) .* sign(r_hat);
    x0_hat(wave_mask==1) = r_hat_ori(wave_mask==1);

    tau_x0 = tau_r;
    %tau_x0(abs(x0_hat)==0) = 0;
    tau_x0 = tau_r * length(x0_hat(abs(x0_hat)>0)) / length(x0_hat);

end


function lambda = input_parameter_est(r_hat, tau_r, lambda, kappa)

    lambda_pre = lambda;

    % do we need some kind of normalization here?
    dim_smp=length(r_hat);
    num_cluster=1;

    block_mat = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        block_mat(:,i) = 0.5/tau_r * (tau_r*lambda(i)-r_hat).^2;
    end

    block_mat_min = [];
    if (num_cluster==1)
        block_mat_min = block_mat;
    else
        block_mat_min = max(block_mat')';
    end

    block_mat = block_mat - block_mat_min; % subtract the minimum value of each row
    block_mat_two_exp = exp(block_mat);
    block_mat_one_exp = exp(-block_mat_min);

    block_mat_erfc = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        block_mat_erfc(:,i) = erfc(sqrt(0.5/tau_r)*(tau_r*lambda(i)-r_hat));
    end

    lambda_tmp_mat_0 = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        lambda_tmp_mat_0(:,i) = lambda(i)/2 * block_mat_two_exp(:,i) .* block_mat_erfc(:,i);
    end

    % compute lambda
    % compute the first order derivative

    der_block = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        der_block(:,i) = ( sqrt(2*tau_r/pi) ./ erfcx(sqrt(0.5/tau_r)*(tau_r*lambda(i)-r_hat)) + r_hat - tau_r*lambda(i) );
    end

    fst_der_lambda = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        fst_der_lambda(:,i) = 1/lambda(i) - der_block(:,i);
    end

    scd_der_lambda = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        scd_der_lambda(:,i) = -1/lambda(i)/lambda(i) + ( tau_r + (r_hat - tau_r*lambda(i)).*der_block(:,i) ) - der_block(:,i).^2;
    end

    lambda_tmp_mat = lambda_tmp_mat_0;
    lambda_tmp_mat_sum = sum(lambda_tmp_mat, 2);
    for (i=1:num_cluster)
        lambda_tmp_mat(:,i) = lambda_tmp_mat(:,i) ./ ( lambda_tmp_mat_sum + eps);
    end

    lambda_tmp_mat_1 = lambda_tmp_mat;
    lambda_tmp_mat_2 = lambda_tmp_mat;
    %for (i=1:num_cluster)
    %    lambda_tmp_mat_1(:,i) = lambda_tmp_mat_1(:,i) .* (fst_der_lambda(:,i) - scd_der_lambda(:,i)*lambda(i));
    %    lambda_tmp_mat_2(:,i) = lambda_tmp_mat_2(:,i) .* scd_der_lambda(:,i);
    %end
    %lambda_new = - sum(lambda_tmp_mat_1) ./ (sum(lambda_tmp_mat_2) + eps);   % to avoid division by 0
    %lambda_new = lambda_new';


    lambda_new = [];
    for (i=1:num_cluster)
        lambda_tmp_mat_1_tmp = lambda_tmp_mat_1(:,i) .* fst_der_lambda(:,i);
        lambda_tmp_mat_2_tmp = lambda_tmp_mat_2(:,i) .* scd_der_lambda(:,i);
        lambda_new_tmp = 0;
        if (sum(lambda_tmp_mat_2_tmp)<0)
            lambda_new_tmp = lambda(i) - sum(lambda_tmp_mat_1_tmp)/sum(lambda_tmp_mat_2_tmp);
        else
            if (sum(lambda_tmp_mat_1_tmp)>0)
                lambda_new_tmp = lambda(i)*1.1;
            else
                lambda_new_tmp = lambda(i)*0.9;
            end
        end
        lambda_new = [lambda_new; lambda_new_tmp];
    end

    lambda_new = max(lambda_new, 1e-12);  % necessary to avoid 0 which leads to NaN

    % lambda_new could be negative, causing problem
    lambda = lambda + kappa * (lambda_new - lambda);

end



function [x_hat, tau_x] = output_function(r_hat, tau_r, omega, theta, phi, gamma, psi)

    dim_smp=length(r_hat);
    num_cluster=length(omega);

    block_mat = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        block_mat(:,i) = (1-gamma) * omega(i) * (tau_r/(phi(i)+tau_r)) * exp( -(abs(theta(i)-r_hat) / sqrt(phi(i)+tau_r)).^2 );
    end

    % compute x_hat
    block_mat_nmr_x = block_mat;
    for (i=1:num_cluster)
        block_mat_nmr_x(:,i) = block_mat_nmr_x(:,i) .* (theta(i) * tau_r + r_hat * phi(i)) / (phi(i) + tau_r);
    end

    block_mat_2 = gamma * (tau_r/(psi+tau_r)) * exp( -(abs(r_hat) / sqrt(psi+tau_r)).^2 );
    block_mat_nrm_x_2 = block_mat_2 .* (r_hat*psi / (psi+tau_r));

    nmr_x = sum(block_mat_nmr_x, 2) + block_mat_nrm_x_2;
    dnm_x = sum(block_mat, 2) + block_mat_2 ;

    x_hat = nmr_x ./ (dnm_x);

    % if dnm_x is zero, set x_hat to r_hat
    dnm_x_zero_idx = (dnm_x==0);
    x_hat(dnm_x_zero_idx) = r_hat(dnm_x_zero_idx);

    % compute tau_x
    block_mat_nmr_x_sq = block_mat;
    for (i=1:num_cluster)
        block_mat_nmr_x_sq(:,i) = block_mat_nmr_x_sq(:,i) .* (  phi(i)*tau_r/(phi(i) + tau_r) + (abs( (theta(i)*tau_r + r_hat*phi(i)) / (phi(i) + tau_r) )).^2  );
    end

    block_mat_nmr_x_sq_2 = block_mat_2 .* ( psi*tau_r/(psi+tau_r) + (abs(r_hat*psi) / (psi+tau_r)).^2 );

    nmr_x_sq = sum(block_mat_nmr_x_sq, 2) + block_mat_nmr_x_sq_2;
    dnm_x_sq = dnm_x;

    tau_x_seq = (nmr_x_sq ./ dnm_x_sq - (abs(x_hat)).^2);  % this is nonnegative in theory

    % if dnm_x is zero, set tau_x_seq to 0
    tau_x_seq(dnm_x_zero_idx) = 0;

    tau_x = mean(tau_x_seq);
    tau_x = max(tau_x, 1e-12);  % just in case

end

function [omega, theta, phi, gamma, psi] = output_parameter_est(r_hat, tau_r, omega, theta, phi, gamma, psi, kappa)

    omega_pre = omega;
    theta_pre = theta;
    phi_pre = phi;
    gamma_pre = gamma;
    psi_pre = psi;

    % do we need some kind of normalization here?
    dim_smp=length(r_hat);
    num_cluster=length(omega);  % the number of Gaussian components excluding the component for outliers, usually one is enough, but you can increase

    lambda_tmp_mat_1 = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        lambda_tmp_mat_1(:,i) = (1-gamma) * omega(i) * 1/(tau_r+phi(i)) * exp( - (abs(r_hat-theta(i)) / sqrt(tau_r+phi(i))).^2 );
    end

    lambda_tmp_2 = gamma * 1/(tau_r+psi) * exp( -(abs(r_hat) / sqrt(tau_r+psi)).^2);

    lambda_tmp_1 = sum(lambda_tmp_mat_1, 2);

    lambda_tmp_sum = lambda_tmp_1 + lambda_tmp_2;
    lambda_tmp_sum_zero_idx = (lambda_tmp_sum==0);    % some of the outliers 


    lambda_block_1 = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        lambda_block_1(:,i) = lambda_tmp_mat_1(:,i) ./ lambda_tmp_sum;
        lambda_block_1(lambda_tmp_sum_zero_idx,i) = 0;
    end

    lambda_block_2 = lambda_tmp_2 ./ lambda_tmp_sum;
    lambda_block_2(lambda_tmp_sum_zero_idx) = 1;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute omega, the weights of the Gaussian mixture components excluding the component for outliers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omega_new = sum(lambda_block_1);
    omega_new = omega_new.';    % this is transpose not conjungate transpose
    omega_new = omega_new / sum(omega_new);

    omega = omega + kappa * (omega_new - omega);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute gamma, the mixture weight of the outlier distribution
    % it is fixed, so we skip here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %gamma_new = sum(lambda_block_2) / (sum(sum(lambda_block_1)) + sum(lambda_block_2));
    %gamma = gamma + kappa * (gamma_new - gamma);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute theta, the mean
    % it is zero mean, so we skip here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    theta_tmp_mat = lambda_block_1;
    theta_tmp_mat_sum_1 = sum(theta_tmp_mat);
    zero_idx=(theta_tmp_mat_sum_1==0);  % find the guassian mixtures with zero weight

    theta = zeros(size(theta));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute phi, the variances of the Gaussian mixture components exclusing the component for outliers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % what happens if theta_tmp_mat is 0/0
    phi_tmp_mat_1 = theta_tmp_mat;
    phi_tmp_mat_2 = theta_tmp_mat;
    for (i=1:num_cluster)
        phi_tmp_mat_1(:,i) = phi_tmp_mat_1(:,i) .* (abs(r_hat-theta(i))).^2;
    end
    phi_new = sum(phi_tmp_mat_1) ./ sum(phi_tmp_mat_2) - tau_r;   % since the coefficients in tau_r are the same, just use tau_r(1)
    phi_new(zero_idx) = phi(zero_idx);  % keep the gaussian mixtures with zero weight fixed

    phi_new_inf_idx = isinf(phi_new);
    phi_new(phi_new_inf_idx) = phi(phi_new_inf_idx); % keep the gaussian mixtures with infinity phi fixed just in case

    phi_new = phi_new';

    for (i=1:num_cluster)
        if (phi_new(i)<0)
            phi_new(i) = phi(i);
        end
    end

    phi = phi + kappa * (phi_new - phi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute psi, the variance of the Gaussian component for outliers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    psi_tmp_mat = lambda_block_2;
    psi_tmp_mat_1 = psi_tmp_mat;
    psi_tmp_mat_2 = psi_tmp_mat;
    psi_new = sum(psi_tmp_mat_1 .* (abs(r_hat).^2)) / sum(psi_tmp_mat_2) - tau_r;    % since the coefficients in tau_r are the same, just use tau_r(1)

    if (psi_new<0)
        psi_new = psi;
    end

    psi = psi + kappa * (psi_new - psi);

end


