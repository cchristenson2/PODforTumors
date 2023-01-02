%{ 
Calibrate tumor for a specific patient using ROM in 2D

Inputs:
    - tumor struct
    - number of calibrated time points
    - number of predicted time points
    - parameter bounds
    - Operator libraries
    - time spacing

Outputs:
    - params: calibrated outputs from LM algorithm
        - diffusivity, proliferation, treatment efficacy (A&C)
    - stats: calibrated statistics against measured data
        - CCC, DICE, Cell % rrror, Vol % error
    - outputs: simulation and prediction results

Contributors: Chase Christenson
%}

function [params, stats, outputs] = ROM_LMCalibration_GlobalKp_comb_3D(tumor, ntp_cal, ntp_pred, bounds, k, A_lib, B_lib, H_lib, T_lib, dt)
%     addpath(genpath('C:/ROM/'))
    
    %% Prep for calibration
    %LM parameters
    e_tol    = 1e-5;    %Calibration SSE goal
    e_conv   = 1e-7;    %Minimum change in SSE for an iteration
    max_it   = 100;     %Maximum iterations
    delta    = 1.001;   %Perturbation magnitude
    pass     = 7;       %Lambda reduction factor for successful updates
    fail     = 9;       %Lambda increase factor for unsuccessful updates
    lambda   = 1;       %Starting lambda
    j_change = 1;       %Build J when equal to J frequency
    j_freq   = 1;       %How many successful updates before updating J
    thresh   = 0.15;
    
    norm = 1;
    
    %Pull out struct variables
    N0 = tumor.N(:,:,:,1);
    N_true = tumor.N(:,:,:,2:ntp_cal+1);
    if(ntp_pred==1)
        N_true_pred = tumor.N(:,:,:,end);
    end
    t = tumor.t_scan(2:1+ntp_cal);
%     ROI = tumor.ROI;
%     idx_ROI = find(ROI);
    bcs = tumor.bcs;
    h = tumor.h;
    
    tx_params.C = tumor.AUC;
    tx_params.beta1 = tumor.beta1;
    tx_params.beta2 = tumor.beta2;
    tx_params.txduration = tumor.t_trx;
    

    %Pull out parameter bounds and set initial guesses
    kp_up  = bounds.kp_bounds(end);
    kp_low = bounds.kp_bounds(1);
    d_up   = bounds.d_bounds(end);
    d_low  = bounds.d_bounds(1);
    alpha_up  = bounds.alpha_bounds(end);
    alpha_low = bounds.alpha_bounds(1);

%     kp_g       = (exp(log(kp_low) + (log(kp_up)-log(kp_low)) * rand(1,1)));
%     kp_g(~ROI) = 0;
%     d_g        = exp(log(d_low) + (log(d_up)-log(d_low)) * rand(1,1));
%     alpha1_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
%     alpha2_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
    
    kp_g     = kp_low * 5;
    d_g      = d_low * 5;
    alpha1_g = alpha_low * 5;
%     alpha2_g = alpha_low * 5;
    
    %% Prep for reduced order modeling
    %Determine projection matrix
    V = getProjectionMatrix(cat(4,N0,N_true), k);
    
    %Reduce operator libraries
    [Ar_lib,Br_lib,Hr_lib,Tr_lib] = reduceOperatorLibrary(A_lib, B_lib, H_lib, T_lib, V, bounds, tx_params);
    
    clear A_lib B_lib H_lib T_lib;
    
    %Reduce MRI data
    N0_r = V'*N0(:);
    N_true_r = zeros(numel(N0_r), ntp_cal);
    for i = 1:ntp_cal
        N_true_r(:,i) = V'*reshape(N_true(:,:,:,i),[],1);
    end
    
    %Initialize Operators
    [A_g,B_g,H_g,Ta_g] = OperatorInterp_wAC_comb(kp_g,d_g,alpha1_g,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
    
    %Initialize SSE
    [N_g_r, TC] = OperatorRXDIF_2D_wAC_comb(N0_r, A_g, B_g, H_g, Ta_g, tx_params, t, dt);
    SSE = sum((N_true_r - N_g_r).^2, 'all');
    
    if(norm == 1)
        norm_kp = kp_up;
        norm_d = d_up;
        norm_alpha = alpha_up;
        
        kp_up = 1;
        d_up = 1;
        alpha_up = 1;
        
        kp_low = kp_low/norm_kp;
        d_low = d_low/norm_d;
        alpha_low = alpha_low/norm_alpha;
        
        kp_g = kp_g./norm_kp;
        d_g  = d_g/norm_d;
        alpha1_g = alpha1_g/norm_alpha;
%         alpha2_g = alpha2_g/norm_alpha;
    end
    
    %% Calibration loop
    iteration = 1;
    while(SSE>e_tol && iteration <= max_it)
        if(j_change == j_freq)
            %Build/Update jacobian
            J = zeros(numel(N0_r)*ntp_cal, 3);
            
            %Perturb parameters
            kp_p = kp_g.*delta;
            dif_kp = kp_p - kp_g;
            d_p = d_g.*delta;
            dif_d = d_p - d_g;
            alpha1_p = alpha1_g.*delta;
            dif_alpha1 = alpha1_p - alpha1_g;
%             alpha2_p = alpha2_g.*delta;
%             dif_alpha2 = alpha2_p - alpha2_g;
            
            if(norm == 1)
                kp_g = kp_g.*norm_kp;
                d_g  = d_g*norm_d;
                alpha1_g = alpha1_g*norm_alpha;
%                 alpha2_g = alpha2_g*norm_alpha;
                
                kp_p = kp_p.*norm_kp;
                d_p  = d_p*norm_d;
                alpha1_p = alpha1_p*norm_alpha;
%                 alpha2_p = alpha2_p*norm_alpha;
            end
            
            %Perturb k
            [A_t, B_t, H_t, Ta_t] = OperatorInterp_wAC_comb(kp_p,d_g,alpha1_g,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
            N_kpp = OperatorRXDIF_2D_wAC_comb(N0_r, A_t, B_t, H_t, Ta_t, tx_params, t, dt);
            J(:,1) = reshape(N_kpp - N_g_r, [], 1)./dif_kp;

            %Perturb d
            [A_t, B_t, H_t, Ta_t] = OperatorInterp_wAC_comb(kp_g,d_p,alpha1_g,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
            N_dp = OperatorRXDIF_2D_wAC_comb(N0_r, A_t, B_t, H_t, Ta_t, tx_params, t, dt);
            J(:,2) = reshape(N_dp - N_g_r, [], 1)./dif_d;

            %Perturb alpha1
            [A_t, B_t, H_t, Ta_t] = OperatorInterp_wAC_comb(kp_g,d_g,alpha1_p,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
            N_alpha1 = OperatorRXDIF_2D_wAC_comb(N0_r, A_t, B_t, H_t, Ta_t, tx_params, t, dt);
            J(:,3) = reshape(N_alpha1 - N_g_r, [], 1)./dif_alpha1;

            %Perturb alpha2
%             [A_t, B_t, H_t, Ta_t, Tc_t] = OperatorInterp_wAC_comb(kp_g,d_g,alpha1_g,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
%             N_alpha2 = OperatorRXDIF_2D_wAC_comb(N0_r, A_t, B_t, H_t, Ta_t, tx_params, t, dt);
%             J(:,4) = reshape(N_alpha2 - N_g_r, [], 1)./dif_alpha2;
            
            if(norm == 1)
                kp_g = kp_g./norm_kp;
                d_g  = d_g/norm_d;
                alpha1_g = alpha1_g/norm_alpha;
%                 alpha2_g = alpha2_g/norm_alpha;
            end
            
            j_change = 0;
        end
        %Residuals & LM Calc
        residuals = reshape(N_true_r - N_g_r, [], 1);
        [update,flags] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*residuals),1e-8,100);

        %Update parameters
        kp_test = kp_g + update(1);
        if(kp_test<kp_low)
            kp_test = kp_low;
        elseif(kp_test>kp_up)
            kp_test = kp_up;
        end
        d_test = d_g + update(2);
        if(d_test<d_low)
            d_test = d_low;
        elseif(d_test>d_up)
            d_test = d_up;
        end
        alpha1_test = alpha1_g + update(3);
        if(alpha1_test<alpha_low)
            alpha1_test = alpha_low;
        elseif(alpha1_test>alpha_up)
            alpha1_test = alpha_up;
        end
%         alpha2_test = alpha2_g + update(4);
%         if(alpha2_test<alpha_low)
%             alpha2_test = alpha_low;
%         elseif(alpha2_test>alpha_up)
%             alpha2_test = alpha_up;
%         end
        
        if(norm == 1)
            kp_test = kp_test.*norm_kp;
            d_test  = d_test*norm_d;
            alpha1_test = alpha1_test*norm_alpha;
%             alpha2_test = alpha2_test*norm_alpha;
        end

        %Run forward_check
        [A_test, B_test, H_test, Ta_test] = OperatorInterp_wAC_comb(kp_test,d_test,alpha1_test,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
        [N_test, TC_test] = OperatorRXDIF_2D_wAC_comb(N0_r, A_test, B_test, H_test, Ta_test, tx_params, t, dt);       
        SSE_test   = sum((N_true_r - N_test).^2,'all');

        if(norm == 1)
            kp_test = kp_test./norm_kp;
            d_test  = d_test/norm_d;
            alpha1_test = alpha1_test/norm_alpha;
%             alpha2_test = alpha2_test/norm_alpha;
        end
        
        if(SSE_test < SSE)
            N_g_r = N_test;
            kp_g = kp_test;
            d_g = d_test;
            alpha1_g = alpha1_test;
%             alpha2_g = alpha2_test;
            TC = TC_test;

            %Check for convergence
            if(SSE-SSE_test < e_conv)
                disp(['Algorithm converged on iteration: ',num2str(iteration)]);
                break;
            end
            SSE = SSE_test;
            
            %Check for tolerance goal
            if(SSE<e_tol)
                disp(['Tolerance reached on iteration: ',num2str(iteration)]);
                break;
            end
            lambda = lambda/pass;
            j_change = j_change+1;
        else
            lambda = lambda*fail;
        end
        iteration = iteration+1; 
        
        
    end
    
    if(norm == 1)
        kp_g = kp_g.*norm_kp;
        d_g  = d_g*norm_d;
        alpha1_g = alpha1_g*norm_alpha;
%         alpha2_g = alpha2_g*norm_alpha;
    end

    %% Make prediction (if needed)
    if(ntp_pred ~= 0)
        %Option 1: predict in reduced state still
        [A_r, B_r, H_r, Ta_r] = OperatorInterp_wAC_comb(kp_g,d_g,alpha1_g,bounds,Ar_lib,Br_lib,Hr_lib,Tr_lib);
        N_pred_r = OperatorRXDIF_2D_wAC_comb(N0_r, A_r, B_r, H_r, Ta_r, tx_params, tumor.t_scan(end), dt);
        N_pred_r_reshape = reshape(V*N_pred_r, size(N0));
        
        %Option 2: predict in full state with found parameters
        N_pred_full = RXDIF_3D_wAC_comb(N0, kp_g*ones(size(N0)), d_g, alpha1_g, tx_params, tumor.t_scan(end), h, dt, bcs);
        outputs.Pred = N_pred_r_reshape;
    end
    %% Finalize results
    %Save calibrated parameter values
    params.kp = kp_g;
    params.d  = d_g;
    params.alpha1 = alpha1_g;
%     params.alpha2 = alpha2_g;
    
    %Resize and shape final simulations from calibration
    [sy,sx,sz] = size(N0);
    N_g_final = zeros(numel(N0),ntp_cal);
    N_sim     = zeros([sy,sx,sz,ntp_cal]);
    for i = 1:ntp_cal
        N_g_final(:,i) = V*N_g_r(:,i);
        N_sim(:,:,:,i) = reshape(N_g_final(:,i), sy,sx,sz);
    end

    %Calculate and save fit statistics (unfinished, add rest later)
    CCC_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        CCC_cal(1,i) = CCC_overlap(N_true(:,:,:,i), N_sim(:,:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        CCC_pred(1) = CCC_overlap(N_true_pred, N_pred_r_reshape, thresh);
        CCC_pred(2) = CCC_overlap(N_true_pred, N_pred_full, thresh);
    else
        CCC_pred = [];
    end
    stats.CCC_cal  = CCC_cal;
    stats.CCC_pred = CCC_pred;
    
    %DICE
    DICE_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        DICE_cal(1,i) = DICE_calc(N_true(:,:,:,i), N_sim(:,:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        DICE_pred(1) = DICE_calc(N_true_pred, N_pred_r_reshape, thresh);
        DICE_pred(1) = DICE_calc(N_true_pred, N_pred_full, thresh);
    else
        DICE_pred = [];
    end
    stats.DICE_cal  = DICE_cal;
    stats.DICE_pred = DICE_pred;
    
    %Vol error
    Vol_err_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        Vol_err_cal(1,i) = getVolError(N_true(:,:,:,i), N_sim(:,:,:,i), tumor.h, tumor.dz, thresh);
    end
    if(ntp_pred ~= 0)
        Vol_err_pred(1) = getVolError(N_true_pred, N_pred_r_reshape, tumor.h, tumor.dz, thresh);
        Vol_err_pred(2) = getVolError(N_true_pred, N_pred_full, tumor.h, tumor.dz, thresh);
    else
        Vol_err_pred = [];
    end
    stats.Vol_err_cal  = Vol_err_cal;
    stats.Vol_err_pred  = Vol_err_pred;
    
    %Cell error
    Cell_err_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        Cell_err_cal(1,i) = getCellError(N_true(:,:,:,i), N_sim(:,:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        Cell_err_pred(1) = getCellError(N_true_pred, N_pred_r_reshape, thresh);
        Cell_err_pred(2) = getCellError(N_true_pred, N_pred_full, thresh);
    else
        Cell_err_pred = [];
    end
    stats.Cell_err_cal  = Cell_err_cal;
    stats.Cell_err_pred  = Cell_err_pred;
    
    %Save final images
    outputs.Sim = N_sim;
    outputs.TC  = TC;
    
    
end