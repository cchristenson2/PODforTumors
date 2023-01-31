%{ 
Calibrate tumor for a specific patient using FOM in 2D

Inputs:
    - tumor struct
    - number of calibrated time points
    - number of predicted time points
    - parameter bounds
    - time spacing

Outputs:
    - params: calibrated outputs from LM algorithm
        - diffusivity, proliferation, treatment efficacy (A&C)
    - stats: calibrated statistics against measured data
        - CCC, DICE, Cell % rrror, Vol % error
    - outputs: simulation and prediction results

Contributors: Chase Christenson
%}

function [params, stats, outputs] = FOM_LMCalibration_GlobalKp_2D(tumor, ntp_cal, ntp_pred, bounds, dt)
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
    N0 = tumor.N(:,:,1);
    N_true = tumor.N(:,:,2:ntp_cal+1);
    if(ntp_pred==1)
        N_true_pred = tumor.N(:,:,end);
    end
    t = tumor.t_scan(2:1+ntp_cal);
    ROI = tumor.ROI;
    idx_ROI = find(ROI);
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

%     kp_g       = (exp(log(kp_low) + (log(kp_up)-log(kp_low)) * rand(1,1))) * ones(size(N0));
%     kp_g(~ROI) = 0;
%     d_g        = exp(log(d_low) + (log(d_up)-log(d_low)) * rand(1,1));
%     alpha1_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
%     alpha2_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
    
    kp_g     = kp_low * 5 *ones(size(N0));
    kp_g(~ROI) = 0;
    d_g      = d_low * 5;
    alpha1_g = alpha_low * 5;
    alpha2_g = alpha_low * 5;
    
    %Initialize SSE
    [N_sim, TC] = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_g, alpha2_g, tx_params, t, h, dt, bcs);
    SSE   = sum((N_true - N_sim).^2,'all');
    
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
        alpha2_g = alpha2_g/norm_alpha;
    end
    J_out = 0;
    res_out=0;
    %% Calibration loop
    %Calibration Loop
    iteration = 1;
    while(iteration<=max_it && SSE>e_tol)
        %Build Jacobian
        if(j_change==j_freq)
            %Perturb and find difference
            kp_p = kp_g*delta;
            dif_kp = kp_p(idx_ROI(1)) - kp_g(idx_ROI(1));
            d_p  = d_g*delta;
            dif_d  = d_p - d_g;
            alpha1_p = alpha1_g*delta;
            dif_alpha1 = alpha1_p - alpha1_g;
            alpha2_p = alpha2_g*delta;
            dif_alpha2 = alpha2_p - alpha2_g;
            
            if(norm == 1)
                kp_g = kp_g.*norm_kp;
                d_g  = d_g*norm_d;
                alpha1_g = alpha1_g*norm_alpha;
                alpha2_g = alpha2_g*norm_alpha;
                
                kp_p = kp_p.*norm_kp;
                d_p  = d_p*norm_d;
                alpha1_p = alpha1_p*norm_alpha;
                alpha2_p = alpha2_p*norm_alpha;
            end
            
            %Simulate perturbed params
            N_kpp    = RXDIF_2D_wAC(N0, kp_p, d_g, alpha1_g, alpha2_g, tx_params, t, h, dt, bcs);
            N_dp     = RXDIF_2D_wAC(N0, kp_g, d_p, alpha1_g, alpha2_g, tx_params, t, h, dt, bcs);
            N_alpha1 = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_p, alpha2_g, tx_params, t, h, dt, bcs);
            N_alpha2 = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_g, alpha2_p, tx_params, t, h, dt, bcs);
            
            %Fill J
            J = zeros(numel(N_true), 4);
            J(:,1) = reshape(N_kpp - N_sim, [], 1)/dif_kp;
            J(:,2) = reshape(N_dp - N_sim, [], 1)/dif_d;
            J(:,3) = reshape(N_alpha1 - N_sim, [], 1)/dif_alpha1;
            J(:,4) = reshape(N_alpha2 - N_sim, [], 1)/dif_alpha2;
            j_change = 0;
            
            if(norm == 1)
                kp_g = kp_g./norm_kp;
                d_g  = d_g/norm_d;
                alpha1_g = alpha1_g/norm_alpha;
                alpha2_g = alpha2_g/norm_alpha;
            end
        end
        
        
        
        %Calculate update with current regularization
        residuals = reshape(N_true - N_sim, [], 1);
%         [update,flags] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*residuals),1e-8,100);

        g = J' * residuals;
        H = (J'*J) + lambda*diag(diag(J'*J));

        gsc = zeros(size(H,1),1);
        Hsc = zeros(size(H,1),size(H,2));
        for i1 = 1:size(H,1)
            gsc(i1,1) = g(i1,1)/sqrt(H(i1,i1));
            for i2 = 1:size(H,2)
                Hsc(i1,i2) = H(i1,i2)/(sqrt(H(i1,i1))*sqrt(H(i2,i2))); 
            end
        end
        [L, U] = lu(Hsc);
        tempsc = L\gsc;
        delsc  = U\tempsc;
        update = zeros(size(delsc));
        for i1 = 1:size(delsc,1)
            update(i1,1) = delsc(i1,1)/sqrt(H(i1,i1));
        end
        
        %Update parameters to test
        kp_test = (kp_g(idx_ROI(1)) + update(1)) * ones(size(N0));
        kp_test(~ROI) = 0;
        if(kp_test(idx_ROI(1))<kp_low)
            kp_test(ROI) = kp_g(ROI) - (kp_g(ROI) - kp_low)./2;
        elseif(kp_test(idx_ROI(1))>kp_up)
            kp_test(ROI) = kp_g(ROI) + (kp_up - kp_G(ROI))./2;
        end
        d_test  = d_g + update(2);
        if(d_test<d_low)
            d_test = d_g - (d_g - d_low)/2;
        elseif(d_test>d_up)
            d_test = d_g + (d_up - d_g)/2;
        end
        alpha1_test  = alpha1_g + update(3);
        if(alpha1_test<alpha_low)
            alpha1_test = alpha1_g - (alpha1_g - alpha_low)/2;
        elseif(alpha1_test>alpha_up)
            alpha1_test = alpha1_g + (alpha_up - alpha1_g)/2;
        end
        alpha2_test  = alpha2_g + update(4);
        if(alpha2_test<alpha_low)
            alpha2_test = alpha2_g - (alpha2_g - alpha_low)/2;
        elseif(alpha2_test>alpha_up)
            alpha2_test = alpha2_g + (alpha_up - alpha2_g)/2;
        end
        
        if(norm == 1)
            kp_test = kp_test.*norm_kp;
            d_test  = d_test*norm_d;
            alpha1_test = alpha1_test*norm_alpha;
            alpha2_test = alpha2_test*norm_alpha;
        end
        
        %Test SSE calculation
        [N_test, TC_test] = RXDIF_2D_wAC(N0, kp_test, d_test, alpha1_test, alpha2_test, tx_params, t, h, dt, bcs);
        SSE_test = sum((N_true - N_test).^2,'all');
        
        if(norm == 1)
            kp_test = kp_test./norm_kp;
            d_test  = d_test/norm_d;
            alpha1_test = alpha1_test/norm_alpha;
            alpha2_test = alpha2_test/norm_alpha;
        end
        
        %Evaluate new error
        if(SSE_test<SSE) %Successful updates
            %Save new parameters
            N_sim = N_test;
            kp_g = kp_test;
            d_g  = d_test;
            alpha1_g = alpha1_test;
            alpha2_g = alpha2_test;
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
            
            J_out = J;
            res_out = residuals;
        else
            lambda = lambda*fail;
        end
        iteration = iteration + 1;
    end
    
    if(norm == 1)
        kp_g = kp_g.*norm_kp;
        d_g  = d_g*norm_d;
        alpha1_g = alpha1_g*norm_alpha;
        alpha2_g = alpha2_g*norm_alpha;
    end
    
    %% Make prediction (if needed)
    if(ntp_pred ~= 0)
        N_pred = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_g, alpha2_g, tx_params, tumor.t_scan(end), h, dt, bcs);
        outputs.Pred = N_pred;
    end
    
    
    
    %% Finalize results
    %Save calibrated parameter values
    params.kp = kp_g(idx_ROI(1));
    params.d  = d_g;
    params.alpha1 = alpha1_g;
    params.alpha2 = alpha2_g;
    
    %Calculate and save fit statistics (unfinished, add rest later)
    %CCC
    CCC_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        CCC_cal(1,i) = CCC_overlap(N_true(:,:,i), N_sim(:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        CCC_pred = CCC_overlap(N_true_pred, N_pred, thresh);
    else
        CCC_pred = [];
    end
    stats.CCC_cal  = CCC_cal;
    stats.CCC_pred = CCC_pred;
    
    %DICE
    DICE_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        DICE_cal(1,i) = DICE_calc(N_true(:,:,i), N_sim(:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        DICE_pred = DICE_calc(N_true_pred, N_pred, thresh);
    else
        DICE_pred = [];
    end
    stats.DICE_cal  = DICE_cal;
    stats.DICE_pred = DICE_pred;
    
    %Vol error
    Vol_err_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        Vol_err_cal(1,i) = getVolError(N_true(:,:,i), N_sim(:,:,i), tumor.h, tumor.dz, thresh);
    end
    if(ntp_pred ~= 0)
        Vol_err_pred = getVolError(N_true_pred, N_pred, tumor.h, tumor.dz, thresh);
    else
        Vol_err_pred = [];
    end
    stats.Vol_err_cal  = Vol_err_cal;
    stats.Vol_err_pred  = Vol_err_pred;
    
    %Cell error
    Cell_err_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        Cell_err_cal(1,i) = getCellError(N_true(:,:,i), N_sim(:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        Cell_err_pred = getCellError(N_true_pred, N_pred, thresh);
    else
        Cell_err_pred = [];
    end
    stats.Cell_err_cal  = Cell_err_cal;
    stats.Cell_err_pred  = Cell_err_pred;
    stats.confidence = ConfInterval(J_out, res_out);
    
    %Save final images
    outputs.Sim = N_sim;
    outputs.TC  = TC;
end