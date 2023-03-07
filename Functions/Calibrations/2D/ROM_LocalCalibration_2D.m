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
    - fig: calibration time course with results for each iteration

Contributors: Chase Christenson
%}

function [params, stats, outputs, fig] = ROM_LocalCalibration_2D(tumor, ntp_cal, ntp_pred, bounds, A_lib, T_lib, dt, model)
    %% Prep for calibration
    %LM parameters
    e_tol    = 1e-6;    %Calibration SSE goal
    e_conv   = 1e-10;    %Minimum change in SSE for an iteration
    max_it   = 1000;     %Maximum iterations
    delta    = 1.001;   %Perturbation magnitude
    delta_d  = 1.01;    %Delta perturb magnitude
    pass     = 7;       %Lambda reduction factor for successful updates
    fail     = 9;       %Lambda increase factor for unsuccessful updates
    lambda   = 1;       %Starting lambda
    j_freq   = 10;       %How many successful updates before updating J
    j_change = j_freq;       %Build J when equal to J frequency
    thresh   = 0.30;
    stuck_check = 0;
    successful_update = 0; guess_changed = 0;
    
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
    kp_up  = bounds.kp_bounds(end) / (delta^2);
    kp_low = bounds.kp_bounds(1) * (delta^2);
    d_up   = bounds.d_bounds(end) / (delta_d^2);
    d_low  = bounds.d_bounds(1) * (delta_d^2);
    alpha_up  = bounds.alpha_bounds(end) / (delta^2);
    alpha_low = bounds.alpha_bounds(1) * (delta^2);
    
%     kp_g       = kp_low * 5 * ones(size(N0));
%     kp_g(~ROI) = 0;

    d_g      = d_up/5;

    if(model~=0)
        alpha1_g = alpha_up/2;
        alpha2_g = alpha_up/2;
    end
    
    %% Prep for reduced order modeling
    start = tic;
    
    %Augment patient data and build kp
    if(model==0)
        [N_aug,kp_aug] = Augment_noTRX_2D(cat(3,N0,N_true), tumor.t_scan(2:end), h, dt, bcs, bounds, ntp_cal);
    elseif(model==1)
        [N_aug,kp_aug] = Augment_wTRX_comb_2D(cat(3,N0,N_true), tumor.t_scan(2:end), h, dt, bcs, bounds, ntp_cal, tx_params);
    else
        [N_aug,kp_aug] = Augment_wTRX_2D(cat(3,N0,N_true), tumor.t_scan(2:end), h, dt, bcs, bounds, ntp_cal, tx_params);
    end
    
    %Determine projection matrix
    [V,k] = getProjectionMatrix(N_aug, 0);
    
    %Reduce operator libraries
    Ar_lib = reduceALibrary(A_lib, V, bounds);
    if(model==0)
        clear T_lib;
    else
        Tr_lib = reduceTLibrary(T_lib, V, bounds, tx_params);
    end
    
    %Build libraries for local proliferation
    kpMode_bounds.low = zeros(1,k);
    kpMode_bounds.up  = zeros(1,k);
    for i = 1:k
        kp_test = zeros(numel(N0),1);
        kp_test(V(:,i)>=0) = kp_up;
        kp_test(V(:,i)<0) = kp_low;
        
        kpMode_bounds.up(i) = V(:,i)' * kp_test;

        kp_test = zeros(numel(N0),1);
        kp_test(V(:,i)>=0) = kp_low;
        kp_test(V(:,i)<0) = kp_up;
        
        kpMode_bounds.low(i) = V(:,i)' * kp_test;
    end
    [Br_lib, Hr_lib] = buildLocalProlifLibrary(N_true(:,:,1), V, k, 20, kpMode_bounds);
    
    t_extra_offline = toc(start);
    
    %Reduce MRI data
    N0_r = V'*N0(:);
    N_true_r = zeros(numel(N0_r), ntp_cal);
    for i = 1:ntp_cal
        N_true_r(:,i) = V'*reshape(N_true(:,:,i),[],1);
    end
    
    %Reduce kp_g and unreduce to start
    kp_gr = V' * kp_aug(:);
    kp_g  = V * kp_gr;
    kp_g(kp_g>kp_up)  = kp_up;
    kp_g(kp_g<kp_low) = kp_low;
    kp_gr = V' * kp_g(:);
    
    kp_gr0 = kp_gr;
    kp_g0 = kp_g;
    
    %Initialize Operators
    A_g = OperatorInterp_A(d_g,bounds,Ar_lib);
    if(model~=0)
        [Ta_g,Tc_g] = OperatorInterp_T(alpha1_g, alpha2_g, bounds, Tr_lib);
    else
        Ta_g = []; Tc_g = [];
    end
    
    B_g = OperatorInterp_local(kp_gr, Br_lib, k);
    H_g = OperatorInterp_local(kp_gr, Hr_lib, k);
    
    
%     B_g_test = zeros(k,k);
%     H_g_test = zeros(k,k^2);
%     for i = 1:numel(N0)
%         B_g_test = B_g_test + V(i,:)'*kp_g(i)*V(i,:);
%         H_g_test = H_g_test + V(i,:)'*kp_g(i)*kron(V(i,:),V(i,:));
%     end


    num_kp = k;
    %Initialize SSE
    if(model==0)
        [N_g_r, TC] = OperatorRXDIF_2D(N0_r, A_g, B_g, H_g, t, dt);
        num_p  = num_kp + 1;
    elseif(model==1)
        [N_g_r, TC] = OperatorRXDIF_2D_wAC_comb(N0_r, A_g, B_g, H_g, Ta_g, tx_params, t, dt);
        num_p  = num_kp + 2;
    else
        [N_g_r, TC] = OperatorRXDIF_2D_wAC(N0_r, A_g, B_g, H_g, Ta_g, Tc_g, tx_params, t, dt);
        num_p  = num_kp + 3;
    end
    SSE = sum((N_true_r - N_g_r).^2, 'all');
    
    if(isempty(gcp('nocreate')))
        myCluster = parcluster('Processes');
        parpool(myCluster.NumWorkers);
    end
    
    %Storage Variables
    kp_store = kp_g(:);
    kp_r_store = kp_gr(:);
    d_store = d_g;
    if(model~=0)
        alpha1_store = alpha1_g;
    end
    if(model==2)
        alpha2_store = alpha2_g;
    end
    SSE_store = SSE;
    lambda_store = lambda;
    
    
    %% Calibration loop
    iteration = 1;
    while(SSE>e_tol && iteration <= max_it)
        if(j_change == j_freq)
            %Build/Update jacobian
            J = zeros(numel(N0_r)*ntp_cal, num_p);

            %Perturb k in full space
            parfor(i = 1:num_kp)
                kp_pr = kp_gr;
                kp_pr(i) = kp_pr(i)*delta;
                dif_kp = kp_pr(i) - kp_gr(i);
                
                B_t = OperatorInterp_local(kp_pr, Br_lib, k);
                H_t = OperatorInterp_local(kp_pr, Hr_lib, k);
                
                if(model==0)
                    N_kp = OperatorRXDIF_2D(N0_r, A_g, B_t, H_t, t, dt);
                elseif(model==1)
                    N_kp = OperatorRXDIF_2D_wAC_comb(N0_r, A_g, B_t, H_t, Ta_g, tx_params, t, dt);
                else
                    N_kp = OperatorRXDIF_2D_wAC(N0_r, A_g, B_t, H_t, Ta_g, Tc_g, tx_params, t, dt);
                end
                
                J(:,i) = reshape(N_kp - N_g_r, [], 1)./dif_kp;
            end

            

            %Perturb d in reduced space
            d_p = d_g.*delta_d;
            dif_d = d_p - d_g;
            A_t = OperatorInterp_A(d_p, bounds, Ar_lib);
            
            if(model==0)
                N_d = OperatorRXDIF_2D(N0_r, A_t, B_g, H_g, t, dt);
                
                J(:,num_kp+1) = reshape(N_d - N_g_r, [], 1)./dif_d;
                
            elseif(model==1)
                N_d = OperatorRXDIF_2D_wAC_comb(N0_r, A_t, B_g, H_g, Ta_g, tx_params, t, dt);
                
                %Perturb alpha1 in reduced space
                alpha1_p = alpha1_g.*delta;
                dif_a1 = alpha1_p - alpha1_g;
                [Ta_t,~] = OperatorInterp_T(alpha1_p, alpha2_g, bounds, Tr_lib);
                N_alpha1 = OperatorRXDIF_2D_wAC_comb(N0_r, A_g, B_g, H_g, Ta_t, tx_params, t, dt);
                
                J(:,num_kp+1) = reshape(N_d - N_g_r, [], 1)./dif_d;
                J(:,num_kp+2) = reshape(N_alpha1 - N_g_r, [], 1)./dif_a1;

            else
                N_d = OperatorRXDIF_2D_wAC(N0_r, A_t, B_g, H_g, Ta_g, Tc_g, tx_params, t, dt);
                
                
                %Perturb alpha1 in reduced space
                alpha1_p = alpha1_g.*delta;
                dif_a1 = alpha1_p - alpha1_g;
                [Ta_t,Tc_t] = OperatorInterp_T(alpha1_p, alpha2_g, bounds, Tr_lib);
                N_alpha1 = OperatorRXDIF_2D_wAC(N0_r, A_g, B_g, H_g, Ta_t, Tc_t, tx_params, t, dt);
                
                
                %Perturb alpha2 in reduced space
                alpha2_p = alpha2_g.*delta;
                dif_a2 = alpha2_p - alpha2_g;
                [Ta_t,Tc_t] = OperatorInterp_T(alpha1_g, alpha2_p, bounds, Tr_lib);
                N_alpha2 = OperatorRXDIF_2D_wAC(N0_r, A_g, B_g, H_g, Ta_t, Tc_t, tx_params, t, dt);
                
                J(:,num_kp+1) = reshape(N_d - N_g_r, [], 1)./dif_d;
                J(:,num_kp+2) = reshape(N_alpha1 - N_g_r, [], 1)./dif_a1;
                J(:,num_kp+3) = reshape(N_alpha2 - N_g_r, [], 1)./dif_a2;
                
            end
            
            j_change = 0;
        end
        %Residuals & LM Calc
        residuals = reshape(N_true_r - N_g_r, [], 1);
        [update,flags] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*residuals),1e-8,100);
        
        %Update parameters
        kp_test_r = kp_gr + update(1:num_kp);
        kp_new = V*kp_test_r; kp_new = kp_new(idx_ROI);
        kp_new(kp_new<kp_low) = kp_g(idx_ROI(kp_new<kp_low)) - ((kp_g(idx_ROI(kp_new<kp_low)) - kp_low)./2);
        kp_new(kp_new>kp_up)  = kp_g(idx_ROI(kp_new>kp_up)) + ((kp_up - kp_g(idx_ROI(kp_new>kp_up)))./2);
        kp_test = zeros(size(kp_g)); kp_test(idx_ROI) = kp_new;
        kp_test_r = V' * kp_test;
        
        d_test = d_g + update(num_kp+1);
        if(d_test<d_low)
            d_test = d_g - (d_g - d_low)/2;
        elseif(d_test>d_up)
            d_test = d_g + (d_up - d_g)/2;
        end
        
        if(model~=0)
            alpha1_test = alpha1_g + update(num_kp+2);
            if(alpha1_test<alpha_low)
                alpha1_test = alpha1_g - (alpha1_g - alpha_low)/2;
            elseif(alpha1_test>alpha_up)
                alpha1_test = alpha1_g + (alpha_up - alpha1_g)/2;
            end
            alpha2_test = alpha2_g;
        end
        
        if(model==2)
            alpha2_test = alpha2_g + update(num_kp+3);
            if(alpha2_test<alpha_low)
                alpha2_test = alpha2_g - (alpha2_g - alpha_low)/2;
            elseif(alpha2_test>alpha_up)
                alpha2_test = alpha2_g + (alpha_up - alpha2_g)/2;
            end
        end

        %Run forward_check
        A_test = OperatorInterp_A(d_test, bounds, Ar_lib);
        if(model~=0)
            [Ta_test,Tc_test] = OperatorInterp_T(alpha1_test, alpha2_test, bounds, Tr_lib);
        else
            Ta_test = []; Tc_test = [];
        end
        B_test = OperatorInterp_local(kp_test_r, Br_lib, k);
        H_test = OperatorInterp_local(kp_test_r, Hr_lib, k);
        
        if(model==0)
            [N_test, TC_test] = OperatorRXDIF_2D(N0_r, A_test, B_test, H_test, t, dt);
        elseif(model==1)
            [N_test, TC_test] = OperatorRXDIF_2D_wAC_comb(N0_r, A_test, B_test, H_test, Ta_test, tx_params, t, dt);
        else
            [N_test, TC_test] = OperatorRXDIF_2D_wAC(N0_r, A_test, B_test, H_test, Ta_test, Tc_test, tx_params, t, dt);
        end
        SSE_test   = sum((N_true_r - N_test).^2,'all');

        if(SSE_test < SSE)
            N_g_r = N_test;
            
            kp_g  = kp_test;
            kp_gr = kp_test_r;
            d_g = d_test;
            if(model~=0)
            alpha1_g = alpha1_test;
            end
            if(model==2)
            alpha2_g = alpha2_test;
            end
            
            TC = TC_test;
            
            A_g = A_test; B_g = B_test; H_g = H_test; 
            if(model~=0)
                Ta_g = Ta_test; 
            end
            if(model==2)
                Tc_g = Tc_test;
            end

            %Check for convergence
            if(SSE-SSE_test < e_conv)
                disp(['ROM algorithm converged on iteration: ',num2str(iteration)]);
                break;
            end
            SSE = SSE_test;
            
            %Check for tolerance goal
            if(SSE<e_tol)
                disp(['ROM tolerance reached on iteration: ',num2str(iteration)]);
                break;
            end
            
            lambda = lambda/pass;
            
            j_change = j_change+1;
            stuck_check = 0;
            successful_update = successful_update + 1;
            
        else
            lambda = lambda*fail;
            if(lambda>1e20)
                lambda = 1e-20;
                j_change = j_freq;
                if(successful_update<=2 && guess_changed == 0) %Started in a local minimum, adjust initial guess
                    kp_g  = kp_up / 5 * ones(size(N0));
                    kp_g(~ROI) = 0;
                    kp_gr = V' * kp_g(:);
                    kp_g = V * kp_gr;
                    guess_changed = 1;
                elseif(stuck_check == 1) %searched whole range with no success
                    disp(['ROM algorithm stuck on iteration: ',num2str(iteration)]);
                    break;
                else %Cannot find an update, start over at minimum lambda
                    stuck_check = 1;
                end
            end
        end
        iteration = iteration+1; 
        
        kp_store = [kp_store, kp_g(:)];
        kp_r_store = [kp_r_store, kp_gr(:)];
        d_store   = [d_store, d_g];
        if(model~=0)
            alpha1_store = [alpha1_store, alpha1_g];
        end
        if(model==2)
            alpha2_store = [alpha2_store, alpha2_g];
        end
        SSE_store = [SSE_store, SSE];
        lambda_store = [lambda_store, lambda];
    end
    [sy,sx] = size(N0);
    kp_g = reshape(kp_g,sy,sx);
    
    %% Make prediction (if needed)
    if(ntp_pred ~= 0)
        %Option 1: predict in reduced state still
        A_r = OperatorInterp_A(d_g, bounds, Ar_lib);
        if(model~=0)
            [Ta_r,Tc_r] = OperatorInterp_T(alpha1_g, alpha2_g, bounds, Tr_lib);
        else
            Ta_r = []; Tc_r = [];
        end
        B_r = OperatorInterp_local(kp_gr, Br_lib, k);
        H_r = OperatorInterp_local(kp_gr, Hr_lib, k);

        
        if(model==0)
            N_pred_r = OperatorRXDIF_2D(N0_r, A_r, B_r, H_r, tumor.t_scan(end), dt);
        elseif(model==1)
            N_pred_r = OperatorRXDIF_2D_wAC_comb(N0_r, A_r, B_r, H_r, Ta_r, tx_params, tumor.t_scan(end), dt);
        else
            N_pred_r = OperatorRXDIF_2D_wAC(N0_r, A_r, B_r, H_r, Ta_r, Tc_r, tx_params, tumor.t_scan(end), dt);
        end
        N_pred_r_reshape = reshape(V*N_pred_r, size(N0));
        
        %Option 2: predict in full state with found parameters
        if(model==0)
            N_pred_full = RXDIF_2D(N0, kp_g, d_g, tumor.t_scan(end), h, dt, bcs);
        elseif(model==1)
            N_pred_full = RXDIF_2D_wAC_comb(N0, kp_g, d_g, alpha1_g, tx_params, tumor.t_scan(end), h, dt, bcs);
        else
            N_pred_full = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_g, alpha2_g, tx_params, tumor.t_scan(end), h, dt, bcs);
        end
        
        outputs.Pred = N_pred_r_reshape;
    end
    %% Finalize results
    [sy,sx] = size(N0);
    %Save calibrated parameter values
    params.kp = reshape(kp_g,sy,sx);
    params.kp_r = kp_gr;
    params.d  = d_g;
    if(model~=0)
        params.alpha1 = alpha1_g;
    end
    if(model==2)
        params.alpha2 = alpha2_g;
    end
    
    %Resize and shape final simulations from calibration
    N_g_final = zeros(numel(N0),ntp_cal);
    N_sim     = zeros([sy,sx,ntp_cal]);
    for i = 1:ntp_cal
        N_g_final(:,i) = V*N_g_r(:,i);
        N_sim(:,:,i) = reshape(N_g_final(:,i), sy,sx);
    end

    %Calculate and save fit statistics (unfinished, add rest later)
    CCC_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        CCC_cal(1,i) = CCC_overlap(N_true(:,:,i), N_sim(:,:,i), thresh);
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
        DICE_cal(1,i) = DICE_calc(N_true(:,:,i), N_sim(:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        DICE_pred(1) = DICE_calc(N_true_pred, N_pred_r_reshape, thresh);
        DICE_pred(2) = DICE_calc(N_true_pred, N_pred_full, thresh);
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
        Cell_err_cal(1,i) = getCellError(N_true(:,:,i), N_sim(:,:,i), thresh);
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
    outputs.iteration = iteration;
    outputs.V   = V;
    outputs.t_extra = t_extra_offline;
    
    outputs.testing.kp_gr0 = kp_gr0;
    outputs.testing.kp_g0  = kp_g0;
    outputs.testing.updates = successful_update;
    
    %Plotting figure
    fig = figure;
    subplot(2,4,1)
    plot(1:iteration,kp_store);
    xlabel('Iteration'); ylabel('Proliferation Rate (1/day)');
    
    subplot(2,4,2)
    plot(1:iteration, d_store);
    xlabel('Iteration'); ylabel('Diffusivity (mm^2/day)');
    
    if(model~=0)
        subplot(2,4,3)
        plot(1:iteration, alpha1_store);
        xlabel('Iteration'); ylabel('Alpha 1 (1/day)');
    end
    
    if(model==2)
        subplot(2,4,4)
        plot(1:iteration, alpha2_store);
        xlabel('Iteration'); ylabel('Alpha 2 (1/day)');
    end
    
    subplot(2,4,5)
    semilogy(1:iteration, SSE_store);
    xlabel('Iteration'); ylabel('SSE');
    
    subplot(2,4,6)
    semilogy(1:iteration, lambda_store);
    xlabel('Iteration'); ylabel('Lambda');
    
    subplot(2,4,7)
    plot(1:iteration, kp_r_store);
    xlabel('Iteration'); ylabel('Proliferation Rate Reduced(1/day)');
    
end