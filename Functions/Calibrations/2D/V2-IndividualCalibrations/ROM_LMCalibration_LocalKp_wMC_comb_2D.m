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

function [params, stats, outputs, fig, temp] = ROM_LMCalibration_LocalKp_wMC_comb_2D(tumor, ntp_cal, ntp_pred, bounds, k, A_lib, B_lib, H_lib, T_lib, dt)
%     addpath(genpath('C:/ROM/'))
    
    %% Prep for calibration
    %LM parameters
    e_tol    = 1e-20;    %Calibration SSE goal
    e_conv   = 1e-7;    %Minimum change in SSE for an iteration
    max_it   = 500;     %Maximum iterations
    delta    = 1.001;   %Perturbation magnitude
    delta_kp = 1.001;
    pass     = 7;       %Lambda reduction factor for successful updates
    fail     = 9;       %Lambda increase factor for unsuccessful updates
    lambda   = 1;       %Starting lambda
    j_change = 3;       %Build J when equal to J frequency
    j_freq   = 3;       %How many successful updates before updating J
    thresh   = 0.15;
    
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
    tissues = tumor.Tissues;
    
    tx_params.C = tumor.AUC;
    tx_params.beta1 = tumor.beta1;
    tx_params.beta2 = tumor.beta2;
    tx_params.txduration = tumor.t_trx;
    

    %Pull out parameter bounds and set initial guesses
    kp_up  = bounds.kp_bounds(end) / (delta^2);
    kp_low = bounds.kp_bounds(1) * (delta^2);
    d_up   = bounds.d_bounds(end) / (delta^2);
    d_low  = bounds.d_bounds(1) * (delta^2);
    alpha_up  = bounds.alpha_bounds(end) / (delta^2);
    alpha_low = bounds.alpha_bounds(1) * (delta^2);
    
    %Set up mechanics variables
    [M, E, nu] =  mech_matrix_build_2D_v5(h, tissues, bcs);
    M = sparse(M);
    % Assemble matrices that calc gradient
    [sx,sy] = size(tissues);
    [matX, matY] = gradN_matrix(sx, sy, h, bcs);
    
    %Augment patient data and build kp
    [N_aug,kp_aug] = Augment_wMC_comb_2D(cat(3,N0,N_true), tumor.t_scan(2:end), h, 1, bcs, bounds, ntp_cal, tx_params, M, E, nu, matX, matY);

%     kp_g       = (exp(log(kp_low) + (log(kp_up)-log(kp_low)) * rand(1,1))) * ones(size(N0));
%     kp_g       = kp_low * 5 * ones(size(N0));
%     kp_g(~ROI) = 0;
%     d_g        = exp(log(d_low) + (log(d_up)-log(d_low)) * rand(1,1));
%     alpha1_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
%     alpha2_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
    d_g      = d_up/2;
%     alpha1_g = alpha_low * 5;
%     alpha2_g = alpha_low * 5;
    alpha1_g = 0.3; alpha2_g = 0.3;
    
    %Get stress maps for augmented data
    Stress_aug = getStressMaps_2D(N_aug, M, E, nu, matX, matY);
    
    
    
    %% Prep for reduced order modeling
    %Determine projection matrix
%     V = getProjectionMatrix(cat(3,N0,N_true), k);
    [V,k,V_full] = getProjectionMatrix(N_aug,0);
    
    [V_s,~,~] = getProjectionMatrix(Stress_aug,k);
    
    %Reduce mechanics variables
    M_r = reduceM(M,V);
    matX_r = V_s' * matX * V_s;
    matY_r = V_s' * matY * V_s;
    E_r = V_s' * E(:);
    
    %Reduce operator libraries
    [Ar_lib,Br_lib,Hr_lib,Tr_lib] = reduceOperatorLibrary(A_lib, B_lib, H_lib, T_lib, V, bounds, tx_params);
    
    %Reduce MRI data
    N0_r = V'*N0(:);
    N_true_r = zeros(numel(N0_r), ntp_cal);
    for i = 1:ntp_cal
        N_true_r(:,i) = V'*reshape(N_true(:,:,i),[],1);
    end
    
        %Reduce kp_g and unreduce to start
%     kp_gr = V' * kp_g(:);
%     kp_g  = V * kp_gr;

%     kp_gr = 0.01.*ones(k,1);
%     kp_g = V * kp_gr;
    kp_g = kp_aug(:);
    kp_gr = V' * kp_g;
    kp_g0 = kp_g;
    

    
    %Initialize Operators
    A_g = OperatorInterp_A(d_g, bounds, Ar_lib);
    [Ta_g,~] = OperatorInterp_T(alpha1_g, alpha2_g, bounds, Tr_lib);
    
    B_f = assembleB(numel(N0), kp_g(:));
    B_g = V'*B_f*V;
    
    H_f = assembleH(numel(N0), kp_g(:));
    H_g = V'*H_f*kron(V,V);
    
%     kron_proj = kron(V,V);
    
    %Initialize SSE
    [N_g_r, TC] = OperatorRXDIF_2D_wMC_comb(N0_r, A_g, d_g, B_g, H_g, Ta_g, tx_params, t, dt, M_r, E_r, nu, matX_r, matY_r);
    
    
    SSE = sum((V*N_true_r - V*N_g_r).^2, 'all');
    
    if(isempty(gcp('nocreate')))
        myCluster = parcluster('Processes');
        parpool(myCluster.NumWorkers);
    end
    
    %Storage Variables
    kp_store = kp_g(:);
    kp_r_store = kp_gr(:);
    d_store = d_g;
    alpha1_store = alpha1_g;
%     alpha2_store = alpha2_g;
    SSE_store = SSE;
    lambda_store = lambda;
    
    num_kp = k;
    num_p  = num_kp + 2;
    
    %% Calibration loop
    iteration = 1;
    while(SSE>e_tol && iteration <= max_it)
        if(j_change == j_freq)
            %Build/Update jacobian
            J = zeros(numel(N0_r)*ntp_cal, num_p);

            %Perturb k in full space
%             start = tic;
            parfor(i = 1:num_kp)
                kp_pr = kp_gr;
                kp_pr(i) = kp_pr(i)*delta_kp;
                dif_kp = kp_pr(i) - kp_gr(i);
                
                kp_p = V * kp_pr;
                difference = kp_p - kp_g;
                idx = find(difference);
                
                B_change = 0; H_change = 0;
                for j = 1:numel(idx)
                    ind = idx(j);
                    V_row = V(ind,:);
                    B_change = B_change + V_row' * difference(ind) * V_row;
                    H_change = H_change + V_row' * difference(ind) * kron(V_row,V_row);
                end
                
                B_t = B_g + B_change; H_t = H_g + H_change;
                
                N_kp = OperatorRXDIF_2D_wMC_comb(N0_r, A_g, d_g, B_t, H_t, Ta_g, tx_params, t, dt, M_r, E_r, nu, matX_r, matY_r);
                
                J(:,i) = reshape(N_kp - N_g_r, [], 1)./dif_kp;
            end
%             stop = toc(start);
%             if(iteration==1)
%                 disp(['kp full perturb ROM = ',num2str(stop),' sec']);
%             end

            %Perturb d in reduced space
            d_p = d_g.*delta;
            dif = d_p - d_g;
            A_t = OperatorInterp_A(d_p, bounds, Ar_lib);
            N_d = OperatorRXDIF_2D_wMC_comb(N0_r, A_t, d_p, B_g, H_g, Ta_g, tx_params, t, dt, M_r, E_r, nu, matX_r, matY_r);
            J(:,num_kp+1) = reshape(N_d - N_g_r, [], 1)./dif;

            %Perturb alpha1 in reduced space
            alpha1_p = alpha1_g.*delta;
            dif = alpha1_p - alpha1_g;
            [Ta_t,~] = OperatorInterp_T(alpha1_p, alpha2_g, bounds, Tr_lib);
            N_alpha1 = OperatorRXDIF_2D_wMC_comb(N0_r, A_g, d_g, B_g, H_g, Ta_t, tx_params, t, dt, M_r, E_r, nu, matX_r, matY_r);
            J(:,num_kp+2) = reshape(N_alpha1 - N_g_r, [], 1)./dif;

            %Perturb alpha2 in reduced space
%             alpha2_p = alpha2_g.*delta;
%             dif = alpha2_p - alpha2_g;
%             [Ta_t,Tc_t] = OperatorInterp_T(alpha1_g, alpha2_p, bounds, Tr_lib);
%             N_alpha2 = OperatorRXDIF_2D_wMC_wAC(N0_r, A_g, d_g, B_g, H_g, Ta_t, Tc_t, tx_params, t, dt, M_r, E_r, nu, matX_r, matY_r);
%             J(:,num_kp+3) = reshape(N_alpha2 - N_g_r, [], 1)./dif;

            j_change = 0;
        end
        %Residuals & LM Calc
        residuals = reshape(N_true_r - N_g_r, [], 1);
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
        if(numel(find(isnan(update)))~=0 || numel(find(update))==0)
            [update,flags] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*residuals),1e-8,100);
        end
        
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
        alpha1_test = alpha1_g + update(num_kp+2);
        if(alpha1_test<alpha_low)
            alpha1_test = alpha1_g - (alpha1_g - alpha_low)/2;
        elseif(alpha1_test>alpha_up)
            alpha1_test = alpha1_g + (alpha_up - alpha1_g)/2;
        end
%         alpha2_test = alpha2_g + update(num_kp+3);
%         if(alpha2_test<alpha_low)
%             alpha2_test = alpha2_g - (alpha2_g - alpha_low)/2;
%         elseif(alpha2_test>alpha_up)
%             alpha2_test = alpha2_g + (alpha_up - alpha2_g)/2;
%         end
        alpha2_test = alpha2_g;

        %Run forward_check
        A_test = OperatorInterp_A(d_test, bounds, Ar_lib);
        [Ta_test,~] = OperatorInterp_T(alpha1_test, alpha2_test, bounds, Tr_lib);

%         B_f = assembleB(numel(N0), kp_test(:));
%         B_test = V'*B_f*V;
%         H_f = assembleH(numel(N0), kp_test(:));
%         H_test = V'*H_f*kron(V,V);
        
        dif_kp = kp_test(:) - kp_g(:);
        B_change = zeros(size(B_g));
        H_change = zeros(size(H_g));
        idx = find(dif_kp);
        for i = 1:numel(idx)
            ind = idx(i);
            B_change = B_change + V(ind,:)' * dif_kp(ind) * V(ind,:);
            H_change = H_change + V(ind,:)' * dif_kp(ind) * kron(V(ind,:),V(ind,:));
        end
        B_test = B_g + B_change;
        H_test = H_g + H_change;
        
        [N_test, TC_test] = OperatorRXDIF_2D_wMC_comb(N0_r, A_test, d_test, B_test, H_test, Ta_test, tx_params, t, dt, M_r, E_r, nu, matX_r, matY_r);       
        SSE_test   = sum((V*N_true_r - V*N_test).^2,'all');

        if(SSE_test < SSE)
            N_g_r = N_test;
            kp_g  = kp_test;
            kp_gr = kp_test_r;
            d_g = d_test;
            alpha1_g = alpha1_test;
%             alpha2_g = alpha2_test;
            TC = TC_test;
            
            A_g = A_test; B_g = B_test; H_g = H_test; Ta_g = Ta_test; %Tc_g = Tc_test;

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
            
            %Temp variables for testing validation
%             [update,flags] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*residuals),1e-8,100);
            
            temp.updates = update;
            temp.lambda = lambda;
            temp.residuals = residuals;
            temp.N_r = N_true_r;
            temp.N_g_r = N_g_r;
            temp.J = J;
            
        else
            lambda = lambda*fail;
            if(lambda>1e30)
                lambda = 1e-30;
            end
        end
        iteration = iteration+1; 
        
        kp_store = [kp_store, kp_g(:)];
        kp_r_store = [kp_r_store, kp_gr(:)];
        d_store   = [d_store, d_g];
        alpha1_store = [alpha1_store, alpha1_g];
%         alpha2_store = [alpha2_store, alpha2_g];
        SSE_store = [SSE_store, SSE];
        lambda_store = [lambda_store, lambda];
    end
    [sy,sx] = size(N0);
    kp_g = reshape(kp_g,sy,sx);
    
    %% Make prediction (if needed)
    if(ntp_pred ~= 0)
        %Option 1: predict in reduced state still
        A_r = OperatorInterp_A(d_g, bounds, Ar_lib);
        [Ta_r,~] = OperatorInterp_T(alpha1_g, alpha2_g, bounds, Tr_lib);

        B_fr = assembleB(numel(N0), kp_g(:));
        B_r = V'*B_fr*V;

        H_fr = assembleH(numel(N0), kp_g(:));
        H_r = V'*H_fr*kron(V,V);
        
        N_pred_r = OperatorRXDIF_2D_wMC_comb(N0_r, A_r, d_g, B_r, H_r, Ta_r, tx_params, tumor.t_scan(end), dt, M_r, E_r, nu, matX_r, matY_r);
        N_pred_r_reshape = reshape(V*N_pred_r, size(N0));
        
        %Option 2: predict in full state with found parameters
        N_pred_full = RXDIF_2D_wMC_comb(N0, kp_g, d_g, alpha1_g, tx_params, tumor.t_scan(end), h, dt, bcs, M, E, nu, matX, matY);
        
        outputs.Pred = N_pred_r_reshape;
    end
    %% Finalize results
    [sy,sx] = size(N0);
    %Save calibrated parameter values
    params.kp = reshape(kp_g,sy,sx);
    params.kp_r = kp_gr;
    params.d  = d_g;
    params.alpha1 = alpha1_g;
%     params.alpha2 = alpha2_g;
    
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
    outputs.V_full = V_full;
%     outputs.Br  = B_g;
%     outputs.Hr  = H_g;
    
    %Plotting figure
    fig = figure;
    subplot(2,4,1)
    plot(1:iteration,kp_store);
    xlabel('Iteration'); ylabel('Proliferation Rate (1/day)');
    subplot(2,4,2)
    plot(1:iteration, kp_r_store);
    xlabel('Iteration'); ylabel('Proliferation Rate Reduced(1/day)');
    subplot(2,4,3)
    plot(1:iteration, alpha1_store);
    xlabel('Iteration'); ylabel('Alpha 1 (1/day)');
%     subplot(2,4,4)
%     plot(1:iteration, alpha2_store);
%     xlabel('Iteration'); ylabel('Alpha 2 (1/day)');
    subplot(2,4,5)
    semilogy(1:iteration, SSE_store);
    xlabel('Iteration'); ylabel('SSE');
    subplot(2,4,6)
    plot(1:iteration, lambda_store);
    xlabel('Iteration'); ylabel('Lambda');
    subplot(2,4,7)
    plot(1:iteration, d_store);
    xlabel('Iteration'); ylabel('Diffusivity (mm^2/day)');
%     savefig(f1, [home,'/Results/',erase(tumor.name,'.mat'),'_ROM_Calibration_local2D']);
end