%{ 
Calibrate tumor for a specific patient using FOM in 2D

Inputs:
    - tumor struct
    - number of calibrated time points
    - number of predicted time points
    - parameter bounds
    - ROM libraries
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

function [params, stats, outputs, fig] = FOM_LMCalibration_LocalKp_noTRX_3D(tumor, ntp_cal, ntp_pred, bounds, dt, workers)
%     addpath(genpath('C:/ROM/'))
    
    %% Prep for calibration
    %LM parameters
    e_tol    = 1e-6;    %Calibration SSE goal
    e_conv   = 1e-10;    %Minimum change in SSE for an iteration
    max_it   = 1000;     %Maximum iterations
    delta    = 1.001;   %Perturbation magnitude
    delta_d  = 1.01;
    pass     = 7;       %Lambda reduction factor for successful updates
    fail     = 9;       %Lambda increase factor for unsuccessful updates
    lambda   = 1;       %Starting lambda
    j_freq   = 10;       %How many successful updates before updating J
    j_change = j_freq;       %Build J when equal to J frequency
    thresh   = 0.15;
    
    stuck_check = 0;
    
    %Pull out struct variables
    N0 = tumor.N(:,:,:,1);
    N_true = tumor.N(:,:,:,2:ntp_cal+1);
    if(ntp_pred==1)
        N_true_pred = tumor.N(:,:,:,end);
    end
    t = tumor.t_scan(2:1+ntp_cal);
    ROI = tumor.ROI;
    idx_ROI = find(ROI);
    bcs = tumor.bcs;
    
    tx_params.C = tumor.AUC;
    tx_params.beta1 = tumor.beta1;
    tx_params.beta2 = tumor.beta2;
    tx_params.txduration = tumor.t_trx;
    
    h = tumor.h;
    
    %Pull out parameter bounds and set initial guesses
    kp_up  = bounds.kp_bounds(end) / (delta^2);
    kp_low = bounds.kp_bounds(1) * (delta^2);
    d_up   = bounds.d_bounds(end) / (delta_d^2);
    d_low  = bounds.d_bounds(1) * (delta_d^2);
%     alpha_up  = bounds.alpha_bounds(end) / (delta^2);
%     alpha_low = bounds.alpha_bounds(1) * (delta^2);

%     kp_g       = (exp(log(kp_low) + (log(kp_up)-log(kp_low)) * rand(1,1))) * ones(size(N0));
    kp_g       = kp_up/2 * ones(size(N0));
    kp_g(~ROI) = 0;
    
%     d_g        = exp(log(d_low) + (log(d_up)-log(d_low)) * rand(1,1));
    d_g      = d_up/5;
    
%     alpha1_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));
%     alpha2_g   = exp(log(alpha_low) + (log(alpha_up)-log(alpha_low)) * rand(1,1));    
%     alpha1_g = alpha_low * 5;
%     alpha2_g = alpha_low * 5;
%     alpha1_g = alpha_up/2;
%     alpha2_g = alpha_up/2;
    
    num_kp = numel(idx_ROI);
    num_p  = num_kp + 1;
    
    %Initialize SSE
    [N_sim, TC] = RXDIF_3D(N0, kp_g, d_g, t, h, dt, bcs);
    SSE   = sum((N_true - N_sim).^2,'all');
    
    if(isempty(gcp('nocreate')))
        if(workers == 0)
            myCluster = parcluster('Processes');
            parpool(myCluster.NumWorkers);
        else
            parpool(workers);
        end
    end
    
    %Storage Variables
    kp_store = kp_g(:);
    d_store = d_g;
%     alpha1_store = alpha1_g;
%     alpha2_store = alpha2_g;
    SSE_store = SSE;
    lambda_store = lambda;
    
    %% Calibration loop
    %Calibration Loop
    iteration = 1;
    cnt = 5;
    while(iteration<=max_it && SSE>e_tol)
        %Build Jacobian
        if(j_change==j_freq)
            J = zeros(numel(N_true), num_p);
            
%             start = tic;
            %Local Params
            parfor(i = 1:num_kp)
                kp_p = kp_g;
                kp_p(idx_ROI(i)) = kp_p(idx_ROI(i))*delta;
                dif_kp = kp_p(idx_ROI(i)) - kp_g(idx_ROI(i));
                N_kp    = RXDIF_3D(N0, kp_p, d_g, t, h, dt, bcs);

                J(:,i) = reshape(N_kp - N_sim, [], 1)/dif_kp;
            end
%             stop = toc(start);
%             if(iteration==1)
%                 disp(['kp full perturb FOM = ',num2str(stop),' sec']);
%             end
            
            %Global Params
            d_p  = d_g*delta;
            dif_d  = d_p - d_g;
%             alpha1_p = alpha1_g*delta;
%             dif_alpha1 = alpha1_p - alpha1_g;
%             alpha2_p = alpha2_g*delta;
%             dif_alpha2 = alpha2_p - alpha2_g;

            N_d     = RXDIF_3D(N0, kp_g, d_p, t, h, dt, bcs);
%             N_alpha1 = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_p, alpha2_g, tx_params, t, h, dt, bcs);
%             N_alpha2 = RXDIF_2D_wAC(N0, kp_g, d_g, alpha1_p, alpha2_g, tx_params, t, h, dt, bcs);
            
            J(:,num_kp+1) = reshape(N_d - N_sim, [], 1)/dif_d;
%             J(:,num_kp+2) = reshape(N_alpha1 - N_sim, [], 1)/dif_alpha1;
%             J(:,num_kp+3) = reshape(N_alpha2 - N_sim, [], 1)/dif_alpha2;
            
            j_change = 0;
%             disp('Built J');
        end
        
        
        
        %Calculate update with current regularization
        residuals = reshape(N_true - N_sim, [], 1);
        [update,flags] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*residuals),1e-8,100);
        
        %Update parameters to test
        kp_new = kp_g(idx_ROI) + update(1:num_kp);
        kp_new(kp_new<kp_low) = kp_g(idx_ROI(kp_new<kp_low)) - ((kp_g(idx_ROI(kp_new<kp_low)) - kp_low)./2);
        kp_new(kp_new>kp_up)  = kp_g(idx_ROI(kp_new>kp_up)) + ((kp_up - kp_g(idx_ROI(kp_new>kp_up)))./2);
        kp_test = zeros(size(kp_g));
        kp_test(idx_ROI) = kp_new;

        d_test  = d_g + update(num_kp+1);
        if(d_test<d_low)
            d_test = d_g - (d_g - d_low)/2;
        elseif(d_test>d_up)
            d_test = d_g + (d_up - d_g)/2;
        end
        
%         alpha1_test  = alpha1_g + update(num_kp+2);
%         if(alpha1_test<alpha_low)
%             alpha1_test = alpha1_g - (alpha1_g - alpha_low)/2;
%         elseif(alpha1_test>alpha_up)
%             alpha1_test = alpha1_g + (alpha_up - alpha1_g)/2;
%         end
%         
%         alpha2_test  = alpha2_g + update(num_kp+3);
%         if(alpha2_test<alpha_low)
%             alpha2_test = alpha2_g - (alpha2_g - alpha_low)/2;
%         elseif(alpha2_test>alpha_up)
%             alpha2_test = alpha2_g + (alpha_up - alpha2_g)/2;
%         end
        
        %Test SSE calculation
        [N_test, TC_test] = RXDIF_3D(N0, kp_test, d_test, t, h, dt, bcs);
        SSE_test = sum((N_true - N_test).^2,'all');
        
        %Evaluate new error
        if(SSE_test<SSE) %Successful updates
            %Save new parameters
            N_sim = N_test;
            kp_g = kp_test;
            d_g  = d_test;
%             alpha1_g = alpha1_test;
%             alpha2_g = alpha2_test;
            TC = TC_test;
            
            %Check for convergence
            if(SSE-SSE_test < e_conv)
                disp(['FOM algorithm converged on iteration: ',num2str(iteration)]);
                break;
            end
            SSE = SSE_test;
            
            %Check for tolerance goal
            if(SSE<e_tol)
                disp(['FOM tolerance reached on iteration: ',num2str(iteration)]);
                break;
            end
            lambda = lambda/pass;
            j_change = j_change+1;
            
            stuck_check = 0;
            
        else
            lambda = lambda*fail;
            if(lambda>1e20)
                lambda = 1e-20;
                j_change = j_freq;
                if(stuck_check == 1)
                    disp(['FOM algorithm stuck on iteration: ',num2str(iteration)]);
                    break;
                else
                    stuck_check = 1;
                end
            end
        end
        iteration = iteration + 1;
        
        kp_store = [kp_store, kp_g(:)];
        d_store   = [d_store, d_g];
%         alpha1_store = [alpha1_store, alpha1_g];
%         alpha2_store = [alpha2_store, alpha2_g];
        SSE_store = [SSE_store, SSE];
        lambda_store = [lambda_store, lambda];
        
        if(mod(iteration,5)==0)
%             disp(['Iteration ', num2str(cnt), ' complete']);
            cnt = cnt + 5;
        end
    end
    
    %% Make prediction (if needed)
    if(ntp_pred ~= 0)
        N_pred = RXDIF_3D(N0, kp_g, d_g, tumor.t_scan(end), h, dt, bcs);
        outputs.Pred = N_pred;
    end
    
    %% Finalize results
    %Save calibrated parameter values
    params.kp = kp_g;
    params.d  = d_g;
%     params.alpha1 = alpha1_g;
%     params.alpha2 = alpha2_g;
    
    %Calculate and save fit statistics (unfinished, add rest later)
    %CCC
    CCC_cal = zeros(1,ntp_cal);
    for i = 1:ntp_cal
        CCC_cal(1,i) = CCC_overlap(N_true(:,:,:,i), N_sim(:,:,:,i), thresh);
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
        DICE_cal(1,i) = DICE_calc(N_true(:,:,:,i), N_sim(:,:,:,i), thresh);
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
        Vol_err_cal(1,i) = getVolError(N_true(:,:,:,i), N_sim(:,:,:,i), tumor.h, tumor.dz, thresh);
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
        Cell_err_cal(1,i) = getCellError(N_true(:,:,:,i), N_sim(:,:,:,i), thresh);
    end
    if(ntp_pred ~= 0)
        Cell_err_pred = getCellError(N_true_pred, N_pred, thresh);
    else
        Cell_err_pred = [];
    end
    stats.Cell_err_cal  = Cell_err_cal;
    stats.Cell_err_pred  = Cell_err_pred;
    
    %Save final images
    outputs.Sim = N_sim;
    outputs.TC  = TC;
    outputs.iteration = iteration;
    
    %Plotting figure
    %Plotting figure
    fig = figure;
    subplot(2,4,1)
    plot(1:iteration,kp_store);
    xlabel('Iteration'); ylabel('Proliferation Rate (1/day)');
    subplot(2,4,2)
    plot(1:iteration, d_store);
    xlabel('Iteration'); ylabel('Diffusivity (mm^2/day)');
%     subplot(2,4,3)
%     plot(1:iteration, alpha1_store);
%     xlabel('Iteration'); ylabel('Alpha 1 (1/day)');
%     subplot(2,4,4)
%     plot(1:iteration, alpha2_store);
%     xlabel('Iteration'); ylabel('Alpha 2 (1/day)');
    subplot(2,4,5)
    plot(1:iteration, SSE_store);
    xlabel('Iteration'); ylabel('SSE');
    subplot(2,4,6)
    plot(1:iteration, lambda_store);
    xlabel('Iteration'); ylabel('Lambda');
%     savefig(f1, [home,'/Results/',erase(tumor.name,'.mat'),'_FOM_Calibration_local2D']);
    
end