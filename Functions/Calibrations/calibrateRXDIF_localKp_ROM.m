%{ 
Calibration for 2D RXDIF model for tumor response

M1 - dn/dt = D*d^2n/h^2 + kp*N(1-N)
M2 - dn/dt = D*d^2n/h^2 + kp*N(1-N) - alpha*C*N*exp(-beta*t)

Input:
    - initial; Cell map from start point
    - data; Cell maps we calibrate to
    - library; Library of reduced operators for each parameter type
    - bounds; parameter bounds, [2 X 2] or [3 X 2] for model 1 and 2 respectively
    - t; time of data acquisition
    - dt; Time spacing
    - model; 1 or 2 depending on which model is desired
    - tx_params; only needed if model = 2

Output:
    - params; 
    - N_cal; Calibrated cell map at same time as data was acquired

Contributors: Chase Christenson

Based off reference:
Hormuth, David A. et al. 
“Mechanically Coupled Reaction-Diffusion Model to Predict Glioma Growth: 
Methodological Details.” Cancer Systems Biology 1711 (2018): 225–241.

%}

function [params,N_cal] = calibrateRXDIF_localKp_ROM(initial, data, library, bounds, kp_mode_bounds, t, dt, model, tx_params, V)
    %Base parameter for levenberg-marquardt non-linear least squares parameter updating
    %LM; (J^T*J + alpha*(diag(J^T*J))) * delta_B = J*T * residuals
    e_tol    = 1e-6;    %Calibration SSE goal
    e_conv   = 1e-7;    %Minimum change in SSE for an iteration
    max_it   = 500;     %Maximum iterations
    delta    = 1.001;   %Perturbation magnitude
    pass     = 7;       %Lambda reduction factor for successful updates
    fail     = 9;       %Lambda increase factor for unsuccessful updates
    lambda   = 1;       %Starting lambda
    num = numel(data);  %Number of data points
    
    %Initial guess for the problem, centered on bounds
    d_guess = mean(bounds(1,:));
    if model == 2
        alpha_guess = mean(bounds(3,:));
    end
    
    %Get initial operators for forward model
    A_guess = OperatorInterp(d_guess,library.Ar_lib);
    
    if model == 2
        T_guess = OperatorInterp(alpha_guess,library.Tr_lib);
    end
    
    %Initial guess and operator for reduced parameters
    [n,k] = size(V); idx = find(sum(V,2));
    kp_target = mean(bounds(2,:));
    if(kp_target == 0)
        kp_target = 1e-6;
    end
    kp_temp = zeros(n,1);
    kp_temp(idx) = kp_target;
    kp_temp_r = V' * kp_temp(:);
    
    A = [V(idx,:); -1.*V(idx,:)];
    B = [bounds(2,2).*ones(numel(idx),1); -1.*bounds(2,1).*ones(numel(idx),1)];

    kp_guess_r = lsqlin(V(idx,:), kp_temp(idx), A, B, [], [], kp_mode_bounds(:,1), kp_mode_bounds(:,2), kp_temp_r,...
            optimoptions('lsqlin','Algorithm','active-set', 'ConstraintTolerance', 0, 'OptimalityTolerance', 1e-8, 'MaxIterations', 20,'Display','off'));
    kp_guess = V * kp_guess_r(:);
    clear kp_target kp_temp kp_temp_r
    
    B_guess = OperatorInterp_local(kp_guess_r,library.Br_lib);
    H_guess = OperatorInterp_local(kp_guess_r,library.Hr_lib);
    
    %Initialize sum squared errors for minimization
    if model == 1
        N_guess = OperatorRXDIF_2D(initial, A_guess, B_guess, H_guess, t, dt);
        num_params = k + 1;
    elseif model == 2
        N_guess = OperatorRXDIF_2D_wAC_comb(initial, A_guess, B_guess, H_guess, T_guess, tx_params, t, dt);
        num_params = k + 2;
    else
        error('Error. Model identifier must be 1 or 2');
    end
    
    SSE = sum((data - N_guess).^2,'all');
    
    %Begin calibration loop
    iter = 1;
    while(iter<=max_it && SSE>e_tol)
        %Build jacobian numerically by perturbing each parameter individually
        J = zeros(num, num_params);
        
        %Global parameter perturbation
        d_perturb  = d_guess*delta;
        d_dif = d_perturb - d_guess;
        A_perturb = OperatorInterp(d_perturb,library.Ar_lib);
        
        if model == 1
            N_d = OperatorRXDIF_2D(initial, A_perturb, B_guess, H_guess, t, dt);
        else % model = 2
            alpha_perturb = alpha_guess*delta;
            alpha_dif = alpha_perturb - alpha_guess;
            T_perturb = OperatorInterp(alpha_perturb,library.Tr_lib);
            
            N_d = OperatorRXDIF_2D_wAC_comb(initial, A_perturb, B_guess, H_guess, T_guess, tx_params, t, dt);
            N_alpha = OperatorRXDIF_2D_wAC_comb(initial, A_guess, B_guess, H_guess, T_perturb, tx_params, t, dt);
            
            J(:,end) = (N_alpha(:) - N_guess(:))./alpha_dif;
        end
        J(:,1) = (N_d(:) - N_guess(:))./d_dif;
        
        %Local parameter perturbation
        for i = 1:k
            kp_perturb_r = kp_guess_r;
            kp_perturb_r(i) = kp_perturb_r(i)*delta;
            kp_dif = kp_perturb_r(i) - kp_guess_r(i);
            B_perturb = OperatorInterp_local(kp_perturb_r,library.Br_lib);
            H_perturb = OperatorInterp_local(kp_perturb_r,library.Hr_lib);
            
            if model == 1
                N_kp = OperatorRXDIF_2D(initial, A_guess, B_perturb, H_perturb, t, dt);
            else % model = 2
                N_kp = OperatorRXDIF_2D_wAC_comb(initial, A_guess, B_perturb, H_perturb, T_guess, tx_params, t, dt);
            end
            
            J(:,i+1) = (N_kp(:) - N_guess(:))./kp_dif;
        end
        
        
        %Use jacobian to calculate current update
        [update,~] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*(data(:) - N_guess(:))),1e-8,100);
        
        d_test  = d_guess + update(1);
        if(d_test<bounds(1,1))
            d_test = d_guess - (d_guess - bounds(1,1))/2;
        elseif(d_test>bounds(1,2))
            d_test = d_guess + (bounds(1,2) - d_guess)/2;
        end
        
        kp_test_r = kp_guess_r + update(2:1+k);
        for i = 1:k
            if(kp_test_r(i) < kp_mode_bounds(i,1))
                kp_test_r(i) = kp_mode_bounds(i,1);
            elseif(kp_test_r(i) > kp_mode_bounds(i,2))
                kp_test_r(i) = kp_mode_bounds(i,2);
            end
        end
        kp_test = V*kp_test_r;
        
        if(~isempty(find(kp_test<bounds(2,1) & abs(kp_test)>1e-6 ,1)) || ~isempty(find(kp_test>bounds(2,2),1))) %if reconstructed kp is out of bounds
            A = [V(idx,:); -1.*V(idx,:)];
            B = [bounds(2,2).*ones(numel(idx),1); -1.*bounds(2,1).*ones(numel(idx),1)];
            
            kp_test_r = lsqlin(V(idx,:), kp_test(idx), A, B, [], [], kp_mode_bounds(:,1), kp_mode_bounds(:,2), kp_test_r,...
                optimoptions('lsqlin','Algorithm','active-set', 'ConstraintTolerance', 0, 'OptimalityTolerance', 1e-8, 'MaxIterations', 50,'Display','off'));
                
            kp_test = V * kp_test_r(:);
        end
        
        if model == 2
            alpha_test = alpha_guess + update(end);
            if(alpha_test<bounds(3,1))
                alpha_test = alpha_guess - (alpha_guess - bounds(3,1))/2;
            elseif(alpha_test>bounds(3,2))
                alpha_test = alpha_guess + (bounds(3,2) - alpha_guess)/2;
            end
        end
        
        %Get operators for test parameters
        A_test = OperatorInterp(d_test,library.Ar_lib);
        if model == 2
            T_test = OperatorInterp(alpha_test,library.Tr_lib);
        end
        
        B_test = OperatorInterp_local(kp_test_r,library.Br_lib);
        H_test = OperatorInterp_local(kp_test_r,library.Hr_lib);
        
        
        %Test new parameters
        if model == 1
            N_test = OperatorRXDIF_2D(initial, A_test, B_test, H_test, t, dt);
        else %model = 2
            N_test = OperatorRXDIF_2D_wAC_comb(initial, A_test, B_test, H_test, T_test, tx_params, t, dt);
        end
        
        SSE_test = sum((data - N_test).^2,'all');
        
        %Check if new SSE is better
        if SSE_test < SSE
            %Yes so save new parameters and simulation
            d_guess  = d_test; A_guess = A_test;
            kp_guess = kp_test; kp_guess_r = kp_test_r; B_guess = B_test; H_guess = H_test;
            if model == 2
                alpha_guess = alpha_test; T_guess = T_test;
            end
            
            N_guess = N_test;
            
            lambda = lambda/pass;
            
            %Check if model has converged
            if(SSE-SSE_test < e_conv)
                disp(['ROM algorithm converged on iteration: ',num2str(iter)]);
                break;
            end
            
            if SSE_test < e_tol
                disp(['ROM tolerance reached on iteration: ',num2str(iter)]);
            end
            
            SSE = SSE_test;
            
        else %Update is not better            
            lambda = lambda*fail; %Increase lambda
        end
        iter = iter + 1;
    end
    
    %Save best guess simulation
    N_cal = N_guess;
    
    %Save parameters
    params.d = d_guess;
    params.kp = kp_guess;
    params.kp_r = kp_guess_r;
    if model == 2
        params.alpha = alpha_guess;
    end

end

