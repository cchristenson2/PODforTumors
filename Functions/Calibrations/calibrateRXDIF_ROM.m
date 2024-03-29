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

function [params,N_cal] = calibrateRXDIF_ROM(initial, data, library, bounds, t, dt, model, tx_params)
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
    kp_guess = mean(bounds(2,:));
    if model == 2
        alpha_guess = mean(bounds(3,:));
    end
    
    %Get initial operators for forward model
    A_guess = OperatorInterp(d_guess,library.Ar_lib);
    B_guess = OperatorInterp(kp_guess,library.Br_lib);
    H_guess = OperatorInterp(kp_guess,library.Hr_lib);
    if model == 2
        T_guess = OperatorInterp(alpha_guess,library.Tr_lib);
    end
    
    %Initialize sum squared errors for minimization
    if model == 1
        N_guess = OperatorRXDIF_2D(initial, A_guess, B_guess, H_guess, t, dt);
        num_params = 2;
    elseif model == 2
        N_guess = OperatorRXDIF_2D_wAC_comb(initial, A_guess, B_guess, H_guess, T_guess, tx_params, t, dt);
        num_params = 3;
    else
        error('Error. Model identifier must be 1 or 2');
    end
    
    SSE = sum((data - N_guess).^2,'all');
    
    %Begin calibration loop
    iter = 1;
    while(iter<=max_it && SSE>e_tol)
        %Build jacobian numerically by perturbing each parameter individually
        J = zeros(num, num_params);
        
        d_perturb  = d_guess*delta;
        d_dif = d_perturb - d_guess;
        A_perturb = OperatorInterp(d_perturb,library.Ar_lib);
        
        kp_perturb = kp_guess*delta;
        kp_dif = kp_perturb - kp_guess;
        B_perturb = OperatorInterp(kp_perturb,library.Br_lib);
        H_perturb = OperatorInterp(kp_perturb,library.Hr_lib);
        
        if model == 1
            N_d = OperatorRXDIF_2D(initial, A_perturb, B_guess, H_guess, t, dt);
            N_kp = OperatorRXDIF_2D(initial, A_guess, B_perturb, H_perturb, t, dt);
        else % model = 2
            alpha_perturb = alpha_guess*delta;
            alpha_dif = alpha_perturb - alpha_guess;
            T_perturb = OperatorInterp(alpha_perturb,library.Tr_lib);
            
            N_d = OperatorRXDIF_2D_wAC_comb(initial, A_perturb, B_guess, H_guess, T_guess, tx_params, t, dt);
            N_kp = OperatorRXDIF_2D_wAC_comb(initial, A_guess, B_perturb, H_perturb, T_guess, tx_params, t, dt);
            N_alpha = OperatorRXDIF_2D_wAC_comb(initial, A_guess, B_guess, H_guess, T_perturb, tx_params, t, dt);
            
            J(:,3) = (N_alpha(:) - N_guess(:))./alpha_dif;
        end
        J(:,1) = (N_d(:) - N_guess(:))./d_dif;
        J(:,2) = (N_kp(:) - N_guess(:))./kp_dif;
        
        
        %Use jacobian to calculate current update
        [update,~] = bicgstab((J'*J + lambda*diag(diag(J'*J))),(J'*(data(:) - N_guess(:))),1e-8,100);
        
        d_test  = d_guess + update(1);
        if(d_test<bounds(1,1))
            d_test = d_guess - (d_guess - bounds(1,1))/2;
        elseif(d_test>bounds(1,2))
            d_test = d_guess + (bounds(1,2) - d_guess)/2;
        end
        
        kp_test = kp_guess + update(2);
        if(kp_test<bounds(2,1))
            kp_test = kp_guess - (kp_guess - bounds(2,1))/2;
        elseif(kp_test>bounds(2,2))
            kp_test = kp_guess + (bounds(2,2) - kp_guess)/2;
        end
        
        if model == 2
            alpha_test = alpha_guess + update(3);
            if(alpha_test<bounds(3,1))
                alpha_test = alpha_guess - (alpha_guess - bounds(3,1))/2;
            elseif(alpha_test>bounds(3,2))
                alpha_test = alpha_guess + (bounds(3,2) - alpha_guess)/2;
            end
        end
        
        %Get operators for test parameters
        A_test = OperatorInterp(d_test,library.Ar_lib);
        B_test = OperatorInterp(kp_test,library.Br_lib);
        H_test = OperatorInterp(kp_test,library.Hr_lib);
        if model == 2
            T_test = OperatorInterp(alpha_test,library.Tr_lib);
        end
        
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
            kp_guess = kp_test; B_guess = B_test; H_guess = H_test;
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
    if model == 2
        params.alpha = alpha_guess;
    end

end

