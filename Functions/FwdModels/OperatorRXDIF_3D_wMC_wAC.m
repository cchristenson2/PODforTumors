%{
Forward model for the reaction diffusion equation with AC treatment in matrix form
dn/dt = D*d2n/dx2 + kN(1-N) - alpha1*exp(-beta1*t)*C*N - alpha2*exp(-beta2*t)*C*N

Input:
    - Cell map from start point
    - A operator (diffusivity)
    - B operator (proliferation Lin)
    - H operator (proliferation Quad)
    - T1 operator (A treatment)
    - T2 operator (C treatment)
    - Treatment parameters
        - beta1 for A
        - beta2 for C
        - treatment duration
    - Time spacing
    - Output times

Output:
    - Cell map at desired times
    - Full cell map time course

Contributors: Chase Christenson
%}

function [N_sim, TC] = OperatorRXDIF_3D_wMC_wAC(N0, A, d, B, H, T1, T2, tx_params, t, dt, M, E, nu, matX, matY, matZ, matX_r, matY_r, matZ_r, V, Vs, reduced, bcs, h, dz)

    freq = 25;

    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    nt_trx = tx_params.txduration./dt; %Indices of treatment times
    trx_on = 0; %Starts off, turns on at first treatment delivery
    trx_cnt = 1;

    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    b1 = tx_params.beta1;
    b2 = tx_params.beta2;
    
    
    % Initialize damped diffusivity
    if(reduced==0)
        damper = get_damper_3D(matX, matY, matZ, N0, M, E, nu);
        S = damper(:);
    else
        N_full = V*N0(:);
        grad_N = Vs' * [matX*N_full(:); matY*N_full(:); matZ*N_full(:)];
        damper = get_damper_reduced_3D(matX, matY, matZ, grad_N, M, E, nu, Vs);
        S = damper(:);
        temp_A = assembleA(bcs(:,:,:,1), d.*S, h, dz, bcs);
        A = V' * temp_A * V;
    end
    
    for k = 2:nt
        %Get time since last treatment
        if(k-1>nt_trx(trx_cnt))
            if(trx_on==0)
                t = 0;
            end
            trx_on = 1;
            if(trx_cnt < numel(nt_trx))
                if(k-1>nt_trx(trx_cnt+1))
                    trx_cnt=trx_cnt+1;
                    t=0;
                end
                t=t+dt;
            else
                t=t+dt;
            end
        end
        %Solve treatment effects
        if(trx_on~=0)
            treat = ((T1.*exp(-b1*t)) + (T2.*exp(-b2*t)))*(N(:,k-1));
        else
            treat = 0;
        end
        
        if(reduced==0)
            X_dot = (matX * N(:,k-1)) .* (matX * (d.*S));
            Y_dot = (matY * N(:,k-1)) .* (matY * (d.*S));
            Z_dot = (matZ * N(:,k-1)) .* (matZ * (d.*S));
            N(:,k) = N(:,k-1) + dt*(S.*(A*N(:,k-1)) + (X_dot + Y_dot + Z_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)) - treat);
        else
%             disp(size(matX_r));
%         disp(size(N(:,k-1)));
%         disp(size(d));
%         disp(size(S));
%         disp(size(V'));
            X_dot = (matX_r * N(:,k-1)) .* (matX_r * (V'*(d.*S)));
            Y_dot = (matY_r * N(:,k-1)) .* (matY_r * (V'*(d.*S)));
            Y_dot = (matZ_r * N(:,k-1)) .* (matZ_r * (V'*(d.*S)));
            N(:,k) = N(:,k-1) + dt*(A*N(:,k-1) + (X_dot + Y_dot + Z_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)) - treat);
        end
        
        
%         N(:,k) = N(:,k-1) + dt*(S.*(A*N(:,k-1)) + (X_dot + Y_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
        
        if mod(k, freq) == 0
            if(reduced==0)
                damper = get_damper_3D(matX, matY, matZ, N0, M, E, nu);
                S = damper(:);
            else
                N_full = V*N0(:);
                grad_N = Vs' * [matX*N_full(:); matY*N_full(:); matZ*N_full(:)];
                damper = get_damper_reduced_3D(matX, matY, matZ, grad_N, M, E, nu, Vs);
                S = damper(:);
                temp_A = assembleA(bcs(:,:,:,1), d.*S, h, dz, bcs);
                A = V' * temp_A * V;
            end
        end
        
    end
    
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end