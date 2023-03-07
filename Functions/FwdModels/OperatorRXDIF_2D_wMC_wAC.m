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

function [N_sim, TC] = OperatorRXDIF_2D_wMC_wAC(N0, d, B, H, T1, T2, tx_params, t, dt, M, E, nu, matX, matY, matX_r, matY_r, V, Vs, Vd, Ar_lib, k, reduced)

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
        damper = get_damper(matX, matY, N0, M, E, nu);
        S = damper(:);
    else
        N_full = V*N0(:);
        grad_N = Vs' * ([matX * N_full(:); matY * N_full(:)]);
        damper = get_damper_reduced(matX, matY, grad_N, M, E, nu, Vs);
        S = damper(:);
        A = OperatorInterp_local(Vd'*S, Ar_lib, k);
        A = d.*A;
    end
    
    for l = 2:nt
        %Get time since last treatment
        if(l-1>nt_trx(trx_cnt))
            if(trx_on==0)
                t = 0;
            end
            trx_on = 1;
            if(trx_cnt < numel(nt_trx))
                if(l-1>nt_trx(trx_cnt+1))
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
            treat = (T1.*exp(-b1*t) + T2.*exp(-b2*t))*(N(:,l-1));
        else
            treat = 0;
        end
        
        if(reduced==0)
            X_dot = (matX * N(:,l-1)) .* (matX * (d.*S));
            Y_dot = (matY * N(:,l-1)) .* (matY * (d.*S));
            N(:,l) = N(:,l-1) + dt*(S.*(A*N(:,l-1)) + (X_dot + Y_dot) + B*N(:,l-1) - H*kron(N(:,l-1), N(:,l-1)) - treat);
        else
            X_dot = (matX_r * N(:,l-1)) .* (matX_r * (V'*(d.*S)));
            Y_dot = (matY_r * N(:,l-1)) .* (matY_r * (V'*(d.*S)));
            N(:,l) = N(:,l-1) + dt*(A*N(:,l-1) + (X_dot + Y_dot) + B*N(:,l-1) - H*kron(N(:,l-1), N(:,l-1)) - treat);
        end
        
        
        if mod(l, freq) == 0
            if(reduced==0)
                damper = get_damper(matX, matY, N(:,l), M, E, nu);
                S = damper(:);
            else
                N_full = V*N0(:);
                grad_N = Vs' * ([matX * N_full(:); matY * N_full(:)]);
                damper = get_damper_reduced(matX, matY, grad_N, M, E, nu, Vs);
                S = damper(:);
                A = OperatorInterp_local(Vd'*S, Ar_lib, k);
                A = d.*A;
            end
        end
        
    end
    
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end