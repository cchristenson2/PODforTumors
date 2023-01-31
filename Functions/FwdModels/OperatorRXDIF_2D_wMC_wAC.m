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

Contributors: Chase Christenson, Graham Pash
%}

function [N_sim, TC] = OperatorRXDIF_2D_wMC_wAC(N0, A1, d, B, H, T1, T2, tx_params, t, dt, M, E, nu, matX, matY)
    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    nt_trx = tx_params.txduration./dt; %Indices of treatment times
    trx_on = 0; %Starts off, turns on at first treatment delivery
    trx_cnt = 1;
    
    freq = 25;
    
%     [sy,sx,sz,n] = size(bcs);
%     if(n==1)
%         sz = 1;
%     end
    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    b1 = tx_params.beta1;
    b2 = tx_params.beta2;
    
    % Initialize damped diffusivity
    damper = get_damper(matX, matY, N(:,1), M, E, nu);
    S = damper(:);
    
    for k = 2:nt
        %Get time since last treatment
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
            treat = (T1.*exp(-b1*t) + T2.*exp(-b2*t))*(N(:,k-1));
        else
            treat = 0;
        end
        
        X_dot = (matX * N(:,k-1)) .* (matX * (d.*S));
        Y_dot = (matY * N(:,k-1)) .* (matY * (d.*S));
        N(:,k) = N(:,k-1) + dt*(S.*(A1*N(:,k-1)) + (X_dot + Y_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)) - treat);
        
        if mod(k, freq) == 0
            damper = get_damper(matX, matY, N(:,k), M, E, nu);
            S = damper(:);
        end
    end
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end