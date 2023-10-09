%{
Forward model for the reaction diffusion equation with AC treatment in matrix form
dn/dt = D*d^2n/h^2 + kp*N(1-N) - alpha1*(exp(-beta1*t) - exp(-beta2*t))*C*N

Input:
    - N0; Cell map from start point
    - A operator (diffusivity)
    - B operator (proliferation Lin)
    - H operator (proliferation Quad)
    - T operator (A/C treatment)
    - tx_params; Treatment parameters struct
        - beta1 for A
        - beta2 for C
        - treatment duration
    - t; Output times
    - dt; Time spacing

Output:
    - N_sim; Cell map at desired times
    - TC; Full cell map time course

Contributors: Chase Christenson, Graham Pash
%}

function [N_sim, TC] = OperatorRXDIF_2D_wAC_comb(N0, A, B, H, T1, tx_params, t, dt)
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
    
    for k = 2:nt
        %Get time since last treatment
        if(k-1>nt_trx(trx_cnt))
            if(trx_on==0)
                t = 0;
                trx_on = 1;
            else
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
        end
        %Solve treatment effects
        if(trx_on~=0)
            treat = T1.*(exp(-b1*t)+exp(-b2*t))*(N(:,k-1));
        else
            treat = 0;
        end
        
        N(:,k) = N(:,k-1) + dt*(A*N(:,k-1) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)) - treat);
    end
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end