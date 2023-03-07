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

function [N_sim, TC] = OperatorRXDIF_2D_wMC(N0, d, B, H, t, dt, M, E, nu, matX, matY, matX_r, matY_r, V, Vs, Vd, Ar_lib, k, reduced)

    freq = 25;

    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    
    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    
    
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
        
        
        if(reduced==0)
            X_dot = (matX * N(:,l-1)) .* (matX * (d.*S));
            Y_dot = (matY * N(:,l-1)) .* (matY * (d.*S));
            N(:,l) = N(:,l-1) + dt*(S.*(A*N(:,l-1)) + (X_dot + Y_dot) + B*N(:,l-1) - H*kron(N(:,l-1), N(:,l-1)));
        else
            X_dot = (matX_r * N(:,l-1)) .* (matX_r * (V'*(d.*S)));
            Y_dot = (matY_r * N(:,l-1)) .* (matY_r * (V'*(d.*S)));
            N(:,l) = N(:,l-1) + dt*(A*N(:,l-1) + (X_dot + Y_dot) + B*N(:,l-1) - H*kron(N(:,l-1), N(:,l-1)));
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