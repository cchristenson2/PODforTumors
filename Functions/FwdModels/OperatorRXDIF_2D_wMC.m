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

function [N_sim, TC] = OperatorRXDIF_2D_wMC(N0, A, d, B, H, t, dt, M, E, nu, matX, matY, matX_r, matY_r, matXY, V, Vs, freq, reduced)

    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    
%     freq = 25;
    
%     [sy,sx,sz,n] = size(bcs);
%     if(n==1)
%         sz = 1;
%     end
    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    A0 = A;
    
    
    % Initialize damped diffusivity
    if(reduced==0)
        damper = get_damper(matX, matY, N0, M, E, nu);
        S = damper(:);
    else
        N_full = V*N0(:);
        grad_N = Vs' * (matXY * [N_full(:); N_full(:)]);
        damper = get_damper_reduced(matX, matY, grad_N, M, E, nu, Vs);
        S = damper(:);
        dif = d.*ones(numel(S),1) - d.*S;
        idx = find(dif>0);
        for i = 1:length(idx)
            A = A + V(idx(i),:)' * dif(idx(i)) * V(idx(i),:);
        end
    end
    
    for k = 2:nt
        
%         disp(size(matX_r));
%         disp(size(N(:,k-1)));
%         disp(size(d));
%         disp(size(S));
%         disp(size(V');
        
        if(reduced==0)
            X_dot = (matX * N(:,k-1)) .* (matX * (d.*S));
            Y_dot = (matY * N(:,k-1)) .* (matY * (d.*S));
            N(:,k) = N(:,k-1) + dt*(S.*(A*N(:,k-1)) + (X_dot + Y_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
        else
%             disp(size(matX_r));
%         disp(size(N(:,k-1)));
%         disp(size(d));
%         disp(size(S));
%         disp(size(V'));
            X_dot = (matX_r * N(:,k-1)) .* (matX_r * (V'*(d.*S)));
            Y_dot = (matY_r * N(:,k-1)) .* (matY_r * (V'*(d.*S)));
            N(:,k) = N(:,k-1) + dt*(A*N(:,k-1) + (X_dot + Y_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
        end
        
        
%         N(:,k) = N(:,k-1) + dt*(S.*(A*N(:,k-1)) + (X_dot + Y_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
        
        if mod(k, freq) == 0
            if(reduced==0)
                damper = get_damper(matX, matY, N(:,k), M, E, nu);
                S = damper(:);
            else
                N_full = V*N(:,k);
                grad_N = Vs' * (matXY * [N_full(:); N_full(:)]);
                damper = get_damper_reduced(matX, matY, grad_N, M, E, nu, Vs);
                S = damper(:);
                dif = d.*ones(numel(S),1) - d.*S;
                A = A0;
                idx = find(dif>0);
                for i = 1:length(idx)
                    A = A + V(idx(i),:)' * dif(idx(i)) * V(idx(i),:);
                end
            end
        end
        
    end
    
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end