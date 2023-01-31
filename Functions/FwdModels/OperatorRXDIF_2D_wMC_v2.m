%{
Forward model for the reaction diffusion equation
dn/dt = D*d2n/dx2 + kN(1-N)

Input:
    - Cell map from start point
    - A operator (diffusivity)
    - B operator (proliferation Lin)
    - H operator (proliferation Quad)
    - Time spacing
    - Output times

Output:
    - Cell map at desired times
    - Full cell map time course

Contributors: Chase Christenson, Graham Pash
%}

function [N_sim, TC] = OperatorRXDIF_2D_wMC_v2(N0, A1, grad, d, B, H, t, dt, h, dz, bcs, M, E, nu, matX, matY, lambda, freq)
    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    
    [sy,sx,sz,n] = size(bcs);
    if(n==1)
        sz = 1;
    end
    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    % Initialize damped diffusivity
    damper = get_damper(matX, matY, N0, M, E, nu, h, lambda);
    S = damper(:);
    
    for k = 2:nt
        
%         disp(size(S));
%         disp(size(A1));
%         disp(size(N(:,k-1)));
%         disp(size(A1*N(:,k-1)));
%         
%         disp(size(d.*S));
%         disp(size(grad*(d.*S)));
%         disp(size(grad*N(:,k-1)));
%         
%         disp(size(B*N(:,k-1)));
        
%         one = S.*(A1*N(:,k-1));
%         two = ((grad*(d.*S)).*(grad*N(:,k)));
%         three = B*N(:,k-1);
%         four = H*kron(N(:,k-1), N(:,k-1));
        
        
%         N(:,k) = N(:,k-1) + dt*(S.*(A1*N(:,k-1)) + ((grad*(d.*S)).*(grad*N(:,k-1))) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
        
        X_dot = (matX * N(:,k-1)) .* (matX * (d.*S));
        Y_dot = (matY * N(:,k-1)) .* (matY * (d.*S));
        N(:,k) = N(:,k-1) + dt*(S.*(A1*N(:,k-1)) + (X_dot + Y_dot) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
%         N(:,k) = N(:,k-1) + dt*(one + two + three + four);
        
        
        if mod(k, freq) == 0
            damper = get_damper(matX, matY, reshape(N(:,k),sy,sx,sz), M, E, nu, h, lambda);
            S = damper(:);
        end
    end
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end