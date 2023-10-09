%{
Forward model for the reaction diffusion equation
dn/dt = D*d^2n/h^2 + kp*N(1-N)

Input:
    - N0; Cell map from start point
    - A operator (diffusivity)
    - B operator (proliferation Lin)
    - H operator (proliferation Quad)
    - t; Outputs times
    - dt; Time spacing

Output:
    - N_sim; Cell map at desired times
    - TC; Full cell map time course (do not use in reduced form)

Contributors: Chase Christenson, Graham Pash
%}

function [N_sim, TC] = OperatorRXDIF_2D(N0, A, B, H, t, dt)
    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    for k = 2:nt
        N(:,k) = N(:,k-1) + dt*(A*N(:,k-1) + B*N(:,k-1) - (H*kron(N(:,k-1), N(:,k-1))));
    end
    TC = N;
%     N_full = N_sim;
    N_sim = N(:,t_ind);
end