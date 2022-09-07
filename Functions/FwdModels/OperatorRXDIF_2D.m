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

function [N_sim, TC] = OperatorRXDIF_2D(N0, A, B, H, t, dt)
    nt = t(end)/dt + 1;
    t_ind = t./dt + 1;
    
    [num,~] = size(N0);
    N = zeros(num, nt);
    N(:,1) = N0;
    
    for k = 2:nt
        N(:,k) = N(:,k-1) + dt*(A*N(:,k-1) + B*N(:,k-1) - H*kron(N(:,k-1), N(:,k-1)));
    end
    TC = squeeze(sum(N,1));
    N_sim = N(:,t_ind);
end