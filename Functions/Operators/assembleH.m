%{
Build operator for quadratic proliferation

dN/dt = H*kron(N,N) = k*N^2

Input:
    - N; Cell map for sizing, 2D or 3D cell map 
    - k; proliferation rate
        - single parameter or local vector with same size as n

Outputs:
    - H operator

Contributors: Graham Pash, Chase Christenson
%}
function [H] = assembleH(N, k)
    n = numel(N);
    H = sparse(n, n^2);
    H(:, 1:n+1:end) = eye(n, 'like', sparse(n,n));
    H = k(:).*sparse(H);
end