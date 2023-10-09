%{
Build operator for linear proliferation term

dN/dt = BN = k*N

Input:
    - N; Cell map for sizing, 2D or 3D cell map 
    - k; proliferation rate
        - single parameter or local vector with same size as n

Outputs:
    - B operator

Contributors: Graham Pash, Chase Christenson
%}
function [B] = assembleB(N, k)
    n = numel(N);
    B = sparse(k(:).*eye(n, 'like', sparse(n,n)));
end