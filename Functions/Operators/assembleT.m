%{
Build operator for treatment effects

dN/dt = -T*N = -alpha*N

Input:
    - N; Cell map for sizing, 2D or 3D cell map 
    - alpha; treatment effect rate

Output:
    - T; Treatment operator

Notes:
    - Concentration map applied to T before input into forward model

Contributors: Graham Pash
%}
function [T] = assembleT(N, alpha)
    n = numel(N);
    T = sparse(alpha*eye(n, 'like', sparse(n,n)));
end