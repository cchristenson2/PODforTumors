%{
Build operator for quadratic proliferation

Input:
    - number of elements in domain
    - proliferation rate
        - single parameter or local vector with same size as n

Outputs:
    - B operator

Contributors: Graham Pash
%}
function [H] = assembleH(n, k)
    H = sparse(n, n^2);
    H(:, 1:n+1:end) = eye(n, 'like', sparse(n,n));
    H = k(:).*sparse(H);
end