%{
Build operator for linear proliferation

Input:
    - number of elements in domain
    - proliferation rate
        - single parameter or local vector with same size as n

Outputs:
    - B operator

Contributors: Graham Pash
%}
function [B] = assembleB(n, k)
    B = sparse(k(:).*eye(n, 'like', sparse(n,n)));
end