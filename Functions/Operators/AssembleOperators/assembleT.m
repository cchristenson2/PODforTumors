%{
Build operator for treatment effects

Input:
    - number of elements in domain
    - alpha for treatment effect rate

Output:
    - Treatment operator

Note:
    - Concentration map applied to T before input into forward model

Contributors: Graham Pash
%}
function [T] = assembleT(n, alpha)
    T = sparse(alpha*eye(n, 'like', sparse(n,n)));
end