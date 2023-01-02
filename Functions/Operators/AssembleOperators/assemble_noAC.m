%{
Build all operators required for no treatment forward model

Input:
    - number of elements in domain
    - diffusivity
    - proliferation
    - grid spacing
    - boundary map

Output:
    - A operator
    - B operator
    - H operator

Contributors: Chase Christenson
%}
function [A,B,H] = assemble_noAC(N, D, kp, h, dz, bcs)
    A = assembleA(N, D, h, dz, bcs);
    
    n = numel(N);
    B = assembleB(n, kp);
    H = assembleH(n, kp);
end