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
function [A,B,H] = assemble_noAC(n, D, kp, h, bcs)
    A = assembleA(n, D, h, bcs);
    B = assembleB(n, kp);
    H = assembleH(n, kp);
end