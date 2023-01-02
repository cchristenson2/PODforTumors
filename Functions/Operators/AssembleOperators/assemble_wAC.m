%{
Build all operators required for A/C treatment forward model

Input:
    - number of elements in domain
    - diffusivity
    - proliferation
    - treatment 1 efficacy
    - treatment 2 efficacy
    - grid spacing
    - boundary map

Output:
    - A operator
    - B operator
    - H operator
    - T1 operator
    - T2 operator

Contributors: Chase Christenson
%}

function [A,B,H,T1,T2] = assemble_wAC(N, D, kp, alpha1, alpha2, h, dz, bcs)
    A = assembleA(N, D, h, dz, bcs);
    
    n = numel(N);
    B = assembleB(n, kp);
    H = assembleH(n, kp);
    T1 = assembleT(n, alpha1);
    T2 = assembleT(n, alpha2);
end