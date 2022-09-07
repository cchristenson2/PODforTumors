%{ 
Builds ROM library for a set of parameter bounds, and grid sizes
Based on forward model with A/C treatment

Inputs:
    - parameter bounds
    - number of grid points
    - grid spacing
    - naming format

Outputs:
    - library for all operators

Contributors: Chase Christenson
%}

function [A_lib, B_lib, H_lib, T_lib] = buildOperatorLibrary(bounds,n,h,fmt,bcs)

%     %Add function paths
%     addpath(genpath('C:/ROM/'))
    
    A_lib = struct;
    B_lib = struct;
    H_lib = struct;
    T_lib = struct;
    
    len = numel(bounds.d_bounds);
    for i = 1:len
        %Diffusivity based operators
        d = bounds.d_bounds(i);
        d_str = replace(num2str(d,fmt),'.','_');
        name = ['d',d_str];
        A = assembleA(n, d, h, bcs);
        eval(['A_lib.',name,'_A = A;']);

        %Proliferation based operators
        kp = bounds.kp_bounds(i);
        kp_str = replace(num2str(kp,fmt),'.','_');
        name = ['kp',kp_str];
        B = assembleB(n,kp);
        H = assembleH(n,kp);
        eval(['B_lib.',name,'_B = B;']);
        eval(['H_lib.',name,'_H = H;']);
        
        %Alpha based operators
        alpha = bounds.alpha_bounds(i);
        alpha_str = replace(num2str(alpha,fmt),'.','_');
        name = ['alpha',alpha_str];
        T = assembleT(n, alpha);
        eval(['T_lib.',name,'_T = T;']);
    end
end