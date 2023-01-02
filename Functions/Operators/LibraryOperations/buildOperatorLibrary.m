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

function [A_lib, B_lib, H_lib, T_lib] = buildOperatorLibrary(bounds,N,h,fmt,bcs,dz)

%     %Add function paths
%     addpath(genpath('C:/ROM/'))
    
    A_lib = struct;
    B_lib = struct;
    H_lib = struct;
    T_lib = struct;
    
    len_d = numel(bounds.d_bounds);
    len_kp = numel(bounds.kp_bounds);
    
    n = numel(N);
    
    for i = 1:len_d
        %Diffusivity based operators
        d = bounds.d_bounds(i);
        d_str = replace(num2str(d,fmt),'.','_');
        name = ['d',d_str];
        eval(['A_lib.',name,'_A = assembleA(N, d, h, dz, bcs);']);

        %Alpha based operators
        alpha = bounds.alpha_bounds(i);
        alpha_str = replace(num2str(alpha,fmt),'.','_');
        name = ['alpha',alpha_str];
        eval(['T_lib.',name,'_T = assembleT(n, alpha);']);
    end
    
    for i = 1:len_kp
        %Proliferation based operators
        kp = bounds.kp_bounds(i);
        kp_str = replace(num2str(kp,fmt),'.','_');
        kp_str = replace(kp_str,'-','n');
        name = ['kp',kp_str];
        eval(['B_lib.',name,'_B = assembleB(n,kp);']);
        eval(['H_lib.',name,'_H = assembleH(n,kp);']);
    end
end