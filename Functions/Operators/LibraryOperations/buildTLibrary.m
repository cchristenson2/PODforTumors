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

function [T_lib] = buildTLibrary(bounds,N,h,fmt,bcs,dz)

%     %Add function paths
%     addpath(genpath('C:/ROM/'))
    
    T_lib = struct;
    
    len_a = numel(bounds.alpha_bounds);
    
    n = numel(N);

    
    for i = 1:len_a
        %Alpha based operators
        alpha = bounds.alpha_bounds(i);
        alpha_str = replace(num2str(alpha,fmt),'.','_');
        name = ['alpha',alpha_str];
        eval(['T_lib.',name,'_T = assembleT(n, alpha);']);
    end

end