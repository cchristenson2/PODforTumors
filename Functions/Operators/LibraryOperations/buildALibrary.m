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

function [A_lib] = buildALibrary(bounds,N,h,fmt,bcs,dz)

%     %Add function paths
%     addpath(genpath('C:/ROM/'))
    
    A_lib = struct;
    
    len_d = numel(bounds.d_bounds);
    
    for i = 1:len_d
        %Diffusivity based operators
        d = bounds.d_bounds(i);
        d_str = replace(num2str(d,fmt),'.','_');
        name = ['d',d_str];
        eval(['A_lib.',name,'_A = assembleA(N, d, h, dz, bcs);']);
    end
end