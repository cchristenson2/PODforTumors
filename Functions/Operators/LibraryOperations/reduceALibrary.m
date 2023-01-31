%{ 
Reduce ROM library for a set of parameter bounds, and given projection matrix
Based on forward model with A/C treatment

Inputs:
    - operator libraries
    - projection matrix
    - parameter bounds
    - naming format
    - Treatment params with concentration map (added to T operators prior to reduction)

Outputs:
    - reduced library for all operators

Contributors: Chase Christenson
%}

function [Ar_lib] = reduceALibrary(A_lib, V, bounds)
    fmt = bounds.fmt;
    num_d = numel(bounds.d_bounds);
    
    Ar_lib = struct;
    
    for i = 1:num_d
        d      = bounds.d_bounds(i);
        d_str     = replace(num2str(d,fmt),'.','_');
        
        %d
        name = ['d',d_str];
        eval(['Ar_lib.',name,'_A = V'' * A_lib.',name,'_A * V;']);
        eval(['clear ',name,';']);
        clear d_str

    end
end