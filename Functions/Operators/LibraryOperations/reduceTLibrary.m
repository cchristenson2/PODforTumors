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

function [Tr_lib] = reduceTLibrary(T_lib, V, bounds, tx_params)
    fmt = bounds.fmt;
    num_a = numel(bounds.alpha_bounds);

    Tr_lib = struct;
    
    for i = 1:num_a
        alpha  = bounds.alpha_bounds(i);
        alpha_str = replace(num2str(alpha,fmt),'.','_');

        %alpha
        name = ['alpha',alpha_str];
        eval(['Tr_lib.',name,'_T = V''*(T_lib.',name,'_T.*diag(tx_params.C(:)))*V;']); %Include patient specific concentration map before reducing
        eval(['clear ',name,';']);
        clear alpha_str
    end
end