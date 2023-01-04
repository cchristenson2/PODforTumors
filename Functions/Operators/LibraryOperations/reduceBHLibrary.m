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

function [Br_lib, Hr_lib] = reduceBHLibrary(B_lib, H_lib, V, bounds)
    fmt = bounds.fmt;
    num_kp = numel(bounds.kp_bounds);
    num_d = numel(bounds.d_bounds);
    
    Br_lib = struct;
    Hr_lib = struct;
    
    for i = 1:num_kp
        kp     = bounds.kp_bounds(i);
        kp_str    = replace(num2str(kp,fmt),'.','_');
        kp_str    = replace(kp_str,'-','n');

        %kp
        name = ['kp',kp_str];
        
%         eval(['disp(size(B_lib.',name,'_B));']);
%         eval(['disp(size(H_lib.',name,'_H));']);
        
        
        eval(['Br_lib.',name,'_B = V'' * B_lib.',name,'_B * V;']);
        eval(['Hr_lib.',name,'_H = V'' * H_lib.',name,'_H * kron(V,V);']);
        eval(['clear ',name,';']);
        clear kp_str
    end

end