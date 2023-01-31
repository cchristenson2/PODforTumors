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

function [Ar_lib, Br_lib, Hr_lib, Tr_lib] = reduceOperatorLibrary(A_lib, B_lib, H_lib, T_lib, V, bounds, tx_params)
    fmt = bounds.fmt;
    num_kp = numel(bounds.kp_bounds);
    num_d = numel(bounds.d_bounds);
    num_a = numel(bounds.alpha_bounds);
    
    Ar_lib = struct;
    Br_lib = struct;
    Hr_lib = struct;
    Tr_lib = struct;
    
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
    
    for i = 1:num_d
        d      = bounds.d_bounds(i);
        d_str     = replace(num2str(d,fmt),'.','_');
        name = ['d',d_str];
        
        eval(['Ar_lib.',name,'_A = V'' * A_lib.',name,'_A * V;']);
        eval(['clear ',name,';']);
        clear d_str

        
        
    end
    
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