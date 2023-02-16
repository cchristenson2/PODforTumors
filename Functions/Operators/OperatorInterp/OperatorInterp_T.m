%{ 
Interpolation of Ta and Tc operator only

Inputs:
    - current efficacies for both treatments
    - parameter bounds
    - T operator library (works with FOM or ROM libraries)

Outputs:
    - New estimated Ta and Tc operatos for current efficacies

Contributors: Chase Christenson
%}
function [Ta, Tc] = OperatorInterp_T(alpha1,alpha2,bounds,T_lib)
    alpha_vec = bounds.alpha_bounds;
    fmt       = bounds.fmt;
    
    %Get closest alpha_a
    alpha1_low = alpha_vec(1); alpha1_up = alpha_vec(2);
    for i = 2:length(alpha_vec)-1
         if(alpha1>=alpha_vec(i))
             alpha1_low = alpha_vec(i);
             alpha1_up = alpha_vec(i+1);
         else
             break;
         end
    end
    %Get closest alpha_c
    alpha2_low = alpha_vec(1); alpha2_up = alpha_vec(2);
    for i = 2:length(alpha_vec)-1
         if(alpha2>=alpha_vec(i))
             alpha2_low = alpha_vec(i);
             alpha2_up = alpha_vec(i+1);
         else
             break;
         end
    end
    
    alpha1_dist = alpha1_up - alpha1_low;
    alpha1_low_dist = alpha1 - alpha1_low;
    alpha1_up_dist  = alpha1_up - alpha1;
    
    alpha2_dist = alpha2_up - alpha2_low;
    alpha2_low_dist = alpha2 - alpha2_low;
    alpha2_up_dist  = alpha2_up - alpha2;
    
    alpha1_low_ratio = 1- alpha1_low_dist/alpha1_dist;
    alpha1_up_ratio  = 1- alpha1_up_dist/alpha1_dist;
    
    alpha2_low_ratio = 1- alpha2_low_dist/alpha2_dist;
    alpha2_up_ratio  = 1- alpha2_up_dist/alpha2_dist;
    
    %1 - low alpha1
    alpha_str  = replace(num2str(alpha1_low,fmt),'.','_');
    name = ['alpha',alpha_str];
    eval(['Ta1 = T_lib.',name,'_T;']);
    %2 - up alpha1
    alpha_str  = replace(num2str(alpha1_up,fmt),'.','_');
    name = ['alpha',alpha_str];
    eval(['Ta2 = T_lib.',name,'_T;']);
    
    %1 - low alpha2
    alpha_str  = replace(num2str(alpha2_low,fmt),'.','_');
    name = ['alpha',alpha_str];
    eval(['Tc1 = T_lib.',name,'_T;']);
    %2 - up alpha2
    alpha_str  = replace(num2str(alpha2_up,fmt),'.','_');
    name = ['alpha',alpha_str];
    eval(['Tc2 = T_lib.',name,'_T;']);
    
    %Interp changing alphas
    Ta = alpha1_low_ratio.*Ta1 + alpha1_up_ratio.*Ta2;
    Tc = alpha2_low_ratio.*Tc1 + alpha2_up_ratio.*Tc2;
    
end