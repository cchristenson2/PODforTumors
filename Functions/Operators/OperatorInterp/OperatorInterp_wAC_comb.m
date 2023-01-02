%{ 
Interpolation of all operators required for A/C treatment forward model

Inputs:
    - current proliferation
    - current diffusivity
    - current A efficacy
    - current C efficacy
    - parameter bounds
    - Operator libraries (works with FOM or ROM libraries)
        - A, B, H, T

Outputs:
    - New estimated A, B, H, Ta, and Tc operators for current parameters

Contributors: Chase Christenson
%}
function [A,B,H,Ta] = OperatorInterp_wAC_comb(k,d,alpha1,bounds,A_lib,B_lib,H_lib,T_lib)
    kp_vec    = bounds.kp_bounds;
    d_vec     = bounds.d_bounds;
    alpha_vec = bounds.alpha_bounds;
    fmt       = bounds.fmt;
    
    %Get closest k
    k_low = kp_vec(1); k_up = kp_vec(2);
    for i = 3:length(kp_vec)-1
         if(k>=kp_vec(i))
             k_low = kp_vec(i);
             k_up = kp_vec(i+1);
         else
             break;
         end
    end
    %Get closest d
    d_low = d_vec(1); d_up = d_vec(2);
    for i = 3:length(d_vec)-1
         if(d>=d_vec(i))
             d_low = d_vec(i);
             d_up = d_vec(i+1);
         else
             break;
         end
    end
    %Get closest alpha_a
    alpha1_low = alpha_vec(1); alpha1_up = alpha_vec(2);
    for i = 3:length(alpha_vec)-1
         if(alpha1>=alpha_vec(i))
             alpha1_low = alpha_vec(i);
             alpha1_up = alpha_vec(i+1);
         else
             break;
         end
    end
%     %Get closest alpha_c
%     alpha2_low = alpha_vec(1); alpha2_up = alpha_vec(2);
%     for i = 3:length(alpha_vec)-1
%          if(alpha2>=alpha_vec(i))
%              alpha2_low = alpha_vec(i);
%              alpha2_up = alpha_vec(i+1);
%          else
%              break;
%          end
%     end
    
    %Get ratio for interpolation scheme
    k_dist = k_up - k_low;
    k_low_dist = k - k_low;
    k_up_dist  = k_up - k;
    
    d_dist = d_up - d_low;
    d_low_dist = d - d_low;
    d_up_dist  = d_up - d;
    
    alpha1_dist = alpha1_up - alpha1_low;
    alpha1_low_dist = alpha1 - alpha1_low;
    alpha1_up_dist  = alpha1_up - alpha1;
    
%     alpha2_dist = alpha2_up - alpha2_low;
%     alpha2_low_dist = alpha2 - alpha2_low;
%     alpha2_up_dist  = alpha2_up - alpha2;
    
    k_low_ratio = 1- k_low_dist/k_dist;
    k_up_ratio  = 1- k_up_dist/k_dist;
    
    d_low_ratio = 1- d_low_dist/d_dist;
    d_up_ratio  = 1- d_up_dist/d_dist;
    
    alpha1_low_ratio = 1- alpha1_low_dist/alpha1_dist;
    alpha1_up_ratio  = 1- alpha1_up_dist/alpha1_dist;
    
%     alpha2_low_ratio = 1- alpha2_low_dist/alpha2_dist;
%     alpha2_up_ratio  = 1- alpha2_up_dist/alpha2_dist;
    
    %Get operators to use for interp
    %1 - low k
    kp_str = replace(num2str(k_low,fmt),'.','_');
    name = ['kp',kp_str];
    eval(['B1 = B_lib.',name,'_B;']);
    eval(['H1 = H_lib.',name,'_H;']);
    %2 - up k
    kp_str = replace(num2str(k_up,fmt),'.','_');
    name = ['kp',kp_str];
    eval(['B2 = B_lib.',name,'_B;']);
    eval(['H2 = H_lib.',name,'_H;']);
    
    %3 - low d
    d_str  = replace(num2str(d_low,fmt),'.','_');
    name = ['d',d_str];
    eval(['A1 = A_lib.',name,'_A;']);
    %4 - up d
    d_str  = replace(num2str(d_up,fmt),'.','_');
    name = ['d',d_str];
    eval(['A2 = A_lib.',name,'_A;']);
    
    %5 - low alpha1
    alpha_str  = replace(num2str(alpha1_low,fmt),'.','_');
    name = ['alpha',alpha_str];
    eval(['Ta1 = T_lib.',name,'_T;']);
    %6 - up alpha1
    alpha_str  = replace(num2str(alpha1_up,fmt),'.','_');
    name = ['alpha',alpha_str];
    eval(['Ta2 = T_lib.',name,'_T;']);
    
%     %7 - low alpha2
%     alpha_str  = replace(num2str(alpha2_low,fmt),'.','_');
%     name = ['alpha',alpha_str];
%     eval(['Tc1 = T_lib.',name,'_T;']);
%     %8 - up alpha2
%     alpha_str  = replace(num2str(alpha2_up,fmt),'.','_');
%     name = ['alpha',alpha_str];
%     eval(['Tc2 = T_lib.',name,'_T;']);
    
    %Interp changing k, low d; 1-2
    B = k_low_ratio.*B1 + k_up_ratio.*B2;
    H = k_low_ratio.*H1 + k_up_ratio.*H2;

    %Interp changing d, low k; 3-4
    A = d_low_ratio.*A1 + d_up_ratio.*A2;
    
    %Interp changing alphas
    Ta = alpha1_low_ratio.*Ta1 + alpha1_up_ratio.*Ta2;
%     Tc = alpha2_low_ratio.*Tc1 + alpha2_up_ratio.*Tc2;
end