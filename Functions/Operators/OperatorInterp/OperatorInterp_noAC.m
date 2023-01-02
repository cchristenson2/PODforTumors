%{ 
Interpolation of all operators required for no treatment forward model

Inputs:
    - current proliferation
    - current diffusivity
    - parameter bounds
    - Operator libraries (works with FOM or ROM libraries)
        - A, B, H

Outputs:
    - New estimated A, B, and H operators for current parameters

Contributors: Chase Christenson
%}
function [A,B,H] = OperatorInterp_noAC(k,d,bounds,A_lib,B_lib,H_lib)
    kp_vec    = bounds.kp_bounds;
    d_vec     = bounds.d_bounds;
    fmt       = bounds.fmt;
    
    %Get closest k/d pairs in library
    k_low = kp_vec(1); k_up = kp_vec(2);
    for i = 3:length(kp_vec)-1
         if(k>=kp_vec(i))
             k_low = kp_vec(i);
             k_up = kp_vec(i+1);
         else
             break;
         end
    end
    d_low = d_vec(1); d_up = d_vec(2);
    for i = 3:length(d_vec)-1
         if(d>=d_vec(i))
             d_low = d_vec(i);
             d_up = d_vec(i+1);
         else
             break;
         end
    end
    
    %Get ratio for interpolation scheme
    k_dist = k_up - k_low;
    k_low_dist = k - k_low;
    k_up_dist  = k_up - k;
    
    d_dist = d_up - d_low;
    d_low_dist = d - d_low;
    d_up_dist  = d_up - d;
    
    k_low_ratio = 1- k_low_dist/k_dist;
    k_up_ratio  = 1- k_up_dist/k_dist;
    
    d_low_ratio = 1- d_low_dist/d_dist;
    d_up_ratio  = 1- d_up_dist/d_dist;
    
    %Get operators to use for interp
    %1 - low k
    kp_str = replace(num2str(k_low,fmt),'.','_');
    kp_str = replace(kp_str,'-','n');
    name = ['kp',kp_str];
    eval(['B1 = B_lib.',name,'_B;']);
    eval(['H1 = H_lib.',name,'_H;']);
    %2 - up k
    kp_str = replace(num2str(k_up,fmt),'.','_');
    kp_str = replace(kp_str,'-','n');
    name = ['kp',kp_str];
    eval(['B2 = B_lib.',name,'_B;']);
    eval(['H2 = H_lib.',name,'_H;']);
    
    %3 - low d
    d_str  = replace(num2str(d_low,fmt),'.','_');
    name = ['d',d_str];
    eval(['A3 = A_lib.',name,'_A;']);
    %4 - up d
    d_str  = replace(num2str(d_up,fmt),'.','_');
    name = ['d',d_str];
    eval(['A4 = A_lib.',name,'_A;']);
    
    
    %Interp changing k, low d; 1-2
    B = k_low_ratio.*B1 + k_up_ratio.*B2;
    H = k_low_ratio.*H1 + k_up_ratio.*H2;
    
    
    %Interp changing d, low k; 3-4
    A = d_low_ratio.*A3 + d_up_ratio.*A4;
end