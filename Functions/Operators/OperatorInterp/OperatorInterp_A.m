%{ 
Interpolation of A operator only

Inputs:
    - current diffusivity
    - parameter bounds
    - A operator library (works with FOM or ROM libraries)

Outputs:
    - New estimated A operator for current diffusivity

Contributors: Chase Christenson
%}
function [A] = OperatorInterp_A(d,bounds,A_lib)
    d_vec     = bounds.d_bounds;
    fmt       = bounds.fmt;
    
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
    d_dist = d_up - d_low;
    d_low_dist = d - d_low;
    d_up_dist  = d_up - d;
    
    d_low_ratio = 1- d_low_dist/d_dist;
    d_up_ratio  = 1- d_up_dist/d_dist;
    
    %Get operators to use for interp
    %1 - low d
    d_str  = replace(num2str(d_low,fmt),'.','_');
    name = ['d',d_str];
    eval(['A1 = A_lib.',name,'_A;']);
    %2 - up d
    d_str  = replace(num2str(d_up,fmt),'.','_');
    name = ['d',d_str];
    eval(['A2 = A_lib.',name,'_A;']);
    
    %Interp changing d, low k; 3-4
    A = d_low_ratio.*A1 + d_up_ratio.*A1;
end