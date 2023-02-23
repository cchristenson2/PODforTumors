
%{
Inputs
    - vec: current local parameter vector, in reduced form
    - Lib: Library of reduced operators for the given parameter
    - k: rank of projection matrix
    - fmt: naming convention for library operators

%}

function [Operator, Modes] = OperatorInterp_local(vec, Lib, k)
    
    
    Modes = struct;
    for j = 1:k %Loop through each local parameter
        curr = vec(j);
        
        eval(['bound = Lib.Mode',num2str(j),'.vec;']);
        bound_low = 1; bound_up = 2;
        
        %Find the outer operators
        for i = 2:(length(bound)-1)
             if(curr>=bound(i))
                 bound_low = i;
                 bound_up = i+1;
             else
                 break;
             end
        end
        
        low = bound(bound_low);
        up  = bound(bound_up);
        
        %Get ratio for interpolation scheme  
        total_dist = abs(up - low);
        low_dist = abs(curr - low);
        up_dist  = abs(up - curr);

        low_ratio = 1 - low_dist/total_dist;
        up_ratio  = 1 - up_dist/total_dist;
        
        %Get operators to use for interp
        %1 - low op
        
        eval(['OP1 = Lib.Mode',num2str(j),'.OP',num2str(bound_low),';']);
        
        %2 - up op
        eval(['OP2 = Lib.Mode',num2str(j),'.OP',num2str(bound_up),';']);
        
        %Interp changing d, low k; 3-4
        OP = low_ratio.*OP1 + up_ratio.*OP2;
        
        eval(['Modes.OP',num2str(j),' = OP;']);
        
    end
    
    Operator = zeros(size(OP));
    for i = 1:k
        eval(['temp = Modes.OP',num2str(i),';']);
        Operator = Operator + temp;
    end
end
