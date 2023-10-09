%{ 
Interpolation of A operator only

Inputs:
    - val; current parameter value
    - lib; operator library (works with FOM or ROM libraries)

Outputs:
    - OP; New estimated operator for current parameter value

Contributors: Chase Christenson
%}
function [OP] = OperatorInterp(val,lib)
    %Get bounds vector
    fields = fieldnames(lib);
    for i = 1:numel(fields)
        vec(i) = str2double(replace(erase(fields{i},'val'),'_','.'));
    end

    if(numel(vec)>2)
        low = vec(1); up = vec(2);
        low_loc = 1; up_loc  = 2;
        for i = 2:length(vec)-1
             if(val>=vec(i))
                 low = vec(i);
                 up  = vec(i+1);
                 low_loc = i; up_loc  = i+1;
             else
                 break;
             end
        end
    else
        low = vec(1);
        up  = vec(2);
        
        low_loc = 1;
        up_loc  = 2;
    end

    %Get ratio for interpolation scheme  
    dist = up - low;
    low_dist = val - low;
    up_dist  = up - val;
    
    low_ratio = 1 - low_dist/dist;
    up_ratio  = 1 - up_dist/dist;
    
    %Get operators to use for interp
    low_str = fields{low_loc};
    up_str  = fields{up_loc};
    
    eval(['OP1 = lib.',low_str,';']);
    eval(['OP2 = lib.',up_str,';']);
    
    %Interp between the library
    OP = low_ratio.*OP1 + up_ratio.*OP2;
end