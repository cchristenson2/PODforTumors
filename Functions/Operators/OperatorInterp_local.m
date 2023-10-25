%{ 
Interpolation of A operator only

Inputs:
    - val; current parameter value
    - lib; operator library (works with FOM or ROM libraries)

Outputs:
    - OP; New estimated operator for current parameter value

Contributors: Chase Christenson
%}
function [Operator] = OperatorInterp_local(vals,lib)
    %Get number of modes
    fields = fieldnames(lib);
    k = numel(fields);
    
    Modes = struct;
    %Interpolate for each mode
    for i = 1:k
        curr = vals(i);
        %Get bounds vector
        eval(['mode_fields = fieldnames(lib.',fields{i},');']);
        for j = 1:numel(mode_fields)
            vec(j) = str2double(replace(replace(erase(mode_fields{j},'val'),'_','.'),'neg','-'));
        end
        %Search for location in library
        if(numel(vec)>2)
            low = vec(1); up = vec(2);
            low_loc = 1; up_loc  = 2;
            for j = 2:length(vec)-1
                 if(curr>=vec(j))
                     low = vec(j);
                     up  = vec(j+1);
                     low_loc = j; up_loc  = j+1;
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
        low_dist = curr - low;
        up_dist  = up - curr;

        low_ratio = 1 - low_dist/dist;
        up_ratio  = 1 - up_dist/dist;
        
        %Get operators to use for interp
        low_str = mode_fields{low_loc};
        up_str  = mode_fields{up_loc};
        
        eval(['OP1 = lib.',fields{i},'.',low_str,';']);
        eval(['OP2 = lib.',fields{i},'.',up_str,';']);
        
        OP = low_ratio.*OP1 + up_ratio.*OP2;
        
        eval(['Modes.OP',num2str(i),' = OP;']);
        
        clear vec;
    end
    
    %Add contribution from each mode
    Operator = zeros(size(OP));
    for i = 1:k
        eval(['temp = Modes.OP',num2str(i),';']);
        Operator = Operator + temp;
    end
end