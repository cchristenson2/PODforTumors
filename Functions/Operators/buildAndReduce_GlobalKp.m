%{ 
Builds ROM library for a set of parameter bounds, and grid sizes
Operator type determined by input

Inputs:
    - bounds; bounds for the parameters of interest
        n x 2 vector at a minimum; giving upper and lower bounds for n parameters
    - N; Cell map for sizing, 2D or 3D cell map 
    - fmt; name formating
    - number of operators to build total, including the bounds (2 minimum)

Outputs:
    - lib; library for all operators in the bounds vector, and each parameter

Contributors: Chase Christenson
%}

function [lib] = buildAndReduce_GlobalKp(bounds,V,fmt,num,param)

    lib = struct;
    
    [n,k] = size(V);
    
    if(num>=2)
        curr_vec = linspace(bounds(1),bounds(2),num);
    else
        disp('Number of operators to build must be 2 or higher');
        lib = NaN;
        return;
    end
    
    for i = 1:numel(curr_vec)
        
        curr = curr_vec(i);
        switch param
            case 'B'
                temp_op = zeros(k,k);
                for j = 1:n
                    temp_op = temp_op + V(j,:)' * curr * V(j,:);
                end
            case 'H'
                temp_op = zeros(k,k^2);
                for j = 1:n
                    temp_op = temp_op + V(j,:)' * curr * kron(V(j,:),V(j,:));
                end
            otherwise
                disp('Param type not accepted for build; ''B'', or ''H''');
                lib = NaN;
                return;
                
        end
        
        name = ['val',replace(num2str(curr,fmt),'.','_')];
        
        
        eval(['lib.',name,' = temp_op;']);
    end
    
end