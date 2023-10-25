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

function [lib,Mode_bounds] = buildAndReduce_LocalKp(bounds,V,fmt,num,param)

    lib = struct;
    [n,k] = size(V);
    
    %Estimate the bounds for the modal coefficients
    Mode_bounds = zeros(k,2);
    for i = 1:k
        %Lower bound
        kp_test = zeros(n,1);
        kp_test(V(:,i)>=0) = bounds(1);
        kp_test(V(:,i)<0) = bounds(2);
        Mode_bounds(i,1) = (V(:,i)' * kp_test)*2;
        
        %Upper bound
        kp_test = zeros(n,1);
        kp_test(V(:,i)>=0) = bounds(2);
        kp_test(V(:,i)<0) = bounds(1);
        Mode_bounds(i,2) = (V(:,i)' * kp_test)*2;
    end
    
    %Build library for each mode
    for i = 1:k
        if(num>=2)
            curr_vec = linspace(Mode_bounds(i,1),Mode_bounds(i,2),num);
        else
            disp('Number of operators to build must be 2 or higher');
            lib = NaN;
            return;
        end
        
        %Build operator for each potential modal coefficient
        for j = 1:num
            curr = curr_vec(j);
            temp_map = V(:,i) * curr;
            switch param
                case 'B'
                    temp_op = zeros(k,k);
                    for z = 1:n
                        temp_op = temp_op + V(z,:)' * temp_map(z) * V(z,:);
                    end
                case 'H'
                    temp_op = zeros(k,k^2);
                    for z = 1:n
                        temp_op = temp_op + V(z,:)' * temp_map(z) * kron(V(z,:),V(z,:));
                    end
                otherwise
                    disp('Param type not accepted for build; ''B'', or ''H''');
                    lib = NaN;
                    return;

            end

            name = ['val',replace(replace(num2str(curr,fmt),'.','_'), '-','neg')];

            eval(['lib.Mode',num2str(i),'.',name,' = temp_op;']);
        end
    end
    
end