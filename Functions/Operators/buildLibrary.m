%{ 
Builds ROM library for a set of parameter bounds, and grid sizes
Operator type determined by input

Inputs:
    - bounds; bounds for the parameters of interest
        n x 2 vector at a minimum; giving upper and lower bounds for n parameters
    - N; Cell map for sizing, 2D or 3D cell map 
    - h, dz; grid spacing, in-plane (h), slice (dz)
    - fmt; name formating
    - params; param type
        - 'A', 'B', 'H', 'T' depending on operator type needed
    - boundary conditions
    - number of operators to build total, including the bounds (2 minimum)

Outputs:
    - lib; library for all operators in the bounds vector, and each parameter

Contributors: Chase Christenson
%}

function [libs] = buildLibrary(bounds,N,h,dz,bcs,params,fmt,num)
    
    n_params = size(params,1);

    libs = struct;
    
    for n = 1:n_params
        param_name = params{n,1};
        bound_loc  = params{n,2};
        
        if(num>=2)
            curr_vec = linspace(bounds(bound_loc,1),bounds(bound_loc,2),num);
        else
            disp('Number of operators to build must be 2 or higher');
            libs = NaN;
            return;
        end

        temp_lib = struct;
        
        for i = 1:numel(curr_vec)
            curr = curr_vec(i);
            curr_str = replace(num2str(curr,fmt),'.','_');
            name = ['val',curr_str];
            
            switch param_name
                case 'A' %Diffusivity, laplacian operator
                    eval(['temp_lib.',name,' = assembleA(N, curr, h, dz, bcs);']);

                case 'B' %Linear proliferation operator
                    eval(['temp_lib.',name,' = assembleB(N, curr);']);

                case 'H' %Quadratic proliferation operator
                    eval(['temp_lib.',name,' = assembleH(N, curr);']);

                case 'T' %Linear treatment efficacy operator
                    eval(['temp_lib.',name,' = assembleT(N, curr);']);
                otherwise
                    disp('Param type not accepted for build; ''A'', ''B'', ''H'', or ''T''');
                    libs = NaN;
                    return;
            end
        end
        eval(['libs.',param_name,'_lib = temp_lib;']);
    end
    
end