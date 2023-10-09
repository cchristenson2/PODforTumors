%{ 
Reduce ROM library for a set of parameter bounds, and given projection matrix
Operator type determined by input

Inputs:
    - lib; operator libraries
    - V; projection matrix
    - bounds; parameter bounds
    - param; param type
        - 'A', 'B', 'H', 'T' depending on operator type needed
    - C; concentration map to locally scale treatment efficacy, only applied if param = 'T'

Outputs:
    - r_lib; reduced library for all operators in the bound vector

Contributors: Chase Christenson
%}

function [r_libs] = reduceLibrary(libs, V, params, C)
    
    n_params = size(params,1);

    r_libs = struct;
    
    for n = 1:n_params
        param_name = params{n,1};
        eval(['lib_f = libs.',param_name,'_lib;']);
        
        fields = fieldnames(lib_f);
        
        temp_lib = struct;
    
        for i = 1:numel(fields)
            %Get current string from lib name
            name = fields{i};
            
            switch param_name
                case 'A' %Diffusivity, laplacian operator
                    eval(['temp_lib.',name,' = V'' * lib_f.',name,' * V;']);

                case 'B' %Linear proliferation operator
                    eval(['temp_lib.',name,' = V'' * lib_f.',name,' * V;']);

                case 'H' %Quadratic proliferation operator
                   eval(['temp_lib.',name,' = V'' * lib_f.',name,' * kron(V,V);']);

                case 'T' %Linear treatment efficacy operator                   
                    eval(['temp_lib.',name,' = V'' * (lib_f.',name,' .* diag(C(:))) * V;']);
                otherwise
                    disp('Param type not accepted for reduction; ''A'', ''B'', ''H'', or ''T''');
                    r_libs = NaN;
                    return;
            end
        end
        
        eval(['r_libs.',param_name,'r_lib = temp_lib;']);
        
    end
    
end