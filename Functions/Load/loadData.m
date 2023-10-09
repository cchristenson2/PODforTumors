%{ 
Loading patient data from file

Inputs:
    - location; File location

Outputs:
    - Tumor struct
        - N; 3D (2D grids) or 4D (3D grids) tumor cell map
        - AUC; concentration map for drug delivery
        - Mask; breast mask for boundary conditions
        - Tissues; segmentation for different tissue types
        - t_scan; days away from 0 for each scan
        - t_trx; days away from 0 for each treatment
        - h; in-plane grid spacing
        - dz; slice grid spacing
        - params; true kp, d, and alpha that made the tumor
            - only exists for virtually created data
        - beta1, beta2; true decay values for A/C chemo
            - only exists for virtually created data

Contributors: Chase Christenson
%}

function [tumor] = loadData(location)
    temp = load(location);
    
    %Stack cell counts by time point in 4th dimension
    names = fields(temp.image_data);
    idx = contains(names, 'NTC');
    names(~idx) = [];
    names = sort(names);
    N = zeros([size(temp.image_data.NTC1), numel(names)]);
    if(size(temp.image_data.NTC1,3)==1) %2D grid
        for i = 1:numel(names)
            eval(['N(:,:,i) = temp.image_data.NTC',num2str(i),';']);
        end
        dims = 2;
    else %3D grid
        for i = 1:numel(names)
            eval(['N(:,:,:,i) = temp.image_data.NTC',num2str(i),';']);
        end
        dims = 3;
    end
    
    theta = 0.7405 * ((temp.schedule_info.imagedims(1))^2 * temp.schedule_info.imagedims(3)) / 4.189e-6;
    
    N = N./theta;
    
    if(dims==2)
        tumor.bcs = buildBoundaries_2D(logical(temp.image_data.BreastMask));
    else
        tumor.bcs = buildBoundaries_3D(logical(temp.image_data.BreastMask));
    end
    
    %Get times from baseline for scheduling
    t_all = cumsum(temp.schedule_info.times);
    t_scan = t_all(temp.schedule_info.schedule=='S');
    t_trx  = t_all(temp.schedule_info.schedule=='A');
    
    tumor.N = N;
    tumor.AUC = temp.image_data.AUC;
    tumor.Mask = temp.image_data.BreastMask;
    tumor.Tissues = temp.image_data.Tissues;
    tumor.t_scan = t_scan;
    tumor.t_trx = t_trx;
    tumor.h = temp.schedule_info.imagedims(1);
    tumor.dz = temp.schedule_info.imagedims(3);
    
    try
        tumor.params.kp = temp.params.kp;
    catch 
    end
    try
        tumor.params.d = temp.params.d;
    catch
    end
    try
        tumor.params.alpha = temp.params.alpha;
    catch
    end
    try
        tumor.beta1 = temp.params.beta1;
    catch
    end
    try
        tumor.beta2 = temp.params.beta2;
    catch
    end
    
end