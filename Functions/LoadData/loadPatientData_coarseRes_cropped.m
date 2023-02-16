%{ 
Loading patient data from file, corase resolution

Inputs:
    - Patient number/name
    - File location

Outputs:
    - Tumor struct
        - N; 4D tumor cell map
        - AUC; concentration map for drug delivery
        - Mask; breast mask for boundary conditions
        - t_scan; days away from 0 for each scan
        - t_trx; days away from 0 for each treatment
        - h; x,y grid spacing
        - dz; z grid spacing
        - CF; correction for downsampling

Contributors: Chase Christenson
%}

function [tumor] = loadPatientData_coarseRes_cropped(location)
    data = load(location);
    
    %Stack cell counts by time point in 4th dimension
    names = fields(data.coarse_res_dat_cropped);
    idx = contains(names, 'NTC');
    names(~idx) = [];
    names = sort(names);
    N = zeros([size(data.coarse_res_dat_cropped.NTC1), numel(names)]);
    for i = 1:numel(names)
        eval(['N(:,:,:,i) = data.coarse_res_dat_cropped.NTC',num2str(i),';']);
    end
    N = N./max(N(:));
    
    %Get times from baseline for scheduling
    t_all = cumsum(data.schedule_info.times);
    t_scan = t_all(data.schedule_info.schedule=='S');
    t_trx  = t_all(data.schedule_info.schedule=='A');
    
    tumor.N = N;
    tumor.AUC = data.coarse_res_dat_cropped.AUC;
    tumor.Mask = data.coarse_res_dat_cropped.BreastMask;
    tumor.Tissues = data.coarse_res_dat_cropped.Tissues;
    tumor.t_scan = t_scan;
    tumor.t_trx = t_trx;
    tumor.h = data.schedule_info.imagedims(1)*data.coarse_res_dat.CF;
    tumor.dz = data.schedule_info.imagedims(3);
    tumor.CF = data.coarse_res_dat_cropped.CF;
end