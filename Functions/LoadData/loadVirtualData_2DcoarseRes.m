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

function [tumor] = loadVirtualData_coarseRes(location)
    data = load(location);
    
    tumor = data.tumor;
end