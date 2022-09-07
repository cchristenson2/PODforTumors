%{ 
Calibrate tumor for a specific patient using both FOM & ROM
Proliferation is local, all other parameters are global
Uses A/C treatment forward model

Inputs:
    - location; file location
    - Cal_flag; 0 for full calibration, 1 for prediction at last time point
    - Parameter bounds
    - savepath for figure outputs

Outputs:
    - FOM
        - Calibrated params
        - Statistics vs. measured data
        - Final simulation and prediction (if needed)
        - time to calibrate
    - ROM
        - Calibrated params
        - Statistics vs. measured data
        - Final simulation and prediction (if needed)
        - time to calibrate
    - tumor used for calibration
 
Contributors: Chase Christenson
%}

function [FOM, ROM, tumor] = calibrationComparison_local_2D(location, cal_flag, bounds, savepath)

%     %Add ROM functions to path
%     addpath(genpath('C:/ROM/'))

    %% Load and Prep Data
    tumor = loadPatientData_coarseRes(location);
    
    %Select middle slice for 2D data
    [~,~,sz,ntp] = size(tumor.N);
    slice = round(sz/2);

    tumor.N    = squeeze(tumor.N(:,:,slice,:));
    tumor.AUC  = tumor.AUC(:,:,slice);
    tumor.Mask = tumor.Mask(:,:,slice);
    tumor.bcs  = buildBoundaries_2D(tumor.Mask);
    
    %Set beta within range for trx response, could switch to a random sampling method
    tumor.beta1 = 0.25;
    tumor.beta2 = 2;

    %% Offline - prebuild library
    start = tic;
    num   = numel(tumor.N(:,:,1));
    [A_lib, B_lib, H_lib, T_lib] = buildOperatorLibrary(bounds, num, tumor.h, bounds.fmt, tumor.bcs);
    time_Offline = toc(start);

    %% Calibrations
    %Set number of calibrated time points
    if(cal_flag==0)
        ntp_cal  = ntp - 1;
        ntp_pred = 0;
    else
        ntp_cal  = ntp - 2;
        ntp_pred = 1;
    end
    %Get ROI for calibrated data
    tumor.ROI = logical(buildROI(tumor, ntp_cal));
    
    dt = 0.50; %Time spacing [days]
    
    %Full order model calibration
    start = tic;
    [params_FOM, stats_FOM, outputs_FOM, fig_FOM] = FOM_LMCalibration_LocalKp_2D(tumor, ntp_cal, ntp_pred, bounds, dt);
    time_FOM = toc(start); %Seconds
    
    %Reduced order model calibration
    start = tic;
    [params_ROM, stats_ROM, outputs_ROM, fig_ROM, testing] = ROM_LMCalibration_LocalKp_2D(tumor, ntp_cal, ntp_pred, bounds, 5, A_lib, B_lib, H_lib, T_lib, dt);
    time_ROM = toc(start); %Seconds
    
    %% Compare calibrations
    FOM.params  = params_FOM;
    FOM.stats   = stats_FOM;
    FOM.outputs = outputs_FOM;
    FOM.ttc     = time_FOM;
    
    ROM.params  = params_ROM;
    ROM.stats   = stats_ROM;
    ROM.outputs = outputs_ROM;
    ROM.ttc     = time_ROM;
    ROM.offline_time = time_Offline;
    ROM.testing = testing;
    
    savefig(fig_FOM, [savepath, '_FOM.fig']);
    savefig(fig_ROM, [savepath, '_ROM.fig']);
    
end