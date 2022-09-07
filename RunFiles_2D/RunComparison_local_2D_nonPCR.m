%{ 
Script to run calibration comparison in 2D with local proliferation
Calibrates for all patients in hardcoded folder
%}
home = fileparts(pwd);
addpath(genpath(home))

bounds.d_bounds = logspace(-6, -3, 10);
bounds.kp_bounds = logspace(-3, -1, 10);
bounds.alpha_bounds = logspace(-6, 0, 10);
bounds.fmt = '%.6f';

cal_flag = 0;

%Non-PCR Patients
location = [home,'/Data/PatientData/nonpcr/'];
files = dir(location);
for i = 3:numel(files)
    name = files(i).name;
    file_loc = [location,name];
    savepath = [home,'/Results/nonpcr/Figs/',erase(name,'.mat'),'_2D_calibration'];
    
    [FOM,ROM, tumor] = calibrationComparison_local_2D(file_loc, cal_flag, bounds, savepath);
    save([home,'/Results/nonpcr/',erase(name,'.mat'),'_LocalCalibration_2D'], 'FOM', 'ROM', 'tumor');
    disp(['Non-PCR Patient ',num2str(i),' complete']);
end