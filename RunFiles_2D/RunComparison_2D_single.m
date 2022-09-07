%{ 
Script to run calibration comparison in 2D with global proliferation
Calibrates for a single hardcoded patient for testing
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


name = files(3).name;
savepath = [home,'/Results/nonpcr/',erase(name,'.mat'),'_2D_calibration'];
file_loc = [location,name];
[FOM,ROM, tumor] = calibrationComparison_local_2D(file_loc, cal_flag, bounds, savepath);
% save([home,'/Results/nonpcr/',erase(name,'.mat'),'_LocalCalibration_2D'], 'FOM', 'ROM', 'tumor');
