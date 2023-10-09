%{ 
Script for generating example virtual data
This script is for a base diffusion + logistic growth equation.
Data is coded to match patient data formatting from processed MRI

Saves:
    - image_data struct
        - NTCi; cell number maps where i is the scan number (i.e., 1 for baseline image)
        - Tissues; segemented map of the different tissue types in the image
        - BreastMask; segmented breast mask
        - AUC; area under the curve used to estimate chemo concentrations
            - zeros if drug was not included in the model
    - schedule_info struct
        - schedule; vector of the different event types a patient has
            - 'A' for chemo delivery, 'S' for imaging visit
        - times; time since the last event in schedule
        - imagedims; [y, x, and z] resolutions in vector format
    - params struct
        - kp; true proliferation rate
        - d; true diffusivity
        - alpha; true drug efficacy

Contributors: Chase Christenson
%}

addpath(genpath(fileparts(pwd)))


kp_range = [1e-3,1e-1];
d_range = [1e-6,1e-1];

%Set grid spacing
h  = 1; %mm
dz = 2; %mm
schedule_info.imagedims = [h, h, dz];

theta = 0.7405 * ((schedule_info.imagedims(1))^2 * schedule_info.imagedims(3)) / 4.189e-6; %Carrying capacity for the voxel size

%Set 100 X 100 grid for cell data
sy = 100;
sx = 100;
N0 = zeros(sy,sx);

%These matrices would typically be provided by imaging data
%Here a simple alternative is provided
Tissues       = ones(sy,sx);
BreastMask    = ones(sy,sx);
AUC           = zeros(sy,sx);

%Set timing of events for imaging vs treatment
schedule_info.schedule = ['S'; 'S'; 'S']; %S for scan, A for chemo

days = [0, 30, 60];
schedule_info.times = [0, diff(days)]; %Days since last event array

%Seed initial tumor
[X,Y] = meshgrid(1:h:sx,1:h:sy);
z = exp(-((X-(sx-1)/2).^2+(Y-(sy-1)/2).^2)/5);
N0(z>0.10) = 0.50; %Initial volume fraction of 50%

%Sample true parameters for the virtual data
kp = kp_range(1) + rand(1,1).*(diff(kp_range));
d  = d_range(1) + rand(1,1).*(diff(d_range));

params.kp = kp;
params.d  = d;

%Run forward evaluation to get true data
bcs = buildBoundaries_2D(BreastMask);
t = cumsum(schedule_info.times(schedule_info.schedule=='S'));

dt = 0.5;

N_true = RXDIF_2D(N0, kp.*ones(sy,sx), d, t, h, dt, bcs);

%Prep variables and store

N_true = N_true.*theta;
for i = 1:size(N_true,3)
    eval(['image_data.NTC',num2str(i),' = N_true(:,:,i);']);
end

image_data.Tissues    = Tissues;
image_data.BreastMask = BreastMask;
image_data.AUC        = AUC;

save('Ex1_patient.mat','image_data','params','schedule_info');
