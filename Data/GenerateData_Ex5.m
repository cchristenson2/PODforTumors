%{ 
Script for generating example virtual data
This script is for a base diffusion + logistic growth equation.
Data is coded to match patient data formatting from processed MRI

This generated proliferation map is locally defined rather than a global variable

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

Random field generation from:
Paul Constantine (2023). Random Field Simulation 
(https://www.mathworks.com/matlabcentral/fileexchange/27613-random-field-simulation),
MATLAB Central File Exchange. Retrieved October 24, 2023. 

%}

addpath(genpath(fileparts(pwd)))

kp_range = [1e-3,1e-1];
d_range = [1e-6,1e-2];
alpha_range = [1e-6, 0.8];

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
AUC           = ones(sy,sx);

%Set timing of events for imaging vs treatment
schedule_info.schedule = ['S'; 'A'; 'A'; 'S'; 'A'; 'A';'S']; %S for scan, A for chemo

days = [0, 10, 20, 30, 40, 50, 60];
schedule_info.times = [0, diff(days)]; %Days since last event array

%Seed initial tumor
[X,Y] = meshgrid(1:h:sx,1:h:sy);
z = exp(-((X-(sx-1)/2).^2+(Y-(sy-1)/2).^2)/150);
N0(z>0.10) = 0.50; %Initial volume fraction of 50%

bcs = buildBoundaries_2D(BreastMask);
t_all = cumsum(schedule_info.times);

t = t_all(schedule_info.schedule=='S');

tx_params.txduration = t_all(schedule_info.schedule=='A');
tx_params.beta1 = 0.25;
tx_params.beta2 = 2.0;
tx_params.C     = AUC;

params.beta1 = tx_params.beta1;
params.beta2 = tx_params.beta2;

dt = 0.5;

% corr.name = 'turbulent';
% corr.c0 = 20; corr.c1 = 5; corr.sigma = 1;
% mesh = [X(:) Y(:)]; % 2-D mesh
% temp = randomfield(corr,mesh);
% kp_mean = kp_range(1) + rand(1,1).*(diff(kp_range));
%     
% kp_map = reshape(normalize(temp,'range',kp_range),sy,sx);

[X_r, Y_r] = meshgrid(linspace(1,sx,round(sx/4)),linspace(1,sy,round(sy/4)));
mask = zeros(sy,sx);
mask(N0~=0) = 1;
mask_r = logical(interp2(X,Y,mask,X_r,Y_r));
kp_r = zeros(size(mask_r));
kp_r(mask_r~=0) = kp_range(1) + rand(numel(mask_r(mask_r~=0)),1)*diff(kp_range);
kp_map = interp2(X_r,Y_r,kp_r,X,Y);
    
run = 1;
while run == 1
    %Sample true parameters for the virtual data
    
    d  = d_range(1) + rand(1,1).*(diff(d_range));
    alpha  = alpha_range(1) + rand(1,1).*(diff(alpha_range));

    %Run forward evaluation to get true data
    N_true = RXDIF_2D_wAC_comb(N0, kp_map, d, alpha, tx_params, t, h, dt, bcs);
    
    %Ensure there is residual tumor at last visit
    if sum(N_true(:,:,end),'all') > sum(N0, 'all')*0.05
        run = 0;
    end
end

N_true(N_true<0.01) = 0;
kp_map(sum(N_true,3)==0) = 0;
%Prep variables and store
params.kp = kp_map;
params.d  = d;
params.alpha = alpha;

figure
for i = 1:size(N_true,3)
    subplot(2,size(N_true,3),i)
    imagesc(N_true(:,:,i),[0,1]); axis image; title(['Visit ',num2str(i)]);
end
subplot(2,size(N_true,3),1+size(N_true,3))
imagesc(kp_map); axis image; colorbar;

N_true = N_true.*theta;
for i = 1:size(N_true,3)
    eval(['image_data.NTC',num2str(i),' = N_true(:,:,i);']);
end

image_data.Tissues    = Tissues;
image_data.BreastMask = BreastMask;
image_data.AUC        = AUC;

save('Ex5_patient.mat','image_data','params','schedule_info');
