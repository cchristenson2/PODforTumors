
%{

This script walks through the creation of a ROM and its calibrate parameters to data
The reaction-diffusion model is used as an example, with logistic growth

dn/dt = D*d^2n/h^2 + kp*N(1-N)

Section 1:
    - Loading in the data
Section 2:
    - Calibrating parameters using the FOM
Section 3:
    - Building a library of reduced models based on snapshots
Section 4:
    - Calibrating parameters using the ROM

Contributors: Chase Christenson
%}


%% Load data
%Add all functions for ROM into the current path
clear; clc; close all;

addpath(genpath(pwd))
if(ispc)
    back = '\';
else
    back = '/';
end

%We use example data that has been processed and has correct formatting for loading
location = [pwd,back,'Data',back,'Ex1-2_patient.mat'];
tumor = loadData(location);

%Extract necessary information for forward model
N0 = tumor.N(:,:,1); %2D initial cell count
N_true = tumor.N(:,:,2:end);

[sy,sx] = size(N0);

kp = tumor.params.kp; %true parameters are known for the virtual patient
d = tumor.params.d;

h = tumor.h; %Grid spacing comes directly from imaging resolution

bcs = tumor.bcs;

t = tumor.t_scan;

%% Calibrate parameters with FOM
dt = 0.50;

bounds = [1e-6, 5e-1; 1e-3, 1e-1]; %Diffusivity; proliferation

start = tic;
[params_FOM, N_cal_FOM] = calibrateRXDIF_FOM(N0, N_true, bounds, t(2:end), h, dt, bcs, 1, []);
t_FOM = toc(start); %time to run in seconds
disp(['FOM calibration time = ',num2str(t_FOM),' sec']);

%Compare calibrated cell maps to the measured data
for i = 1:size(N_true,3)
    CCC = CCC_calc(N_true(:,:,i), N_cal_FOM(:,:,i));
    disp(['Visit ',num2str(i+1),' CCC between FOM and measured = ',num2str(CCC,'%.3f')]);
end
fprintf('\n');

%Compare output parameters to the true values
disp(['Prolfieration % error for FOM = ',num2str(100*(params_FOM.kp - tumor.params.kp)/tumor.params.kp, '%.3f')]);
disp(['Diffusivity % error for FOM = ',num2str(100*(params_FOM.d - tumor.params.d)/tumor.params.d, '%.3f')]);
fprintf('\n');

%% Build ROM
start = tic;
%First build operators that represent the bounds of the potential parameters
param_types = {'A',1; 'B',2};
fmt = '%.6f'; %naming convention that gets used for operator indexing

Lib = buildLibrary(bounds,N0,h,[],bcs,param_types,fmt,2); 
%Notes - Number of operators to build for the range spanned by bounds is last value
%      - [] is where dz, or slice spacing would go for 3D models

%Step 2: prep the snapshot data for SVD by filling gaps in imaging data
%Since we don't know parameters we fill by averaging and smoothing
N_augmented = augmentCellMaps_2D(cat(3,N0,N_true), t(2:end), 4); %Last parameter is depth of augmentation

%%Step 3: Basis is built using SVD
[V,k] = getProjectionMatrix(N_augmented, 0);

%Step 4: Reduce the library with the new basis
Lib_r = reduceLibrary(Lib, V, param_types, []);
%Notes - Last [] represents the concentration map, which is needed for treatment models

% H-library is built individually since multiple N X N^2 matrices are memory intensive
%Approximate method builds reduced operators directly, without middle step to build the full operator
Lib_r.Hr_lib = buildAndReduce_GlobalKp(bounds(2,:),V,fmt,2,'H');

%Reduce initial conditions
N0_r = V' * N0(:);

%Reduce data for calibration;
N_true_r = zeros(k,size(N_true,3));
for i = 1:size(N_true,3)
    N_true_r(:,i) = V' * reshape(N_true(:,:,i),[],1);
end

t_offline = toc(start);

%Visualize the modes for reduction
figure
for i = 1:k
    subplot(1,k,i)
    imagesc(reshape(V(:,i),sy,sx),[min(V(:)), max(V(:))]);
    axis image; axis off; colorbar;
    title(['Mode shape ',num2str(i)]);
end

disp(['ROM offline build time = ',num2str(t_offline),' sec']);
fprintf('\n');
%% Calibrate parameters with ROM
start = tic;
[params_ROM, N_cal_ROM] = calibrateRXDIF_ROM(N0_r, N_true_r, Lib_r, bounds, t(2:end), dt, 1, []);
t_ROM = toc(start); %time to run in seconds
disp(['ROM calibration time = ',num2str(t_ROM),' sec']);

%Compare calibrated cell maps to the measured data
for i = 1:size(N_true,3)
    CCC = CCC_calc(N_true(:,:,i), reshape(V*N_cal_ROM(:,i),sy,sx));
    disp(['Visit ',num2str(i+1),' CCC between ROM and measured = ',num2str(CCC,'%.3f')]);
end
fprintf('\n');

%Compare output parameters to the true values
disp(['Prolfieration % error for ROM = ',num2str(100*(params_ROM.kp - tumor.params.kp)/tumor.params.kp, '%.3f')]);
disp(['Diffusivity % error for ROM = ',num2str(100*(params_ROM.d - tumor.params.d)/tumor.params.d, '%.3f')]);
fprintf('\n');

disp(['Total time % change = ',num2str(100*((t_offline+t_ROM) - t_FOM)/t_FOM),'%']);

%% Visualize calibration
nt = size(tumor.N, 3);
figure
for i = 1:nt
    subplot(3,nt,i)
    imagesc(tumor.N(:,:,i), [0,1]); 
    axis image; colorbar; title(['Visit ',num2str(i)]);
    if i==1,  ylabel('Measured data'), end
    
    if i~=1
        subplot(3,nt,i+nt)
        imagesc(N_cal_FOM(:,:,i-1), [0,1]);
        axis image; colorbar;
        if i==2,  ylabel('FOM calibration'), end
        
        subplot(3,nt,i+2*nt)
        imagesc(reshape(V*N_cal_ROM(:,i-1),sy,sx), [0,1]);
        axis image; colorbar;
        if i==2,  ylabel('ROM calibration'), end
    end
end


%{
Notes: The total time for ROM now decreases compared to the total FOM time

The proliferation parameter consistently calibrates with high accuracy, but 
the ROM struggles to capture the true diffusivity. Diffusion is still 
visible in the ROM outputs but is primarily captured through mode shapes.
%}