
%{

This script walks through the creation of a ROM and its calibrate parameters to data
The reaction-diffusion model is used as an example, with logistic growth
- The proliferation map is now a spatially varying field

dn/dt = D*d^2n/h^2 + kp(x)*N(1-N)

Section 1:
    - Loading in the data
Section 2:
    - Building a library of reduced models based on snapshots
Section 3:
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
location = [pwd,back,'Data',back,'Ex5_patient.mat'];
tumor = loadData(location);

%Extract necessary information for forward model
N0 = tumor.N(:,:,1); %2D initial cell count
N_true = tumor.N(:,:,2:end);

[sy,sx] = size(N0);

kp = tumor.params.kp; %true parameters are known for the virtual patient
d = tumor.params.d;
alpha = tumor.params.alpha;

tx_params.txduration = tumor.t_trx;
tx_params.beta1 = tumor.beta1;
tx_params.beta2 = tumor.beta2;
tx_params.C     = tumor.AUC;

h = tumor.h; %Grid spacing comes directly from imaging resolution

bcs = tumor.bcs;

t = tumor.t_scan;

%% Calibrate parameters with FOM
%Not shown, lengthy calibration without HPC and highly parallel code

%% Build ROM
dt = 0.50;
bounds = [1e-6, 1e-2; 1e-3, 1e-1; 1e-6, 0.8]; %Diffusivity; proliferation; alpha

start = tic;
%First build operators that represent the bounds of the potential parameters (Global parameters only)
param_types = {'A',1; 'T',3};
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
Lib_r = reduceLibrary(Lib, V, param_types, tx_params.C);
%Notes - Last [] represents the concentration map, which is needed for treatment models

% B and H libraries represent spatially varying parameters. To account for this,
% libraries are built for each mode in the basis, V
[Lib_r.Br_lib, kp_mode_bounds] = buildAndReduce_LocalKp(bounds(2,:),V,fmt,3,'B');
Lib_r.Hr_lib = buildAndReduce_LocalKp(bounds(2,:),V,fmt,3,'H');

%Step 5: Reduce data for fwd evaluations and comparisons
%Reduce initial conditions
N0_r = V' * N0(:);

%Reduce data for calibration;
N_true_r = zeros(k,size(N_true,3));
for i = 1:size(N_true,3)
    N_true_r(:,i) = V' * reshape(N_true(:,:,i),[],1);
end

t_offline = toc(start);

disp(['ROM offline build time = ',num2str(t_offline),' sec']);
fprintf('\n');
%% Calibrate parameters with ROM
start = tic;
[params_ROM, N_cal_ROM] = calibrateRXDIF_localKp_ROM(N0_r, N_true_r, Lib_r, bounds, kp_mode_bounds, t(2:end), dt, 2, tx_params, V);
t_ROM = toc(start); %time to run in seconds
disp(['ROM calibration time = ',num2str(t_ROM),' sec']);

%Compare calibrated cell maps to the measured data
for i = 1:size(N_true,3)
    CCC = CCC_calc(N_true(:,:,i), reshape(V*N_cal_ROM(:,i),sy,sx));
    disp(['Visit ',num2str(i+1),' CCC between ROM and measured = ',num2str(CCC,'%.3f')]);
end
fprintf('\n');

%Compare output parameters to the true values
disp(['Diffusivity % error for ROM = ',num2str(100*(params_ROM.d - tumor.params.d)/tumor.params.d, '%.3f')]);
disp(['Alpha % error for ROM = ',num2str(100*(params_ROM.alpha - tumor.params.alpha)/tumor.params.alpha, '%.3f')]);
idx = intersect(find(abs(kp)>1e-3), find(abs(params_ROM.kp)>1e-3));
disp(['Proliferation CCC in overlap for ROM = ',num2str(CCC_calc(kp(idx), params_ROM.kp(idx)), '%.3f')]);
fprintf('\n');

%% Visualize calibration
nt = size(tumor.N, 3);
figure
for i = 1:nt
    subplot(2,nt,i)
    imagesc(tumor.N(:,:,i), [0,1]); 
    axis image; colorbar; title(['Visit ',num2str(i)]);
    if i==1,  ylabel('Measured data'), end
    
    if i~=1

        subplot(2,nt,i+nt)
        imagesc(reshape(V*N_cal_ROM(:,i-1),sy,sx), [0,1]);
        axis image; colorbar;
        if i==2,  ylabel('ROM calibration'), end
    end
end

figure
subplot(1,3,1)
imagesc(kp, bounds(2,:)); axis image; title('True Proliferation');
subplot(1,3,2)
imagesc(reshape(params_ROM.kp,size(N0)), bounds(2,:)); axis image; title('ROM Proliferation');
subplot(1,3,3)
imagesc(kp - reshape(params_ROM.kp,size(N0))); axis image; title('True - ROM'); colorbar;


%{
Notes: The total time for ROM now decreases compared to the total FOM time

The proliferation parameter consistently calibrates with high accuracy, but 
the ROM struggles to capture the true diffusivity. Diffusion is still 
visible in the ROM outputs but is primarily captured through mode shapes.
%}