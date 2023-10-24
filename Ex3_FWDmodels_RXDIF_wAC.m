%{
This script walks through the creation of a ROM and its usage to solve a forward model
The reaction-diffusion model is used as an example, with logistic growth
and cell death due to chemotherapy

dn/dt = D*d^2n/h^2 + kp*N(1-N) - alpha*C*N*sum_{i=1}^{T}exp(-t*beta_i)

Section 1:
    - Loading in the data
Section 2:
    - Running the forward model with the full PDE
Section 3:
    - Building a reduced model based on snapshots
Section 4:
    - Running the reduced model

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
location = [pwd,back,'Data',back,'Ex3-4_patient.mat'];
tumor = loadData(location);

%Extract necessary information for forward model
N0 = tumor.N(:,:,1); %2D initial cell count
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

%% Run FOM with known parameters
dt = 0.50;

start = tic;
[N_FOM,TC_FOM] = RXDIF_2D_wAC_comb(N0, kp.*ones(sy,sx), d, alpha, tx_params, t(2:end), h, dt, bcs);
t_FOM = toc(start); %time to run in seconds

disp(['FOM runtime = ',num2str(t_FOM),' sec']);
fprintf('\n');
%% Build ROM
start = tic;
%Step one is to prep the snapshot data for SVD by filling gaps in imaging data
%Here we will do that with the known parameters
N_augmented = augmentCellMaps_2D_knownParams_wAC(tumor.N, tumor.t_scan, h, dt, kp, d, alpha, tx_params, bcs);

%Basis is built using SVD
[V,k] = getProjectionMatrix(N_augmented, 0);

%Use to basis to build the ROM
%Need operators that represent the diffusion and proliferation terms
A = assembleA(N0, d, h, [], bcs);
B = assembleB(N0, kp);
% H = assembleH(N0, kp);
T = assembleT(N0, alpha);

%Reduce these with a Galerkin Projection
A_r = V' * A * V;
B_r = V' * B * V;
T_r = V' * T * V;

% H_r = V' * H * kron(V,V); %This is the exact method, but is a large matrix that typically can't be solved
%Can approximate with 100% accuracy since we have linear structure in the full operator
H_r = zeros(k,k^2);
for i = 1:size(V,1)
    H_r = H_r + V(i,:)' * kp * kron(V(i,:), V(i,:));
end

%Reduce initial conditions
N0_r = V' * N0(:);

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
%% Run ROM with known parameters
start = tic;
N_ROM_r = OperatorRXDIF_2D_wAC_comb(N0_r, A_r, B_r, H_r, T_r, tx_params, t(2:end), dt);
t_online = toc(start);

%Must resize ROM simulations to compare to full size
N_ROM = zeros(size(N_FOM));
for i = 1:numel(t(2:end))
    N_ROM(:,:,i) = reshape(V*N_ROM_r(:,i),sy,sx);
end

% Compare to measured data and FOM simulation
nt = size(tumor.N, 3);
figure
for i = 1:nt
    subplot(3,nt,i)
    imagesc(tumor.N(:,:,i), [0,1]); 
    axis image; colorbar; title(['Visit ',num2str(i)]);
    if i==1,  ylabel('Measured data'), end
    
    if i~=1
        subplot(3,nt,i+nt)
        imagesc(N_FOM(:,:,i-1), [0,1]);
        axis image; colorbar;
        if i==2,  ylabel('FOM simualtion'), end
        
        subplot(3,nt,i+2*nt)
        imagesc(N_ROM(:,:,i-1), [0,1]);
        axis image; colorbar;
        if i==2,  ylabel('ROM simualtion'), end
    end
end

disp(['ROM runtime time = ',num2str(t_online),' sec']);
fprintf('\n');

%Compare timing from the two methods
disp(['Runtime % change = ',num2str(100*(t_online - t_FOM)/t_FOM),'%']);
disp(['Total time % change = ',num2str(100*((t_offline+t_online) - t_FOM)/t_FOM),'%']);


%{
Notes: due to the offline stage the ROM is slower than a single full order solve

If we have to perform at least 2 runs of the 2D FOM, the time added by building
the ROM, is saved on the back-end due to the speed of the ROM

The number of forward runs required for ROM to surpass FOM is larger in 3D
%}