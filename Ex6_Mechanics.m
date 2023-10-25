%{

This script walks through the a reduced formualtion for determining the
mechanical stress in the tumor microenvironment, based on linear elastic
tissue assumptions

grad(G*grad(U))+grad((G/(1-2*nu))*grad(U)) = lambda*grad(N)
G - shear modulus
U - Tissue displacements
nu - Poisson's ratio
lambda - coupling factor
N - Cell maps

Section 1:
    - Loading in the data
Section 2:
    - Full order stress solve
Section 2:
    - Building reduced order model for mechanics
Section 3:
    - Reduced order stress solve

Contributors: Chase Christenson

Mechanics formulation references:
(1) Weis, Jared A et al. 
“A Mechanically Coupled Reaction-Diffusion Model for Predicting the 
Response of Breast Tumors to Neoadjuvant Chemotherapy.” 
Physics in medicine & biology 58.17 (2013): 5851–5866.

(2) Hormuth, David A. et al. 
“A Mechanically Coupled Reaction–diffusion Model That Incorporates 
Intra-Tumoural Heterogeneity to Predict in Vivo Glioma Growth.” 
Journal of the Royal Society interface 14.128 (2017): 20161010–20161010.

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
N_true = tumor.N(:,:,2:end); %Need full time course to build ROM

[sy,sx] = size(N0);

h = tumor.h; %Grid spacing comes directly from imaging resolution

bcs = tumor.bcs;
Tissues = tumor.Tissues;

t = tumor.t_scan;


%% Full order stress solve
start = tic;
[M,E,nu] = mech_matrix_build_2D(h, Tissues, bcs);
[d_dX, d_dY] = grad_matrix(h, bcs);
t_FOM_build = toc(start);

start = tic;
[~,sig_vm] = get_damper(d_dX, d_dY, N0, M, E, nu);
t_FOM_calc = toc(start);

disp(['FOM build time = ',num2str(t_FOM_build),' sec']);
disp(['FOM stress calculation time = ',num2str(t_FOM_calc),' sec']);
fprintf('\n');
%% Build reduced model
start = tic;
%Step 1: prep the snapshot data for SVD by filling gaps in imaging data
%Since we don't know parameters we fill by averaging and smoothing
N_augmented = augmentCellMaps_2D(cat(3,N0,N_true), t(2:end), 4); %Last parameter is depth of augmentation

%Step 2: Get relevant outputs from stress calculation at each time point in the augmented dataset
% Displacement maps (each direction), strain maps (principal and shear), stress maps (principal and shear)
[~, ~, Ux_aug, Uy_aug, Exx_aug, Eyy_aug, Exy_aug, Sxx_aug, Syy_aug, Sxy_aug] = getMechanicsMaps_2D(N_augmented, M, E, nu, d_dX, d_dY);

%Step 3: Get projection matrix for each portion of the stress calculation
[V,k] = getProjectionMatrix(N_augmented, 0); %rank is determined by cell map reduction

[V_u, V_Ux, V_Uy] = getDisplacementProjection_2D(Ux_aug, Uy_aug, k); 
%V_Ux and V_Uy are needed later, so must be saved individually as well as bunched into V_u
clear Ux_aug Uy_aug;

V_e = getStrainProjection_2D(Exx_aug, Eyy_aug, Exy_aug, k);
clear Exx_aug Eyy_aug Exy_aug;

V_s = getStressProjection_2D(Sxx_aug, Syy_aug, Sxy_aug, k);
clear Sxx_aug Syy_aug Sxy_aug;

%Step 4: Build and reduce mechanics matrices for solver
[M,E,nu] = mech_matrix_build_2D(h, Tissues, bcs);
[d_dX, d_dY] = grad_matrix(h, bcs);
S_mat = buildStrainMat(N0,E,nu);

%Reduced operator for displacement solve
M_r = V_u' * M * V_u; 

%Reduced operator for cell gradient calculation
Vu_gradXY_V = V_u' * [d_dX, zeros(size(d_dX)); zeros(size(d_dX)), d_dY] * [V, zeros(size(V)); zeros(size(V)), V]; 

%Reduced operator for strain solve
Ve_gradXYY_Vu = V_e' * [d_dX, zeros(size(d_dX)), zeros(size(d_dX)); ...
                        zeros(size(d_dX)), d_dY, zeros(size(d_dX)); ...
                        zeros(size(d_dX)), zeros(size(d_dX)), d_dY] * ...
                        [V_Ux, zeros(size(V_Uy)), zeros(size(V_Ux)); zeros(size(V_Ux)), V_Uy, zeros(size(V_Ux)); zeros(size(V_Ux)), zeros(size(V_Uy)), V_Ux];

%Reduced operator for stress solve
Vs_SMAT_Ve = V_s' * S_mat * V_e;
t_ROM_build = toc(start);


clear S_mat V_Ux V_Uy V_u V_e M; %None of these are needed for the forward solver

%% Reduced order stress solve
%Input to reduced stress solver is reduced cell map
%Output is full sized stress map
N0_r = V' * N0(:);

start = tic;
[~,sig_vm_ROM] = get_damper_reduced(Vu_gradXY_V, Ve_gradXYY_Vu, Vs_SMAT_Ve, N0_r, M_r, V_s);
t_ROM_calc = toc(start);

disp(['ROM build time = ',num2str(t_ROM_build),' sec']);
disp(['ROM stress calculation time = ',num2str(t_ROM_calc),' sec']);
fprintf('\n');


%% Comparison
%Visualization
figure
subplot(1,3,1)
imagesc(reshape(sig_vm,size(N0)));
axis image; axis off; title('Sigma_v_m FOM');
subplot(1,3,2)
imagesc(reshape(sig_vm_ROM,size(N0))); 
axis image; axis off; title('Sigma_v_m ROM');
subplot(1,3,3)
scatter(sig_vm, sig_vm_ROM, 50);
axis square; xlabel('FOM Calculation'); ylabel('ROM Calculation');
refline(1,0);
title(['CCC = ',num2str(CCC_calc(sig_vm, sig_vm_ROM),'%.3f')]);

disp(['Build time % change = ',num2str(100*(t_ROM_build - t_FOM_build)/t_FOM_build),'%']);
disp(['Calculation time % change = ',num2str(100*(t_ROM_calc - t_FOM_calc)/t_FOM_calc),'%']);
fprintf('\n');

%%% Can try with different times from the augmented dataset or true data %%%

%{
Notes: The offline stage for the mechanics model is long compared to the build
for the non-mechanically coupled counter part

The calculation is around 30X faster in 2D, this has a larger benefit in 3D

For a single run, FOM is still faster but when we are trying to calibrate
parameters, the ROM is still the optimal choice. A similar calibration to
Example 5, but with mechanics, is 40X faster in 2D, and "BLANK"X faster in 3D
when using ROM. 

%}