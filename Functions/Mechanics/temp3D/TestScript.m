data = loadPatientData_coarseRes('data_282.mat');


N = data.N(:,:,:,1);
[sy,sx,sz] = size(N);
bcs = buildBoundaries_3D(data.Mask);
h = data.h;
dz = data.dz;
tissue = data.Tissues;

[matX,matY,matZ] = gradN_matrix_3D(sx,sy,sz,h,dz,bcs);
[M, E, nu] = mech_matrix_build_3D_v1(h, dz, tissue, bcs);


[damper, s_vm, Ux, Uy, Uz] = get_damper_3D(matX, matY, matZ, N, M, E, nu);

slice = round(sz/2);
figure
subplot(2,3,1)
imagesc(damper(:,:,slice)); axis image; axis off; title('Damper');
subplot(2,3,2)
imagesc(s_vm(:,:,slice)); axis image; axis off; title('Sig_v_m');
% subplot(2,3,3)
% imagesc(s_zz(:,:,slice)); axis image; axis off; title('Ux');
subplot(2,3,4)
imagesc(Ux(:,:,slice)); axis image; axis off; title('U_x');
subplot(2,3,5)
imagesc(Uy(:,:,slice)); axis image; axis off; title('U_y');
subplot(2,3,6)
imagesc(Uz(:,:,slice)); axis image; axis off; title('U_z');


figure
spy(M);


[M2D, ~,~] =  mech_matrix_build_2D_v5(h,tissue(:,:,8), buildBoundaries_2D(data.Mask(:,:,8)));


figure
subplot(1,2,1)
spy(M); title('3D');

subplot(1,2,2)
spy(M2D); title('2D');