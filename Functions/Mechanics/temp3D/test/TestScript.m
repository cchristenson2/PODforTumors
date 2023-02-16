data = loadPatientData_coarseRes('data_282.mat');


N = data.N(:,:,:,1);
[sy,sx,sz] = size(N);
bcs = buildBoundaries_3D(data.Mask);
h = data.h;
dz = data.dz;
tissue = data.Tissues;

[matX,matY,matZ] = gradN_matrix_3D(sx,sy,sz,h,dz,bcs);
[M, E, nu] = mech_matrix_build_3D_v1(h, dz, tissue, bcs);


[damper, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz] = get_damper_3D(matX, matY, matZ, N, M, E, nu, h, dz);

slice = round(sz/2);
figure
subplot(2,3,1)
imagesc(s_xx(:,:,slice)); axis image; axis off; title('Sig_x_x');
subplot(2,3,2)
imagesc(s_yy(:,:,slice)); axis image; axis off; title('Sig_y_y');
subplot(2,3,3)
imagesc(s_zz(:,:,slice)); axis image; axis off; title('Sig_z_z');
subplot(2,3,4)
imagesc(s_xy(:,:,slice)); axis image; axis off; title('Sig_x_y');
subplot(2,3,5)
imagesc(s_xz(:,:,slice)); axis image; axis off; title('Sig_x_z');
subplot(2,3,6)
imagesc(s_yz(:,:,slice)); axis image; axis off; title('Sig_y_z');


figure
spy(M);


[M2D, ~,~] =  mech_matrix_build_2D_v5(h,tissue(:,:,8), buildBoundaries_2D(data.Mask(:,:,8)));


figure
subplot(1,2,1)
spy(M); title('3D');

subplot(1,2,2)
spy(M2D); title('2D');