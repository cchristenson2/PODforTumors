data = loadPatientData_coarseRes('data_282.mat');


N = data.N(:,:,:,1);
[sy,sx,sz] = size(N);
bcs = buildBoundaries_3D(data.Mask);
h = data.h;
dz = data.dz;
tissue = data.Tissues;

[matX,matY,matZ] = gradN_matrix(sx,sy,sz,h,dz,bcs);
[M, E, nu] = mech_matrix_build_3D_v1(h, dz, tissue, bcs);


[damper, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz] = get_damper(matX, matY, matZ, N, M, E, nu, h, dz);


figure
spy(M);


[M2D, ~,~] =  mech_matrix_build_2D_v5(h,tissue(:,:,8), buildBoundaries_2D(data.Mask(:,:,8)));


figure
subplot(1,2,1)
spy(M); title('3D');

subplot(1,2,2)
spy(M2D); title('2D');