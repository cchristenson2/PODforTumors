function [damper, s_vm, Ux, Uy, Uz] = get_damper_3D(matX, matY, matZ, N, M, E, nu)
    [sy,sx,sz] = size(N);

    lambda1 = 2.5e-3;
    lambda2 = 2.5e-3;
    
    n_vect = reshape(N,[],1);
    % calculate the gradient using matX and matY
    gradX = matX * n_vect;
    gradY = matY * n_vect;
    gradZ = matZ * n_vect;
    grad_N = [gradX; gradY; gradZ];

    % Calculate displacement
%     Uv = M\(lambda1*grad_N);

%     disp(size(N));
%     disp(size(M));
%     disp(size(grad_N));

    [Uv, flags] = bicgstab(M, lambda1*grad_N, 1e-6);

    % Reshape U into matrix
%     Um = zeros(sy,sx,sz,3);
%     Um(:,:,1) = reshape(Uv(1:sx*sy*sz), sy,sx,sz);
%     Um(:,:,2) = reshape(Uv(sx*sy*sz+1:sx*sy*sz*2),sy,sx,sz);
%     Um(:,:,3) = reshape(Uv(sx*sy*sz*2+1:end),sy,sx,sz);
    

    Um = zeros(sy*sx*sz,3);
    Um(:,1) = Uv(1:sx*sy*sz);            Ux = reshape(Um(:,1), sy, sx, sz);
    Um(:,2) = Uv(sx*sy*sz+1:sx*sy*sz*2); Uy = reshape(Um(:,2), sy, sx, sz);
    Um(:,3) = Uv(sx*sy*sz*2+1:end);      Uz = reshape(Um(:,3), sy, sx, sz);



%     [e_xx,e_yy,e_zz,e_xy,e_xz,e_yz] = strains(Um, h, h, dz);
    e_xx = reshape(matX * Um(:,1),sy,sx,sz);
    e_yy = reshape(matY * Um(:,2),sy,sx,sz);
    e_zz = reshape(matZ * Um(:,3),sy,sx,sz);
    e_xy = reshape((matY * Um(:,1) + matX * Um(:,2))/2,sy,sx,sz);
    e_xz = reshape((matZ * Um(:,1) + matX * Um(:,3))/2,sy,sx,sz);
    e_yz = reshape((matZ * Um(:,2) + matY * Um(:,3))/2,sy,sx,sz);
    
    [s_xx,s_yy,s_zz,s_xy,s_xz,s_yz,s_vm] = stresses_3D(e_xx,e_yy,e_zz,e_xy,e_xz,e_yz, E, nu);
    
%     s_xx = reshape(s_xx, sy,sx,sz); s_yy = reshape(s_yy, sy,sx,sz); s_zz = reshape(s_zz, sy,sx,sz);
%     s_xy = reshape(s_xy, sy,sx,sz); s_xz = reshape(s_xz, sy,sx,sz); s_yz = reshape(s_yz, sy,sx,sz);
%     s_vm = reshape(s_vm, sy,sx,sz);
    
    
%     slice = round(sz/2);
%     figure
%     subplot(4,3,1)
%     imagesc(Ux(:,:,slice)); title('U_x'); axis image;
%     subplot(4,3,2)
%     imagesc(Uy(:,:,slice)); title('U_y'); axis image;
%     subplot(4,3,3)
%     imagesc(Uz(:,:,slice)); title('U_z'); axis image;
%     
%     subplot(4,3,4)
%     imagesc(s_xx(:,:,slice)); title('sig_x_x'); axis image;
%     subplot(4,3,5)
%     imagesc(s_yy(:,:,slice)); title('sig_y_y'); axis image;
%     subplot(4,3,6)
%     imagesc(s_zz(:,:,slice)); title('sig_z_z'); axis image;
%     
%     subplot(4,3,7)
%     imagesc(s_xy(:,:,slice)); title('sig_x_y'); axis image;
%     subplot(4,3,8)
%     imagesc(s_xz(:,:,slice)); title('sig_x_z'); axis image;
%     subplot(4,3,9)
%     imagesc(s_yz(:,:,slice)); title('sig_y_z'); axis image;
%     
%     subplot(4,3,10)
%     imagesc(s_vm(:,:,slice)); title('sig_v_m'); axis image;
%     subplot(4,3,11)
%     imagesc(N(:,:,slice)); title('Cell Count'); axis image;
%     
%     figure
%     subplot(1,3,1)
%     imagesc(max(abs(Ux),[],3)); title('Max U_x'); axis image;
%     subplot(1,3,2)
%     imagesc(max(abs(Uy),[],3)); title('Max U_y'); axis image;
%     subplot(1,3,3)
%     imagesc(max(abs(Uz),[],3)); title('Max U_z'); axis image;
    
    damper = exp(-1*s_vm * lambda2);
end