function [damper, s_vm, Ux, Uy, Uz] = get_damper_reduced_3D(matX, matY, matZ, grad_N, M_r, E, nu, Vs)
    lambda1 = 2.5e-3;
    lambda2 = 2.5e-3;
    

    % Calculate displacement
    [Uv_r, flags] = bicgstab(M_r, lambda1*grad_N, 1e-6);

    % Reshape U into matrix
    Uv = Vs * Uv_r;

    num = numel(Uv);
    Um = zeros(num/3,3);
    Um(:, 1) = Uv(1:num/3);           Ux = Um(:,1);
    Um(:, 2) = Uv(num/3+1:2*(num/3)); Uy = Um(:,2);
    Um(:, 3) = Uv(2*(num/3)+1:end);   Uz = Um(:,3);


    e_xx = matX * Um(:,1);
    e_yy = matY * Um(:,2);
    e_zz = matZ * Um(:,3);
    e_xy = matY * Um(:,1);
    e_xz = matZ * Um(:,1);
    e_yz = matZ * Um(:,2);
%     e_xy = (matY * Um(:,1) + matX * Um(:,2))/2;
%     e_xz = (matZ * Um(:,1) + matX * Um(:,3))/2;
%     e_yz = (matZ * Um(:,2) + matY * Um(:,3))/2;
    
    [s_xx,s_yy,s_zz,s_xy,s_xz,s_yz,s_vm] = stresses_3D(e_xx,e_yy,e_zz,e_xy,e_xz,e_yz, E(:), nu);

    damper = exp(-1*s_vm * lambda2);
end