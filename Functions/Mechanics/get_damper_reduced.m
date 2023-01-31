function [damper, s_vm, Ux, Uy] = get_damper_reduced(matX, matY, grad_N, M_r, E, nu, Vs, transVx_Vxx, transVy_Vyy, transVx_Vxy, Vxx, Vyy, Vxy, matX_rxx, matY_ryy, matY_rxy)
%     [sy,sx] = size(N);

    lambda1 = 2.5e-3;
    lambda2 = 2.5e-3;
    
%     n_vect = reshape(N_r,[],1);
    % calculate the gradient using matX and matY
%     grad_N = matXY_r * n_vect;

    % Calculate displacement
    start = tic;

    Uv_r = M_r\(lambda1.*grad_N);
    
    back_solve = toc(start);
    
    %Resize back to full order
    
    start = tic;

%     Uv = Vs * Uv_r;
    
    Uv = Uv_r;
    
    num = numel(Uv);
    Um = zeros(num/2,2);
    Um(:, 1) = Uv(1:num/2); Ux = Um(:,1);
    Um(:, 2) = Uv(num/2+1:end); Uy = Um(:,2);
    
    displacement = toc(start);


%     [e_xx,e_yy,e_xy] = strains(Um, h, h);
    
    start = tic;

%     e_xx = matX*Um(:,1);
%     e_yy = matY*Um(:,2);
%     e_xy = matY*Um(:,1);
    
    e_xx = matX_rxx * (transVx_Vxx' * Um(:,1));
    e_yy = matY_ryy * (transVy_Vyy' * Um(:,2));
    e_xy = matY_rxy * (transVx_Vxy' * Um(:,1));
    
    e_xx = Vxx * e_xx;
    e_yy = Vyy * e_yy;
    e_xy = Vxy * e_xy;
    
    strain = toc(start);
    
    start = tic;

    [s_xx, s_yy, s_xy, s_vm] = stresses(e_xx,e_yy, e_xy, E(:), nu);
    
    stress = toc(start);
    
    disp(['Solve time = ', num2str(back_solve)]);
    disp(['Displacement time = ', num2str(displacement)]);
    disp(['Strain time = ', num2str(strain)]);
    disp(['Stress time = ', num2str(stress)]);
    
    damper = exp(-1*s_vm * lambda2);
end