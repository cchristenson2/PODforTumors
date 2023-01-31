function [s_x,s_y,s_z,s_xy,s_xz,s_yz,s_vm] = stresses(e_xx,e_yy,e_zz,e_xy,e_xz,e_yz, E, v)
% Factor we need to multiply by
mult_fac = E/(1+v)/(1-2*v);
% Matrix we need to multiply by
strain_mat = [1-v, v, v, 0,0,0; v, 1-v, v, 0,0,0; v, v, 1-v, 0,0,0; 0,0,0, 1-2*v, 0,0; 0,0,0,0, 1-2*v, 0; 0,0,0,0,0, 1-2*v];

% get our sizes
[sy,sx,sz] = size(e_xx);

% Initialize
s_x = zeros(sy,sx,sz);
s_y = zeros(sy,sx,sz);
s_z = zeros(sy,sx,sz);
s_xy = zeros(sy,sx,sz);
s_xz = zeros(sy,sx,sz);
s_yz = zeros(sy,sx,sz);
s_vm = zeros(sy,sx,sz);

% Loop through all points
for k = 1:sz
    for i = 1:sx
        for j = 1:sy
            % Obtain our stress vector
            s_vec = mult_fac(j,i,k) * strain_mat * [e_xx(j,i,k); e_yy(j,i,k); e_zz(j,i,k); e_xy(j,i,k); e_xz(j,i,k); e_yz(j,i,k)];
            % Store each component
            s_x(j,i,k) = s_vec(1);
            s_y(j,i,k) = s_vec(2);
            s_z(j,i,k) = s_vec(3);
            s_xy(j,i,k) = s_vec(4);
            s_xz(j,i,k) = s_vec(5);
            s_yz(j,i,k) = s_vec(6);

            % Calculate and store von mises stress
            s_vm(j,i,k) = sqrt(0.5*((s_vec(1)-s_vec(2))^2 + (s_vec(1)-s_vec(3))^2 + (s_vec(2)-s_vec(3))^2 + 6*(s_vec(4)^2 + s_vec(5)^2 + s_vec(6)^2)));
        end
    end
end

end