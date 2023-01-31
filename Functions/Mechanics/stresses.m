function [s_x, s_y, s_xy, s_vm] = stresses(e_xx,e_yy, e_xy, E, v)
% Factor we need to multiply by
mult_fac = E/(1+v)/(1-2*v);
% Matrix we need to multiply by
strain_mat = [1-v, v, 0; v, 1-v, 0; 0, 0, 1-2*v];

% get our sizes
[sy,sx] = size(e_xx);

% Initialize
s_x = zeros(sy,sx);
s_y = zeros(sy,sx);
s_xy = zeros(sy,sx);
s_vm = zeros(sy,sx);

% Loop through all points
for i = 1:sx
    for j = 1:sy
        % Obtain our stress vector
        s_vec = mult_fac(j,i) * strain_mat * [e_xx(j,i); e_yy(j,i); e_xy(j,i)];
        % Store each component
        s_x(j,i) = s_vec(1);
        s_y(j,i) = s_vec(2);
        s_xy(j,i) = s_vec(3);
        % Calculate and store von mises stress
        s_vm(j,i) = sqrt(s_vec(1)^2-s_vec(1)*s_vec(2)+s_vec(2)^2+3*s_vec(3)^2);
    end
end
end