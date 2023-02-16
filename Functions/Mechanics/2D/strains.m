function [e_xx,e_yy,e_xy] = strains(Um, dx, dy)
% Split U matrix into x and y displacements
Ux = Um(:,:,1);
Uy = Um(:,:,2);

% Get size
[sy,sx] = size(Ux);

% Initialize matrices to store strain
e_xx = zeros(sy,sx);
e_yy = zeros(sy,sx);
e_xy = zeros(sy,sx);

% Loop through all points
for i = 1:sx
    for j = 1:sy
        % Perform forward difference with order of accuracy equal to our
        % central difference, so using two points froward (and for
        % backwards use two points backwards)
        
        if i == 1 % forwards difference
            duxdx = (-3*Ux(j,i) + 4*Ux(j, i+1) - Ux(j, i+2))/(2*dx);
            duydx = (-3*Uy(j,i) + 4*Uy(j, i+1) - Uy(j, i+2))/(2*dx);
        elseif i == sx % backwards difference
           duxdx = (3*Ux(j,i) - 4*Ux(j, i-1) + Ux(j, i-2))/(2*dx);
            duydx = (3*Uy(j,i) - 4*Uy(j, i-1) + Uy(j, i-2))/(2*dx); 
        else % central difference
            duxdx = (Ux(j,i+1) - Ux(j, i-1))/(2*dx);
            duydx = (Uy(j,i+1) - Uy(j, i-1))/(2*dx);
        end

        if j == 1 % backwards difference
            duxdy = (3*Ux(j,i) - 4*Ux(j+1, i) + Ux(j+2, i))/(2*dy);
            duydy = (3*Uy(j,i) - 4*Uy(j+1, i) + Uy(j+2, i))/(2*dy);
        elseif j == sy % forwards difference
            duxdy = (-3*Ux(j,i) + 4*Ux(j-1, i) - Ux(j-2, i))/(2*dy);
            duydy = (-3*Uy(j,i) + 4*Uy(j-1, i) - Uy(j-2, i))/(2*dy); 
        else % central difference
            duxdy = (Ux(j+1,i) - Ux(j-1, i))/(2*dy);
            duydy = (Uy(j+1,i) - Uy(j-1, i))/(2*dy);
        end
        
        % form strains
        e_xx(j,i) = duxdx;
        e_yy(j,i) = duydy;
        e_xy(j,i) = (duxdy + duydx) / 2;

    end
end

end