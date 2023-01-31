function [M, E, nu] =  mech_matrix_build_3D_v1(h, dz, tissue, bcs)
% h is our step in x and y ( h= h)
% dz is our step in z
% tissue is our tissue type and tumor map for this patient
% bcs is the map of boundary conditions 
    %  Boundary matrix bcs(y,x,z,3) = [y,x,z, [y boundary type, x boundary type, z boundary type]]
    %             Value = -1 if boundary is behind, left, above
    %             Value = 0 if no boundary is present
    %             Value = 1 if boundary is in front, right, below
    %             Value = 2 if outside of mask

% a bit of setup
[sy,sx,sz] = size(tissue);
M = zeros(sy*sx*sz*2,sy*sx*sz*2); % it big

% First thing, use our tissue map to set mechanics params
% we use one nu value, but multiple E values
nu = 0.45; 
E_adipose = 2e3;
E_fibro = 2*2e3;
E_tumor = 10*2e3;

E = zeros(size(tissue));
E(tissue == 1) = E_adipose;
E(tissue == 2) = E_fibro;
E(tissue == 0) = E_tumor;

% E = imgaussfilt(E,2);

E(bcs(:,:,:,1) ==2 ) = 0;



G = E / (2*(1+nu));
% And set up some variables for mechanics that we'll use often
% Ke = E/((1+nu)*(1-2*nu));
% Ks = (1-2*nu)/2;
% Km = 1-nu;
k1 = 2*(1-nu) / (1-2*nu);
k2 = 2*nu / (1-2*nu);

% Now we want to set up variables to make it easy to access where we should
% be in our matrix M based on our x,y coordinates. This will make a matrix
% that looks like
%[1,4,7
% 2,5,8
% 3,6,9]
% on an example 3x3 domain such that if we're looking for the y=1, x=3
% position in the matrix
%[1,2,3
% 4,5,6
% 7,8,9]
% flattened with matlab's standard reshape to [1,4,7,...] then we can take
% itm(1,3) and we obtain 7, which works because the entry at y=1,x=3 is at
% the 7th index in its flattened form
count = 0;
for z = 1:sz
    for x = 1:sx
        for y = 1:sy
            count = count + 1;
            itm(y,x,z) = count;
        end
    end
end

% set up mult factors for BCs
for z = 1:sz
    for y = 1:sy
        for x = 1:sx
            count = itm(y,x,z);
            county = itm(y,x,z)+sy*sx;
            countz = itm(y,x,z)+sy*sx*sz;
            %Ke = Ke_mat(y,x);
            % ARE WE IN THE BREAST
            if bcs(y,x,z,2) == 2
                M(count,count) = 1; % forces 0 displacement
                M(county,county) = 1; % forces 0 displacement
                M(countz,countz) = 1; % forces 0 displacement
                continue
            end

            % X DIRECTION -- COMPLETE
            % First, check to see if we are on an x boundary
            % If we are, set to 1
            if bcs(y,x,z,2) == 1 || bcs(y,x,z,2) == -1 %x == 1 || x == sx
                M(count,count) = 1; % forces 0 displacement
            else
                
                % X ONLY DERIVATIVES -- COMPLETE
                % TERM MULT BY K1
                % dG/dx * dUx/dx  +  G * d^2Ux/dx^2
                M(count,itm(y,x+1,z)) = k1 * (G(y,x+1,z)-G(y,x-1,z))/(4*h^2) + k1 * G(y,x,z)/h^2;

                M(count,itm(y,x,z))   = -2 * k1 * G(y,x,z)/h^2;

                M(count,itm(y,x-1,z)) = k1 * (G(y,x-1,z)-G(y,x+1,z))/(4*h^2) + k1 * G(y,x,z)/h^2;
                
                % Y & X DERIVATIVES -- COMPLETE
                % We are on top boundary (y-1 does not exist) -- COMPLETE
                if bcs(y,x,z,1) == -1 % y == 1
                    % TERM MULT BY K2
                    % dG/dx * dUy/dy  +  G * d^2Uy/(dx*dy)
                    M(count,itm(y+1,x,z)+(sy*sx)) = k2 * (G(y,x+1,z)-G(y,x-1,z))/(2*h^2);

                    M(count,itm(y,x,z)+(sy*sx)) = k2 * (G(y,x-1,z)-G(y,x+1,z))/(2*h^2);


                    M(count,itm(y+1,x+1,z)+(sy*sx)) = k2 * G(y,x,z) / (2*h^2);

                    M(count,itm(y+1,x-1,z)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(count,itm(y,x+1,z)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(count,itm(y,x-1,z)+(sy*sx)) = k2 * G(y,x,z) / (2*h^2);
                    
                    
                    % TERM MULT BY 2
                    % dG/dy * dUx/dy  +  G * d^2Ux/dy^2
                    M(count,itm(y,x,z))   = M(count,itm(y,x,z)) + 2 * (G(y+1,x,z)-G(y,x,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);
                    
                    M(count,itm(y+1,x,z)) = 2 * (G(y,x,z)-G(y+1,x,z))/(4*h^2) - 4 * G(y,x,z)/(h^2);

                    M(count,itm(y+2,x,z)) =  2 * G(y,x,z)/(h^2);
                    
                % We are on bottom boundary (y+1 does not exist)  -- COMPLETE  
                elseif bcs(y,x,z,1) == 1 % y == sy
                    % TERM MULT BY K2
                    % dG/dx * dUy/dy  +  G * d^2Uy/(dx*dy)
                    M(count,itm(y,x,z)+(sy*sx)) = k2 * (G(y,x+1,z)-G(y,x-1,z))/(2*h^2);

                    M(count,itm(y-1,x,z)+(sy*sx)) = k2 * (G(y,x-1,z)-G(y,x+1,z))/(2*h^2);


                    M(count,itm(y,x+1,z)+(sy*sx)) = k2 * G(y,x,z) / (2*h^2);

                    M(count,itm(y,x-1,z)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(count,itm(y-1,x+1,z)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(count,itm(y-1,x-1,z)+(sy*sx)) = k2 * G(y,x,z) / (2*h^2);
                    
                    
                    % TERM MULT BY 2
                    % dG/dy * dUx/dy  +  G * d^2Ux/dy^2
                    M(count,itm(y,x,z))   = M(count,itm(y,x,z)) + 2 * (G(y,x,z)-G(y-1,x,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                    M(count,itm(y-1,x,z)) = 2 * (G(y-1,x,z)-G(y,x,z))/(4*h^2) - 4 * G(y,x,z)/(h^2);
                    
                    M(count,itm(y-2,x,z)) = 2 * G(y,x,z)/(h^2);
                    
                else % y is on interior -- COMPLETE
                    % TERM MULT BY K2
                    % dG/dx * dUy/dy  +  G * d^2Uy/(dx*dy)
                    M(count,itm(y+1,x,z)+(sy*sx)) = k2 * (G(y,x+1,z)-G(y,x-1,z))/(4*h^2);

                    M(count,itm(y-1,x,z)+(sy*sx)) = k2 * (G(y,x-1,z)-G(y,x+1,z))/(4*h^2);


                    M(count,itm(y+1,x+1,z)+(sy*sx)) = k2 * G(y,x,z) / (4*h^2);

                    M(count,itm(y+1,x-1,z)+(sy*sx)) = -1 * k2 * G(y,x,z) / (4*h^2);

                    M(count,itm(y-1,x+1,z)+(sy*sx)) = -1 * k2 * G(y,x,z) / (4*h^2);

                    M(count,itm(y-1,x-1,z)+(sy*sx)) = k2 * G(y,x,z) / (4*h^2);
                    
                    % TERM MULT BY 2
                    % dG/dy * dUx/dy  +  G * d^2Ux/dy^2
                    M(count,itm(y+1,x,z)) = 2 * (G(y+1,x,z)-G(y-1,x,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                    M(count,itm(y,x,z))   = M(count,itm(y,x,z)) - 4 * G(y,x,z)/(h^2);

                    M(count,itm(y-1,x,z)) = 2 * (G(y-1,x,z)-G(y+1,x,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);
                    
                end
                
                % Z & X DERIVATIVES -- COMPLETE
                % We are on above boundary (z-1 does not exist) -- COMPLETE
                if bcs(y,x,z,3) == -1 %z == 1
                    % TERM MULT BY K2
                    % dG/dx * dUz/dz  +  G * d^2Uz/(dx*dz)
                    M(count,itm(y,x,z+1)+(2*sy*sx)) = k2 * (G(y,x+1,z)-G(y,x-1,z))/(2*h*dz);

                    M(count,itm(y,x,z)+(2*sy*sx)) = k2 * (G(y,x-1,z)-G(y,x+1,z))/(2*h*dz);


                    M(count,itm(y,x+1,z+1)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);

                    M(count,itm(y,x-1,z+1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(count,itm(y,x+1,z)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(count,itm(y,x-1,z)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dz * dUx/dz  +  G * d^2Ux/dz^2
                    M(count,itm(y,x,z))   = M(count,itm(y,x,z)) + 2 * (G(y,x,z)-G(y,x,z+1))/(dz^2) + 2 * G(y,x,z)/(dz^2);
                    
                    M(count,itm(y,x,z+1)) = 2 * (G(y,x,z+1)-G(y,x,z))/(dz^2) - 4 * G(y,x,z)/(dz^2);

                    M(count,itm(y,x,z+2)) =  2 * G(y,x,z)/(dz^2);
   
                % We are on below boundary (z+1 does not exist) -- COMPLETE
                elseif bcs(y,x,z,3) == 1 %z == sz 
                    % TERM MULT BY K2
                    % dG/dx * dUz/dz  +  G * d^2Uz/(dx*dz)
                    M(count,itm(y,x,z)+(2*sy*sx)) = k2 * (G(y,x+1,z)-G(y,x-1,z))/(2*h*dz);

                    M(count,itm(y,x,z-1)+(2*sy*sx)) = k2 * (G(y,x-1,z)-G(y,x+1,z))/(2*h*dz);


                    M(count,itm(y,x+1,z)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);

                    M(count,itm(y,x-1,z)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(count,itm(y,x+1,z-1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(count,itm(y,x-1,z-1)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    
                    % TERM MULT BY 2
                    % dG/dz * dUx/dz  +  G * d^2Ux/dz^2
                    M(count,itm(y,x,z))   = M(count,itm(y,x,z)) + 2 * (G(y,x,z)-G(y,x,z-1))/(dz^2) + 2 * G(y,x,z)/(dz^2);

                    M(count,itm(y,x,z-1)) = 2 * (G(y,x,z-1)-G(y,x,z))/(dz^2) - 4 * G(y,x,z)/(dz^2);
                    
                    M(count,itm(y,x,z-2)) = 2 * G(y,x,z)/(dz^2);

                else %z is on interior -- COMPLETE
                    % TERM MULT BY K2
                    % dG/dx * dUz/dz  +  G * d^2Uz/(dx*dz)
                    M(count,itm(y,x,z+1)+(2*sy*sx)) = k2 * (G(y,x+1,z)-G(y,x-1,z))/(4*h*dz);

                    M(count,itm(y,x,z-1)+(2*sy*sx)) = k2 * (G(y,x-1,z)-G(y,x+1,z))/(4*h*dz);


                    M(count,itm(y,x+1,z+1)+(2*sy*sx)) = k2 * G(y,x,z) / (4*h*dz);

                    M(count,itm(y,x-1,z+1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(count,itm(y,x+1,z-1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(count,itm(y,x-1,z-1)+(2*sy*sx)) = k2 * G(y,x,z) / (4*h*dz);
                    
                    
                    % TERM MULT BY 2
                    % dG/dz * dUx/dz  +  G * d^2Ux/dz^2
                    M(count,itm(y,x,z+1)) = 2 * (G(y,x,z+1)-G(y,x,z-1))/(dz^2) + 2 * G(y,x,z)/(dz^2);

                    M(count,itm(y,x,z))   = M(count,itm(y,x,z)) - 4 * G(y,x,z)/(dz^2);

                    M(count,itm(y,x,z-1)) = 2 * (G(y,x,z-1)-G(y,x,z+1))/(dz^2) + 2 * G(y,x,z)/(dz^2);
                    
                end
            end


            % Y DIRECTION -- COMPLETE
            if bcs(y,x,z,1) == 1 || bcs(y,x,z,1) == -1% y == 1 || y == sy
                M(county,county) = 1;
            else
                % Y ONLY DERIVATIVES -- COMPLETE
                % TERM MULT BY K1
                % dG/dy * dUy/dy  +  G * d^2Uy/dy^2
                M(county,itm(y+1,x,z)+(sy*sx)) = k1 * (G(y+1,x,z)-G(y-1,x,z))/(4*h^2) + k1 * G(y,x,z)/h^2;

                M(county,itm(y,x,z)+(sy*sx))   = -2 * k1 * G(y,x,z)/h^2;

                M(county,itm(y-1,x,z)+(sy*sx)) = k1 * (G(y-1,x,z)-G(y+1,x,z))/(4*h^2) + k1 * G(y,x,z)/h^2;
                
                % X & Y DERIVATIVES -- COMPLETE
                % We are on left boundary (x-1 does not exist)
                if bcs(y,x,z,2) == -1 % x == 1
                    % TERM MULT BY K2
                    % dG/dy * dUx/dx  +  G * d^2Ux/(dx*dy)
                    M(county,itm(y,x+1,z)) = k2 * (G(y+1,x,z)-G(y-1,x,z))/(2*h^2);

                    M(county,itm(y,x,z)) = k2 * (G(y-1,x,z)-G(y+1,x,z))/(2*h^2);


                    M(county,itm(y+1,x+1,z)) = k2 * G(y,x,z) / (2*h^2);

                    M(county,itm(y-1,x+1,z)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(county,itm(y+1,x,z)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(county,itm(y-1,x,z)) = k2 * G(y,x,z) / (2*h^2);
                    
                    % TERM MULT BY 2
                    % dG/dx * dUy/dx  +  G * d^2Uy/dx^2
                    M(county,itm(y,x,z)+(sy*sx))   = M(county,itm(y,x,z)) + 2 * (G(y,x,z)-G(y,x+1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);
                    
                    M(county,itm(y,x+1,z)+(sy*sx)) = 2 * (G(y,x+1,z)-G(y,x,z))/(4*h^2) - 4 * G(y,x,z)/(h^2);
                    
                    M(county,itm(y,x+2,z)+(sy*sx)) =  2 * G(y,x,z)/(h^2);
                    
                % We are on right boundary (x+1 does not exist) -- COMPLETE
                elseif bcs(y,x,z,2) == 1 % x == sx
                    % TERM MULT BY K2
                    % dG/dy * dUx/dx  +  G * d^2Ux/(dx*dy)
                    M(county,itm(y,x,z)) = k2 * (G(y+1,x,z)-G(y-1,x,z))/(2*h^2);

                    M(county,itm(y,x-1,z)) = k2 * (G(y-1,x,z)-G(y+1,x,z))/(2*h^2);


                    M(county,itm(y+1,x,z)) = k2 * G(y,x,z) / (2*h^2);

                    M(county,itm(y-1,x,z)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(county,itm(y+1,x-1,z)) = -1 * k2 * G(y,x,z) / (2*h^2);

                    M(county,itm(y-1,x-1,z)) = k2 * G(y,x,z) / (2*h^2);
                    
                    % TERM MULT BY 2
                    % dG/dx * dUy/dx  +  G * d^2Uy/dx^2
                    M(county,itm(y,x,z)+(sy*sx))   = M(county,itm(y,x,z)) + 2 * (G(y,x,z)-G(y,x-1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                    M(county,itm(y,x-1,z)+(sy*sx)) = 2 * (G(y,x-1,z)-G(y,x,z))/(4*h^2) - 4 * G(y,x,z)/(h^2);
                    
                    M(county,itm(y,x-2,z)+(sy*sx)) =  2 * G(y,x,z)/(h^2);
                    
                % x is on interior -- COMPLETE
                else
                    % TERM MULT BY K2
                    % dG/dy * dUx/dx  +  G * d^2Ux/(dx*dy)
                    M(county,itm(y,x+1,z)) = k2 * (G(y+1,x,z)-G(y-1,x,z))/(4*h^2);

                    M(county,itm(y,x-1,z)) = k2 * (G(y-1,x,z)-G(y+1,x,z))/(4*h^2);


                    M(county,itm(y+1,x+1,z)) = k2 * G(y,x,z) / (4*h^2);

                    M(county,itm(y-1,x+1,z)) = -1 * k2 * G(y,x,z) / (4*h^2);

                    M(county,itm(y+1,x-1,z)) = -1 * k2 * G(y,x,z) / (4*h^2);

                    M(county,itm(y-1,x-1,z)) = k2 * G(y,x,z) / (4*h^2);
                    
                    % TERM MULT BY 2
                    % dG/dx * dUy/dx  +  G * d^2Uy/dx^2
                    M(county,itm(y,x+1,z)+(sy*sx)) = 2 * (G(y,x+1,z)-G(y,x-1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                    M(county,itm(y,x,z)+(sy*sx))   = M(county,itm(y,x,z)) - 4 * G(y,x,z)/(h^2);

                    M(county,itm(y,x-1,z)+(sy*sx)) = 2 * (G(y,x-1,z)-G(y,x+1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                end
                
                % Z & Y DERIVATIVES -- COMPLETE
                % We are on above boundary (z-1 does not exist) -- COMPLETE
                if bcs(y,x,z,3) == -1 % z == 1
                    % TERM MULT BY K2
                    % dG/dy * dUz/dz  +  G * d^2Uz/(dy*dz)
                    M(county,itm(y,x,z+1)+(2*sy*sx)) = k2 * (G(y+1,x,z)-G(y-1,x,z))/(2*h*dz);

                    M(county,itm(y,x,z)+(2*sy*sx)) = k2 * (G(y-1,x,z)-G(y+1,x,z))/(2*h*dz);


                    M(county,itm(y+1,x,z+1)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);

                    M(county,itm(y-1,x,z+1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(county,itm(y+1,x,z)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(county,itm(y-1,x,z)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dz * dUy/dz  +  G * d^2Uy/dz^2
                    M(county,itm(y,x,z)+(sy*sx))   = M(county,itm(y,x,z)) + 2 * (G(y,x,z)-G(y,x,z+1))/(dz^2) + 2 * G(y,x,z)/(dz^2);
                    
                    M(county,itm(y,x,z+1)+(sy*sx)) = 2 * (G(y,x,z+1)-G(y,x,z))/(dz^2) - 4 * G(y,x,z)/(dz^2);

                    M(county,itm(y,x,z+2)+(sy*sx)) = 2 * G(y,x,z)/(dz^2);
                    
                % We are on below boundary (z+1 does not exist) -- COMPLETE
                elseif bcs(y,x,z,3) == 1 % z == sz
                    % TERM MULT BY K2
                    % dG/dy * dUz/dz  +  G * d^2Uz/(dy*dz)
                    M(county,itm(y,x,z)+(2*sy*sx)) = k2 * (G(y+1,x,z)-G(y-1,x,z))/(2*h*dz);

                    M(county,itm(y,x,z-1)+(2*sy*sx)) = k2 * (G(y-1,x,z)-G(y+1,x,z))/(2*h*dz);


                    M(county,itm(y+1,x,z)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);

                    M(county,itm(y-1,x,z)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(county,itm(y+1,x,z-1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(county,itm(y-1,x,z-1)+(2*sy*sx)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dz * dUy/dz  +  G * d^2Uy/dz^2
                    M(county,itm(y,x,z)+(sy*sx))   = M(county,itm(y,x,z)) + 2 * (G(y,x,z)-G(y,x,z-1))/(dz^2) + 2 * G(y,x,z)/(dz^2);

                    M(county,itm(y,x,z-1)+(sy*sx)) = 2 * (G(y,x,z-1)-G(y,x,z))/(dz^2) - 4 * G(y,x,z)/(dz^2);
                    
                    M(county,itm(y,x,z-2)+(sy*sx)) = 2 * G(y,x,z)/(dz^2);
                    
                % z is on interior -- COMPLETE
                else
                    % TERM MULT BY K2
                    % dG/dy * dUz/dz  +  G * d^2Uz/(dy*dz)
                    M(county,itm(y,x,z+1)+(2*sy*sx)) = k2 * (G(y+1,x,z)-G(y-1,x,z))/(4*h*dz);

                    M(county,itm(y,x,z-1)+(2*sy*sx)) = k2 * (G(y-1,x,z)-G(y+1,x,z))/(4*h*dz);


                    M(county,itm(y+1,x,z+1)+(2*sy*sx)) = k2 * G(y,x,z) / (4*h*dz);

                    M(county,itm(y-1,x,z+1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(county,itm(y+1,x,z-1)+(2*sy*sx)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(county,itm(y-1,x,z-1)+(2*sy*sx)) = k2 * G(y,x,z) / (4*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dz * dUy/dz  +  G * d^2Uy/dz^2
                    M(county,itm(y,x,z+1)+(sy*sx)) = 2 * (G(y,x,z+1)-G(y,x,z-1))/(dz^2) + 2 * G(y,x,z)/(dz^2);

                    M(county,itm(y,x,z)+(sy*sx))   = M(county,itm(y,x,z)) - 4 * G(y,x,z)/(dz^2);

                    M(county,itm(y,x,z-1)+(sy*sx)) = 2 * (G(y,x,z-1)-G(y,x,z+1))/(dz^2) + 2 * G(y,x,z)/(dz^2);
                        
                end
            end 

           
            % Z DIRECTION
            if bcs(y,x,z,3) == 1 || bcs(y,x,z,3) == -1%y == 1 || y == sy
                M(countz,countz) = 1;
            else
                % Z ONLY DERIVATIVES -- COMPLETE
                % TERM MULT BY K1
                % dG/dz * dUz/dz  +  G * d^2Uz/dz^2
                M(countz,itm(y,x,z+1)+(2*sy*sx)) = k1 * (G(y,x,z+1)-G(y,x,z-1))/(4*dz^2) + k1 * G(y,x,z)/dz^2;

                M(countz,itm(y,x,z)+(2*sy*sx))   = -2 * k1 * G(y,x,z)/dz^2;

                M(countz,itm(y,x,z-1)+(2*sy*sx)) = k1 * (G(y,x,z-1)-G(y,x,z+1))/(4*dz^2) + k1 * G(y,x,z)/dz^2;      
                        
                % X & Z DERIVATIVES -- COMPLETE
                % We are on left boundary (x-1 does not exist) -- COMPLETE
                if bcs(y,x,z,2) == -1 %x == 1
                    % TERM MULT BY K2
                    % dG/dz * dUx/dx  +  G * d^2Ux/(dx*dz)
                    M(countz,itm(y,x+1,z)) = k2 * (G(y,x,z+1)-G(y,x,z-1))/(2*h*dz);

                    M(countz,itm(y,x,z)) = k2 * (G(y,x,z-1)-G(y,x,z+1))/(2*h*dz);


                    M(countz,itm(y,x+1,z+1)) = k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y,x+1,z-1)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x,z+1)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x,z-1)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dx * dUz/dx  +  G * d^2Uz/dx^2
                    M(countz,itm(y,x,z)+(2*sy*sx))   = M(countz,itm(y,x,z)+(2*sy*sx)) + 2 * (G(y,x,z)-G(y,x+1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);
                    
                    M(countz,itm(y,x+1,z)+(2*sy*sx)) = 2 * (G(y,x+1,z)-G(y,x,z))/(4*h^2) - 4 * G(y,x,z)/(h^2);

                    M(countz,itm(y,x+2,z)+(2*sy*sx)) = 2 * G(y,x,z)/(h^2);
                        
                % We are on right boundary (x+1 does not exist) -- COMPLETE
                elseif bcs(y,x,z,2) == 1 %x == sx
                    % TERM MULT BY K2
                    % dG/dz * dUx/dx  +  G * d^2Ux/(dx*dz)
                    M(countz,itm(y,x,z)) = k2 * (G(y,x,z+1)-G(y,x,z-1))/(2*h*dz);

                    M(countz,itm(y,x-1,z)) = k2 * (G(y,x,z-1)-G(y,x,z+1))/(2*h*dz);


                    M(countz,itm(y,x,z+1)) = k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x,z-1)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x-1,z+1)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x-1,z-1)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dx * dUz/dx  +  G * d^2Uz/dx^2
                    M(countz,itm(y,x,z)+(2*sy*sx))   = M(countz,itm(y,x,z)+(2*sy*sx)) + 2 * (G(y,x,z)-G(y,x-1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                    M(countz,itm(y,x-1,z)+(2*sy*sx)) = 2 * (G(y,x-1,z)-G(y,x,z))/(4*h^2) - 4 * G(y,x,z)/(h^2);
                    
                    M(countz,itm(y,x-2,z)+(2*sy*sx)) = 2 * G(y,x,z)/(h^2);

                % x is on interior -- COMPLETE
                else
                    % TERM MULT BY K2
                    % dG/dz * dUx/dx  +  G * d^2Ux/(dx*dz)
                    M(countz,itm(y,x+1,z)) = k2 * (G(y,x,z+1)-G(y,x,z-1))/(4*h*dz);

                    M(countz,itm(y,x-1,z)) = k2 * (G(y,x,z-1)-G(y,x,z+1))/(4*h*dz);


                    M(countz,itm(y,x+1,z+1)) = k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y,x+1,z-1)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y,x-1,z+1)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y,x-1,z-1)) = k2 * G(y,x,z) / (4*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dx * dUz/dx  +  G * d^2Uz/dx^2
                    M(countz,itm(y,x+1,z)+(2*sy*sx)) = 2 * (G(y,x+1,z)-G(y,x-1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                    M(countz,itm(y,x,z)+(2*sy*sx))   = M(countz,itm(y,x,z)+(2*sy*sx)) - 4 * G(y,x,z)/(h^2);

                    M(countz,itm(y,x-1,z)+(2*sy*sx)) = 2 * (G(y,x-1,z)-G(y,x+1,z))/(4*h^2) + 2 * G(y,x,z)/(h^2);

                end
                
                % Y & Z DERIVATIVES -- COMPLETE
                % We are on top boundary (y-1 does not exist) -- COMPLETE
                if bcs(y,x,z,1) == -1 %y == 1
                    % TERM MULT BY K2
                    % dG/dz * dUy/dy  +  G * d^2Uy/(dy*dz)
                    M(countz,itm(y+1,x,z)+(sy*sx)) = k2 * (G(y,x,z+1)-G(y,x,z-1))/(2*h*dz);

                    M(countz,itm(y,x,z)+(sy*sx)) = k2 * (G(y,x,z-1)-G(y,x,z+1))/(2*h*dz);


                    M(countz,itm(y+1,x,z+1)+(sy*sx)) = k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y+1,x,z-1)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x,z+1)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x,z-1)+(sy*sx)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dy * dUz/dy  +  G * d^2Uz/dy^2
                    M(countz,itm(y,x,z)+(2*sy*sx))   = M(countz,itm(y,x,z)) + 2 * (G(y,x,z)-G(y+1,x,z))/(h^2) + 2 * G(y,x,z)/(h^2);
                    
                    M(countz,itm(y+1,x,z)+(2*sy*sx)) = 2 * (G(y+1,x,z)-G(y,x,z))/(h^2) - 4 * G(y,x,z)/h^2;

                    M(countz,itm(y+2,x,z)+(2*sy*sx)) = 2 * G(y,x,z)/h^2;
                        
                % We are on bottom boundary (y+1 does not exist) -- COMPLETE
                elseif bcs(y,x,z,1) == 1 %y == sy
                    % TERM MULT BY K2
                    % dG/dz * dUy/dy  +  G * d^2Uy/(dy*dz)
                    M(countz,itm(y,x,z)+(sy*sx)) = k2 * (G(y,x,z+1)-G(y,x,z-1))/(2*h*dz);

                    M(countz,itm(y-1,x,z)+(sy*sx)) = k2 * (G(y,x,z-1)-G(y,x,z+1))/(2*h*dz);


                    M(countz,itm(y,x,z+1)+(sy*sx)) = k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y,x,z-1)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y-1,x,z+1)+(sy*sx)) = -1 * k2 * G(y,x,z) / (2*h*dz);

                    M(countz,itm(y-1,x,z-1)+(sy*sx)) = k2 * G(y,x,z) / (2*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dy * dUz/dy  +  G * d^2Uz/dy^2
                    

                    M(countz,itm(y,x,z)+(2*sy*sx))   = M(countz,itm(y,x,z)) + 2 * (G(y,x,z)-G(y-1,x,z))/(h^2) + 2 * G(y,x,z)/(h^2);

                    M(countz,itm(y-1,x,z)+(2*sy*sx)) = 2 * (G(y-1,x,z)-G(y,x,z))/(h^2) - 4 * G(y,x,z)/h^2;
                    
                    M(countz,itm(y-2,x,z)+(2*sy*sx)) = 2 * G(y,x,z)/h^2;

                % y is on interior -- COMPLETE
                else
                    % TERM MULT BY K2
                    % dG/dz * dUy/dy  +  G * d^2Uy/(dy*dz)
                    M(countz,itm(y+1,x,z)+(sy*sx)) = k2 * (G(y,x,z+1)-G(y,x,z-1))/(4*h*dz);

                    M(countz,itm(y-1,x,z)+(sy*sx)) = k2 * (G(y,x,z-1)-G(y,x,z+1))/(4*h*dz);


                    M(countz,itm(y+1,x,z+1)+(sy*sx)) = k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y+1,x,z-1)+(sy*sx)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y-1,x,z+1)+(sy*sx)) = -1 * k2 * G(y,x,z) / (4*h*dz);

                    M(countz,itm(y-1,x,z-1)+(sy*sx)) = k2 * G(y,x,z) / (4*h*dz);
                    
                    % TERM MULT BY 2
                    % dG/dy * dUz/dy  +  G * d^2Uz/dy^2
                    M(countz,itm(y+1,x,z)+(2*sy*sx)) = 2 * (G(y+1,x,z)-G(y-1,x,z))/(h^2) + 2 * G(y,x,z)/h^2;

                    M(countz,itm(y,x,z)+(2*sy*sx))   = M(countz,itm(y,x,z)) - 4 * G(y,x,z)/(h^2);

                    M(countz,itm(y-1,x,z)+(2*sy*sx)) = 2 * (G(y-1,x,z)-G(y+1,x,z))/(h^2) + 2 * G(y,x,z)/h^2;

                end
            end
        end
    end
end

end %End of functions


