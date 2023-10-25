%{
Mechanics matrix build for 2D domains

inputs:
    - h: in-plane spacing
    - tissues: map of tissue segmentations, 1 = adipose, 2 = fibroglandular
        0 = tumor (sy,sx)
    - bcs: Boundary condition map (sy,sx,2)
outputs:
    - M: matrix for mechanical displacement solve
    - E: Map of elastic modulus from literature
    - nu: Poisson's ratio

Contributors: Chase Christenson, Casey Stowers

Based off reference:
Hormuth, David A. et al. 
“Mechanically Coupled Reaction-Diffusion Model to Predict Glioma Growth: 
Methodological Details.” Cancer Systems Biology 1711 (2018): 225–241.


%}

function [M, E, nu] =  mech_matrix_build_2D(h, tissues, bcs)
% h is our step in x and y ( h= h)
% tissue is our tissue type and tumor map for this patient
% bcs is the map of boundary conditions 
    %  Boundary matrix bcs(y,x,2) = [y,x, [y boundary type, x boundary type]]
    %             Value = -1 if boundary is behind, left
    %             Value = 0 if no boundary is present
    %             Value = 1 if boundary is in front, right
    %             Value = 2 if outside of mask

% a bit of setup
[sy,sx] = size(tissues);
M = sparse(sy*sx*2,sy*sx*2); % it big

% First thing, use our tissue map to set mechanics params
% we use one nu value, but multiple E values
nu = 0.45; 
E_adipose = 2e3;
E_fibro = 2*2e3;
E_tumor = 10*2e3;

E = zeros(size(tissues));
E(tissues == 1) = E_adipose;
E(tissues == 2) = E_fibro;
E(tissues == 0) = E_tumor;

% E = imgaussfilt(E,2);

E(bcs(:,:,1) ==2 ) = 0;



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
for x = 1:sx
    for y = 1:sy
        count = count + 1;
        itm(y,x) = count;
    end
end

% set up mult factors for BCs

for y = 1:sy
    for x = 1:sx
        count = itm(y,x);
        county = itm(y,x)+sy*sx;
        %Ke = Ke_mat(y,x);
        % ARE WE IN THE BREAST
        if bcs(y,x,2) == 2
            M(count,count) = 1; % forces 0 displacement
            M(county,county) = 1; % forces 0 displacement
        else

            % X DIRECTION
            % First, check to see if we are on an x boundary
            % If we are, set to 1
            if bcs(y,x,2) == 1 || bcs(y,x,2) == -1 %x == 1 || x == sx
                M(count,count) = 1; % forces 0 displacement
            else
                % Check to see if we are on a top boundary
                if bcs(y,x,1) == -1 %y == 1
                    % Can't do y-1 if we're on top now can we

                   % TERM MULT BY k1
                   % no y derivatives!! so change nothing :)
                   % dG/dx * dUx/dx + G d^2Ux/dx^2
                   % Ux(y,x+1)
                        M(count,itm(y,x+1)) = k1 * ((G(y,x+1)-G(y,x-1))/(4*h^2) + G(y,x) / h^2);
                   % Ux(y,x)
                        M(count,itm(y,x)) = k1 * ( -2 * G(y,x) / h^2);
                   % Ux(y,x-1)
                        M(count,itm(y,x-1)) = k1 * (-1 * (G(y,x+1)-G(y,x-1))/(4*h^2) + G(y,x) / h^2);

                    % TERM MULT BY 2
                    % dG/dy * dUx/dy + G d^2Ux/dy^2
                    % gotta change like everything D':
                    try
                    % Ux(y,x)
                        M(count,itm(y,x)) = 2 * (-1*(G(y+1,x) - G(y,x))/(2*h^2) + G(y,x) / h^2) + M(count,itm(y,x));
                    % Ux(y+1,x)
                        M(count,itm(y+1,x)) = 2 * ((G(y+1,x) - G(y,x))/(2*h^2) -2*G(y,x) / h^2);
                    % Ux(y+2,x)
                        M(count,itm(y+2,x)) = 2 * ( G(y,x) / h^2);
                    catch
                    % Ux(y,x)
                        M(count,itm(y,x)) = 2 * (-1*(G(y+1,x) - G(y,x))/(2*h^2) + G(y,x) / h^2) + M(count,itm(y,x));
                    % Ux(y+1,x)
                        M(count,itm(y+1,x)) = 2 * ((G(y+1,x) - G(y,x))/(2*h^2) - G(y,x) / h^2);
                    end

                    % TERM MULT BY k2 - first derivatives
                    % dG/dx * dUy/dy
                    % need to remove y-1 on Uy
                    % Uy(y+1,x)
                        M(count,sx*sy+itm(y+1,x)) = k2 * (G(y,x+1) - G(y,x-1)) / (2*h^2);
                    % Uy(y,x)
                        M(count,sx*sy+itm(y,x)) = k2 * -1*(G(y,x+1) - G(y,x-1)) / (2*h^2);

                    % TERM MULT BY k2 - mixed derivative
                    % again remove y-1s
                    % G*d^2Uy/dxdy
                    % Uy(y+1,x+1)
                        M(count,sx*sy+itm(y+1,x+1)) = k2 * G(y,x) / (2*h^2);
                    % Uy(y,x+1)
                        M(count,sx*sy+itm(y,x+1)) = -k2 * G(y,x) / (2*h^2);
                    % Uy(y+1,x-1)
                        M(count,sx*sy+itm(y+1,x-1)) = -k2 * G(y,x) / (2*h^2);
                    % Uy(y,x-1) 
                        M(count,sx*sy+itm(y,x-1)) = k2 * G(y,x) / (2*h^2); 

                % Check to see if we are on a bottom boundary
                elseif bcs(y,x,1) == 1%y == sy
                   % Means no y+1, need to change to a one sided differencing

                   % TERM MULT BY k1
                   % No changes b/c no y derivatives
                   % dG/dx * dUx/dx + G d^2Ux/dx^2
                   % Ux(y,x+1)
                        M(count,itm(y,x+1)) = k1 * ((G(y,x+1)-G(y,x-1))/(4*h^2) + G(y,x) / h^2);
                   % Ux(y,x)
                        M(count,itm(y,x)) = k1 * ( -2 * G(y,x) / h^2);
                   % Ux(y,x-1)
                        M(count,itm(y,x-1)) = k1 * (-1 * (G(y,x+1)-G(y,x-1))/(4*h^2) + G(y,x) / h^2);

                    % TERM MULT BY 2
                    % dG/dy * dUx/dy + G d^2Ux/dy^2
                    % change basically everything just rip
                    
                    
                    try
                    % Ux(y,x) 
                        M(count,itm(y,x)) = 2 * ((G(y,x) - G(y-1,x))/(2*h^2) + G(y,x) / h^2) + M(count,itm(y,x));

                    % Ux(y-1, x)
                        M(count, itm(y-1,x)) = 2 * (-1*(G(y,x) - G(y-1,x))/(2*h^2) -2 * G(y,x) / h^2);

                    % Ux(y-2, x)
                        M(count,itm(y-2,x)) = 2 * ( G(y,x) / h^2);
                    catch
                    % Ux(y,x) 
                        M(count,itm(y,x)) = 2 * ((G(y,x) - G(y-1,x))/(2*h^2) + G(y,x) / h^2) + M(count,itm(y,x));

                    % Ux(y-1, x)
                        M(count, itm(y-1,x)) = 2 * (-1*(G(y,x) - G(y-1,x))/(2*h^2) - G(y,x) / h^2); 
                    end


                    % TERM MULT BY k2 - first derivatives
                    % dG/dx * dUy/dy
                    % change for Uy
                    % Uy(y,x)
                        M(count,sx*sy+itm(y,x)) = k2 * (G(y,x+1) - G(y,x-1)) / (2*h^2);
                    % Uy(y-1,x)
                        M(count,sx*sy+itm(y-1,x)) = k2 * -1*(G(y,x+1) - G(y,x-1)) / (2*h^2);

                    % TERM MULT BY k2 - mixed derivative
                    % remove all our y+1s
                    % G*d^2Uy/dxdy
                    % Uy(y,x+1)
                        M(count,sx*sy+itm(y,x+1)) = k2 * G(y,x) / (2*h^2);
                    % Uy(y-1,x+1)
                        M(count,sx*sy+itm(y-1,x+1)) = -k2 * G(y,x) / (2*h^2);
                    % Uy(y,x-1)
                        M(count,sx*sy+itm(y,x-1)) = -k2 * G(y,x) / (2*h^2);
                    % Uy(y-1,x-1) 
                        M(count,sx*sy+itm(y-1,x-1)) = k2 * G(y,x) / (2*h^2); 

                % Otherwise, we are on the interior    
                else
                   % TERM MULT BY k1
                   % dG/dx * dUx/dx + G d^2Ux/dx^2
                   % Ux(y,x+1)
                        M(count,itm(y,x+1)) = k1 * ((G(y,x+1)-G(y,x-1))/(4*h^2) + G(y,x) / h^2);
                   % Ux(y,x)
                        M(count,itm(y,x)) = k1 * ( -2 * G(y,x) / h^2);
                   % Ux(y,x-1)
                        M(count,itm(y,x-1)) = k1 * (-1 * (G(y,x+1)-G(y,x-1))/(4*h^2) + G(y,x) / h^2);

                    % TERM MULT BY 2
                    % dG/dy * dUx/dy + G d^2Ux/dy^2
                    % Ux(y+1,x)
                        M(count,itm(y+1,x)) = 2 * ((G(y+1,x) - G(y-1,x))/(4*h^2) + G(y,x) / h^2);
                    % Ux(y,x)
                        M(count,itm(y,x)) = 2 * ( -2 * G(y,x) / h^2) + M(count,itm(y,x));
                    % Ux(y-1,x)
                        M(count,itm(y-1,x)) = 2 * (-1*(G(y+1,x) - G(y-1,x))/(4*h^2) + G(y,x) / h^2);

                    % TERM MULT BY k2 - first derivatives
                    % dG/dx * dUy/dy
                    % Uy(y+1,x)
                        M(count,sx*sy+itm(y+1,x)) = k2 * (G(y,x+1) - G(y,x-1)) / (4*h^2);
                    % Uy(y-1,x)
                        M(count,sx*sy+itm(y-1,x)) = k2 * -1*(G(y,x+1) - G(y,x-1)) / (4*h^2);

                    % TERM MULT BY k2 - mixed derivative
                    % G*d^2Uy/dxdy
                    % Uy(y+1,x+1)
                        M(count,sx*sy+itm(y+1,x+1)) = k2 * G(y,x) / (4*h^2);
                    % Uy(y-1,x+1)
                        M(count,sx*sy+itm(y-1,x+1)) = -k2 * G(y,x) / (4*h^2);
                    % Uy(y+1,x-1)
                        M(count,sx*sy+itm(y+1,x-1)) = -k2 * G(y,x) / (4*h^2);
                    % Uy(y-1,x-1) 
                        M(count,sx*sy+itm(y-1,x-1)) = k2 * G(y,x) / (4*h^2); 
                end

            end


            % Y DIRECTION
            if bcs(y,x,1) == 1 || bcs(y,x,1) == -1%y == 1 || y == sy
                M(county,county) = 1;
            else
                % Check and see if we are at left bound
               if bcs(y,x,2) == -1%x == 1
                   % If we're on a left bound, we can't do x-1
                   % So switch to forwards differencing in those cases

                   % TERM MULT BY k1 - only x derivatives, so no changes
                   % (dg/dy * dUy/dy + G d2Uy/dy^2)
                   %Uy(y+1,x)
                        M(county,sx*sy+itm(y+1,x)) = k1 * ((G(y+1,x)-G(y-1,x))/(4*h^2) + G(y,x) / h^2); 
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = k1 * (-2 * G(y,x) / h^2);
                   % Uy(y-1,x)
                        M(county,sx*sy+itm(y-1,x)) = k1 * (-1*(G(y+1,x)-G(y-1,x))/(4*h^2) + G(y,x) / h^2);

                   % FIRST DERIVATIVES IN k2 TERM - need to changes Ux to
                   % forwards
                   % dg/dy * dUx/dx
                   % Ux(y,x+1)
                        M(county,itm(y,x+1)) = k2 * (G(y+1,x)-G(y-1,x)) / (2*h^2);
                   % Ux(y,x)
                        M(county,itm(y,x)) = k2 * -1 * (G(y+1,x)-G(y-1,x)) / (2*h^2);

                   % FIRST DERIVATIVES IN 2 TERM - need to change G to forwards
                   % dG/dx * dUx/dy
                   % Ux(y+1,x)
                        M(county,itm(y+1,x)) = 2 * (G(y,x+1)-G(y,x)) / (2*h^2);
                   % Ux(y-1,x)
                        M(county,itm(y-1,x)) = 2 * -1 * (G(y,x+1)-G(y,x)) / (2*h^2);

                   % MIXED DERIVATIVE TERMS - change to forwards in x
                   % There is a mixed derivative on Ux with both k2 and 2 
                   k3 = (k2 + 2) * G(y,x);

                   % G* d^2Ux/dydx
                   % Ux(y+1,x+1)
                        M(county,itm(y+1,x+1)) = k3 / (2*h^2);
                   % Ux(y-1,x+1)
                        M(county,itm(y-1,x+1)) = k3 * -1 / (2*h^2);
                   % Ux(y+1,x)
                        M(county,itm(y+1,x)) = k3 * -1 / (2*h^2) + M(county,itm(y+1,x));
                   % Ux(y-1,x)
                        M(county,itm(y-1,x)) = k3 / (2*h^2) + M(county,itm(y-1,x));


               % or are we on a right bound
               elseif bcs(y,x,2) == 1
                   % if we are on a right bound, then we can't do x+1
                   % So, we want to switch to backwards differencing on the x
                   % derivatives

                   % TERM MULT BY k1 - y derivatives only, so don't need to
                   % make any changes
                   % (dg/dy * dUy/dy + G d2Uy/dy^2)
                   % Uy(y+1,x)
                        M(county,sx*sy+itm(y+1,x)) = k1 * ((G(y+1,x)-G(y-1,x))/(4*h^2) + G(y,x) / h^2); 
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = k1 * (-2 * G(y,x) / h^2);
                   % Uy(y-1,x)
                        M(county,sx*sy+itm(y-1,x)) = k1 * (-1*(G(y+1,x)-G(y-1,x))/(4*h^2) + G(y,x) / h^2);

                   % FIRST DERIVATIVES IN k2 TERM - need to change to backwards
                   % on Ux
                   % dg/dy * dUx/dx
                   % Ux(y,x)
                        M(county,itm(y,x)) = k2 * (G(y+1,x)-G(y-1,x)) / (2*h^2);
                   % Ux(y,x-1)
                        M(county,itm(y,x-1)) = k2 * -1 * (G(y+1,x)-G(y-1,x)) / (2*h^2);

                   % FIRST DERIVATIVES IN 2 TERM - need to change to backwards
                   % on G
                   % dG/dx * dUx/dy
                   % Ux(y+1,x)
                        M(county,itm(y+1,x)) = 2 * (G(y,x)-G(y,x-1)) / (2*h^2);
                   % Ux(y-1,x)
                        M(county,itm(y-1,x)) = 2 * -1 * (G(y,x)-G(y,x-1)) / (2*h^2);

                   % MIXED DERIVATIVE TERMS - change to backwards
                   % There is a mixed derivative on Ux with both k2 and 2 
                   k3 = (k2 + 2) * G(y,x);

                   % G* d^2Ux/dydx
                   % Ux(y+1,x)
                        M(county,itm(y+1,x)) = k3 / (2*h^2) + M(county,itm(y+1,x));
                   % Ux(y-1,x)
                        M(county,itm(y-1,x)) = k3 * -1 / (2*h^2) + M(county,itm(y-1,x));
                   % Ux(y+1,x-1)
                        M(county,itm(y+1,x-1)) = k3 * -1 / (2*h^2);
                   % Ux(y+1,x+1)
                        M(county,itm(y-1,x-1)) = k3 / (2*h^2);

               % we are in the interior   
               else
                   % TERM MULT BY k1 
                   % (dg/dy * dUy/dy + G d2Uy/dy^2)
                   %Uy(y+1,x)
                        M(county,sx*sy+itm(y+1,x)) = k1 * ((G(y+1,x)-G(y-1,x))/(4*h^2) + G(y,x) / h^2); 
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = k1 * (-2 * G(y,x) / h^2);
                   % Uy(y-1,x)
                        M(county,sx*sy+itm(y-1,x)) = k1 * (-1*(G(y+1,x)-G(y-1,x))/(4*h^2) + G(y,x) / h^2);

                   % FIRST DERIVATIVES IN k2 TERM
                   % dg/dy * dUx/dx
                   % Ux(y,x+1)
                        M(county,itm(y,x+1)) = k2 * (G(y+1,x)-G(y-1,x)) / (4*h^2);
                   % Ux(y,x-1)
                        M(county,itm(y,x-1)) = k2 * -1 * (G(y+1,x)-G(y-1,x)) / (4*h^2);

                   % FIRST DERIVATIVES IN 2 TERM
                   % dG/dx * dUx/dy
                   % Ux(y+1,x)
                        M(county,itm(y+1,x)) = 2 * (G(y,x+1)-G(y,x-1)) / (4*h^2);
                   % Ux(y-1,x)
                        M(county,itm(y-1,x)) = 2 * -1 * (G(y,x+1)-G(y,x-1)) / (4*h^2);

                   % MIXED DERIVATIVE TERMS
                   % There is a mixed derivative on Ux with both k2 and 2 
                   k3 = (k2 + 2) * G(y,x);

                   % G* d^2Ux/dydx
                   % Ux(y+1,x+1)
                        M(county,itm(y+1,x+1)) = k3 / (4*h^2);
                   % Ux(y-1,x+1)
                        M(county,itm(y-1,x+1)) = k3 * -1 / (4*h^2);
                   % Ux(y+1,x-1)
                        M(county,itm(y+1,x-1)) = k3 * -1 / (4*h^2);
                   % Ux(y+1,x+1)
                        M(county,itm(y-1,x-1)) = k3 / (4*h^2);

                end
            end
        end

    end
end



end