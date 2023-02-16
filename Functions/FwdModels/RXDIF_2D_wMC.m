%{ 
Forward model for the reaction diffusion equation with AC treatment in grid form
dn/dt = D*d2n/h2 + kN(1-N)

Input:
    - Cell map from start point
    - diffusivity
    - proliferation map
    - Time spacing
    - Output times
    - boundary conditions

Output:
    - Cell map at desired times
    - Full cell map time course

Contributors: Chase Christenson
%}

function [N_sim, TC] = RXDIF_2D_wMC(initial, kp, d0, t, h, dt, bcs, M, E, nu, matX, matY)
    
    theta = 1; %If using volume fractions
    freq = 25;
    
    t_ = (t./dt) + 1; %Indices of densities to output
    
    %Intialize solution matrix
    [sy,sx] = size(initial);
    nt = length(0:dt:t(end));
    Sim = zeros(sy,sx,nt);
    Sim(:,:,1) = initial;
    
    
    
    % Initialize damped diffusivity
    damper = get_damper(matX, matY, initial, M, E, nu);
    d = reshape(d0.*damper, sy, sx);
    
    %Time stepping
    for k = 1:nt-1
        temp = zeros(sy,sx);
        N = Sim(:,:,k);
        
        %Space stepping
        for y = 1:sy
            for x = 1:sx
                boundary = bcs(y,x,:);

                %FDM in Y direction
                if(boundary(1)==0)
                    y1 = (d(y+1,x)-d(y-1,x))/(2*h);
                    y2 = (N(y+1,x)-N(y-1,x))/(2*h);
                    lap_y = d(y,x)*(N(y+1,x)-2*N(y,x)+N(y-1,x))/(h^2);
                    inv_y = lap_y + (y1*y2);

                elseif(boundary(1)==1)
                    inv_y = d(y,x)*(-2*N(y,x)+2*N(y-1,x))/(h^2);

                elseif(boundary(1)==-1)
                    inv_y = d(y,x)*(-2*N(y,x)+2*N(y+1,x))/(h^2);
                else
                    inv_y = 0;
                end

                %FDM in X direction
                if(boundary(2)==0)
                    x1 = (d(y,x+1)-d(y,x-1))/(2*h);
                    x2 = (N(y,x+1)-N(y,x-1))/(2*h);
                    lap_x = d(y,x)*(N(y,x+1)-2*N(y,x)+N(y,x-1))/(h^2);
                    inv_x = lap_x + (x1*x2);

                elseif(boundary(2)==1)
                    inv_x = d(y,x)*(-2*N(y,x)+2*N(y,x-1))/(h^2);

                elseif(boundary(2)==-1)
                    inv_x = d(y,x)*(-2*N(y,x)+2*N(y,x+1))/(h^2);
                else
                    inv_x = 0;
                end

                invasion = inv_y + inv_x;
                prolif   = N(y,x)*kp(y,x)*(1-(N(y,x)/theta));

                temp(y,x) = N(y,x) + dt*(invasion + prolif);
            end
        end
        
        Sim(:,:,k+1) = temp;
        
        if mod(k, freq) == 0
            damper = get_damper(matX, matY, temp, M, E, nu);
            d = reshape(d0.*damper, sy, sx);
        end
        
        
    end
    TC = squeeze(sum(sum(Sim,2),1));
    N_sim = Sim(:,:,t_);
end