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

function [N_sim, TC] = RXDIF_2D(initial, kp, d, t, h, dt, bcs)
    
    theta = 1; %If using volume fractions

    t_ = (t/dt) + 1; %Indices of densities to output
    
    %Intialize solution matrix
    [sy,sx] = size(initial);
    nt = length(0:dt:t(end));
    Sim = zeros(sy,sx,nt);
    Sim(:,:,1) = initial;
    
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
                    inv_y = d*(N(y+1,x)-2*N(y,x)+N(y-1,x))/(h^2);

                elseif(boundary(1)==1)
                    inv_y = d*(-2*N(y,x)+2*N(y-1,x))/(h^2);

                elseif(boundary(1)==-1)
                    inv_y = d*(-2*N(y,x)+2*N(y+1,x))/(h^2);
                else
                    inv_y = 0;
                end

                %FDM in X direction
                if(boundary(2)==0)
                    inv_x = d*(N(y,x+1)-2*N(y,x)+N(y,x-1))/(h^2);

                elseif(boundary(2)==1)
                    inv_x = d*(-2*N(y,x)+2*N(y,x-1))/(h^2);

                elseif(boundary(2)==-1)
                    inv_x = d*(-2*N(y,x)+2*N(y,x+1))/(h^2);
                else
                    inv_x = 0;
                end

                invasion = inv_y + inv_x;
                prolif   = N(y,x)*kp(y,x)*(1-(N(y,x)/theta));

                temp(y,x) = N(y,x) + dt*(invasion + prolif);
            end
        end
        
        Sim(:,:,k+1) = temp;
    end
    TC = squeeze(sum(sum(Sim,2),1));
    N_sim = Sim(:,:,t_);
end