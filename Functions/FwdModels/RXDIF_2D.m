%{ 
Forward model for the reaction diffusion equation with AC treatment in grid form
dn/dt = D*d^2n/h^2 + kp*N(1-N)

Input:
    - initial; Cell map from start point
    - kp; proliferation map
    - d; diffusivity value
    - t; Output times
    - h; in-plane spacing
    - dt; Time spacing
    - bcs; boundary conditions

Output:
    - N_sim; Cell map at desired times
    - TC; Full cell map time course
    - kp; final proliferation map after accounting for movement into new voxels

Contributors: Chase Christenson
%}

function [N_sim, TC, kp] = RXDIF_2D(initial, kp, d, t, h, dt, bcs)
    
    theta = 1; %If using volume fractions
    thresh = 1e-3;

    t_ = (t./dt) + 1; %Indices of densities to output
    
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
                
                if(kp(y,x)==0 && N(y,x) >= thresh)
                    try n1 = kp(y+1,x-1); catch, n1 = 0; end
                    try n2 = kp(y+1,x);   catch, n2 = 0; end
                    try n3 = kp(y+1,x+1); catch, n3 = 0; end
                    try n4 = kp(y,x+1);   catch, n4 = 0; end
                    try n5 = kp(y-1,x+1); catch, n5 = 0; end
                    try n6 = kp(y-1,x);   catch, n6 = 0; end
                    try n7 = kp(y-1,x-1); catch, n7 = 0; end
                    try n8 = kp(y,x-1);   catch, n8 = 0; end
                    
                    nn = [n1, n2, n3, n4, n5, n6, n7, n8];
                    nn = nn(nn~=0);
                    kp(y,x) = sum(nn)/numel(nn);
                    if(isnan(kp(y,x)))
                        kp(y,x) = 0;
                    end
                end
                
                boundary = bcs(y,x,:);

                %FDM in Y direction
                if(boundary(1)==0)
                    lap_y = d*(N(y+1,x)-2*N(y,x)+N(y-1,x))/(h^2);
                    inv_y = lap_y;

                elseif(boundary(1)==1)
                    inv_y = d*(-2*N(y,x)+2*N(y-1,x))/(h^2);

                elseif(boundary(1)==-1)
                    inv_y = d*(-2*N(y,x)+2*N(y+1,x))/(h^2);
                else
                    inv_y = 0;
                end

                %FDM in X direction
                if(boundary(2)==0)
                    lap_x = d*(N(y,x+1)-2*N(y,x)+N(y,x-1))/(h^2);
                    inv_x = lap_x;

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