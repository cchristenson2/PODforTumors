%{ 
Forward model for the reaction diffusion equation with AC treatment in grid form
dn/dt = D*d2n/h2 + kN(1-N) - alpha1*exp(-beta1*t)*C*N - alpha2*exp(-beta2*t)*C*N

Input:
    - Cell map from start point
    - diffusivity
    - proliferation map
    - alpha1 for A
    - alpha2 for C
    - Treatment parameters
        - beta1 for A
        - beta2 for C
        - treatment duration
    - Time spacing
    - Output times
    - boundary conditions

Output:
    - Cell map at desired times
    - Full cell map time course

Contributors: Chase Christenson
%}

function [N_sim, TC] = RXDIF_3D_wAC_comb(initial, kp, d, alpha1, tx_params, t, h, dz, dt, bcs)
    
    theta = 1; %If using volume fractions

    t_ind = (t/dt) + 1; %Indices of densities to output
    nt_trx = tx_params.txduration./dt; %Indices of treatment times
    trx_on = 0; %Starts off, turns on at first treatment delivery
    trx_cnt = 1;
    
    %Intialize solution matrix
    [sy,sx,sz] = size(initial);
    nt = length(0:dt:t(end));
    Sim = zeros(sy,sx,sz,nt);
    Sim(:,:,:,1) = initial;
    
    %Time stepping
    for k = 2:nt
        temp = zeros(sy,sx,sz);
        N = Sim(:,:,:,k-1);
        
        %Get time since last treatment
        if(k-1>nt_trx(trx_cnt))
            if(trx_on==0)
                t = 0;
            end
            trx_on = 1;
            if(trx_cnt < numel(nt_trx))
                if(k-1>nt_trx(trx_cnt+1))
                    trx_cnt=trx_cnt+1;
                    t=0;
                end
                t=t+dt;
            else
                t=t+dt;
            end
        end
        
        %Space stepping
        for z = 1:sz
            for y = 1:sy
                for x = 1:sx
                    boundary = bcs(y,x,z,:);

                    %FDM in Y direction
                    if(boundary(1)==0)
                        inv_y = d*(N(y+1,x,z)-2*N(y,x,z)+N(y-1,x,z))/(h^2);

                    elseif(boundary(1)==1)
                        inv_y = d*(-2*N(y,x,z)+2*N(y-1,x,z))/(h^2);

                    elseif(boundary(1)==-1)
                        inv_y = d*(-2*N(y,x,z)+2*N(y+1,x,z))/(h^2);
                    else
                        inv_y = 0;
                    end

                    %FDM in X direction
                    if(boundary(2)==0)
                        inv_x = d*(N(y,x+1,z)-2*N(y,x,z)+N(y,x-1,z))/(h^2);

                    elseif(boundary(2)==1)
                        inv_x = d*(-2*N(y,x,z)+2*N(y,x-1,z))/(h^2);

                    elseif(boundary(2)==-1)
                        inv_x = d*(-2*N(y,x,z)+2*N(y,x+1,z))/(h^2);
                    else
                        inv_x = 0;
                    end
                    
                    %FDM in Z direction
                    if(boundary(3)==0)
                        inv_z = d*(N(y,x,z+1)-2*N(y,x,z)+N(y,x,z-1))/(dz^2);

                    elseif(boundary(3)==1)
                        inv_z = d*(-2*N(y,x,z)+2*N(y,x,z-1))/(dz^2);

                    elseif(boundary(3)==-1)
                        inv_z = d*(-2*N(y,x,z)+2*N(y,x,z+1))/(dz^2);
                    else
                        inv_z = 0;
                    end

                    invasion = inv_y + inv_x + inv_z;
                    prolif   = N(y,x,z)*kp(y,x,z)*(1-(N(y,x,z)/theta));

                    %Treatment calculation
                    if(trx_on~=0)
                        treat = alpha1*(exp(-tx_params.beta1*t)+exp(-tx_params.beta2*t))*tx_params.C(y,x,z)*N(y,x,z);
                    else
                        treat = 0;
                    end

                    temp(y,x,z) = N(y,x,z) + dt*(invasion + prolif - treat);
                end
            end
        end
        
        Sim(:,:,:,k) = temp;
    end
    TC = squeeze(sum(sum(Sim,2),1));
%     TC = Sim;
    N_sim = Sim(:,:,:,t_ind);
end