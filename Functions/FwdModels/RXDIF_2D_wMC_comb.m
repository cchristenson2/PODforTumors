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

function [N_sim, TC] = RXDIF_2D_wMC_comb(initial, kp, d0, alpha1, tx_params, t, h, dt, bcs, M, E, nu, matX, matY)
    
    theta = 1; %If using volume fractions
    freq = 25;

    t_ = (t./dt) + 1; %Indices of densities to output
    
    nt_trx = tx_params.txduration./dt; %Indices of treatment times
    trx_on = 0; %Starts off, turns on at first treatment delivery
    trx_cnt = 1;
    
    %Intialize solution matrix
    [sy,sx] = size(initial);
    nt = length(0:dt:t(end));
    Sim = zeros(sy,sx,nt);
    Sim(:,:,1) = initial;
    
    
    
    % Initialize damped diffusivity
    damper = get_damper(matX, matY, initial, M, E, nu);
    d = reshape(d0.*damper,sy,sx);
    
    %Time stepping
    for k = 1:nt-1
        temp = zeros(sy,sx);
        N = Sim(:,:,k);
        
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
                
                %Treatment calculation
                if(trx_on~=0)
%                     treat = (alpha1*exp(-tx_params.beta1*t) + alpha2*exp(-tx_params.beta2*t))*tx_params.C(y,x)*N(y,x);
                    treat = alpha1*(exp(-tx_params.beta1*t) + exp(-tx_params.beta2*t))*tx_params.C(y,x)*N(y,x);
                else
                    treat = 0;
                end

                invasion = inv_y + inv_x;
                prolif   = N(y,x)*kp(y,x)*(1-(N(y,x)/theta));

                temp(y,x) = N(y,x) + dt*(invasion + prolif - treat);
            end
        end
        
        Sim(:,:,k+1) = temp;
        
        if mod(k, freq) == 0
            damper = get_damper(matX, matY, N, M, E, nu);
            d = reshape(d0 .* damper, sy,sx);
        end
        
        
    end
    TC = squeeze(sum(sum(Sim,2),1));
    N_sim = Sim(:,:,t_);
end