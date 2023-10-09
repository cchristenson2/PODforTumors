%{ 
Forward model for the reaction diffusion equation with AC treatment in grid form
dn/dt = D*d^2n/h^2 + kp*N(1-N) - alpha1*(exp(-beta1*t) - exp(-beta2*t))*C*N

Input:
    - N0; Cell map from start point
    - kp; proliferation map
    - d; diffusivity
    - alpha; efficacy for A/C
    - tx_params; Treatment parameters
        - beta1 for A
        - beta2 for C
        - txduration; treatment duration
    - Time spacing
    - Output times
    - boundary conditions

Output:
    - Cell map at desired times
    - Full cell map time course

Contributors: Chase Christenson
%}

function [N_sim, TC] = RXDIF_2D_wAC_comb(initial, kp, d, alpha1, tx_params, t, h, dt, bcs)
    
    theta = 1; %If using volume fractions
    thresh = 1e-3;
    
    t_ind = (t/dt) + 1; %Indices of densities to output
    nt_trx = tx_params.txduration./dt; %Indices of treatment times
    trx_on = 0; %Starts off, turns on at first treatment delivery
    trx_cnt = 1;
    
    %Intialize solution matrix
    [sy,sx] = size(initial);
    nt = length(0:dt:t(end));
    Sim = zeros(sy,sx,nt);
    Sim(:,:,1) = initial;
    
    %Time stepping
    for k = 2:nt
        temp = zeros(sy,sx);
        N = Sim(:,:,k-1);
        
        %Get time since last treatment
        if(k-1>nt_trx(trx_cnt))
            if(trx_on==0)
                t = 0;
                trx_on = 1;
            else
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
        end
        
        %Space stepping
        for y = 1:sy
            for x = 1:sx
                boundary = bcs(y,x,:);
                
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

                %Treatment calculation
                if(trx_on~=0)
                    treat = alpha1*(exp(-tx_params.beta1*t)+exp(-tx_params.beta2*t))*tx_params.C(y,x)*N(y,x);
                else
                    treat = 0;
                end

                temp(y,x) = N(y,x) + dt*(invasion + prolif - treat);
            end
        end
        
        Sim(:,:,k) = temp;
    end
    TC = squeeze(sum(sum(Sim,2),1));
%     TC = Sim;
    N_sim = Sim(:,:,t_ind);
end