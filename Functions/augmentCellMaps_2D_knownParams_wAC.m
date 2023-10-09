%{ 
Augmentation for data between measured time points

Inputs:
    - N; cell maps we want to augment between
    - t; time of image acquisition in N
    - h; grid spacing for the data
    - dt; time step we want to use
    - kp; true proliferation if available
    - d; true diffusivity if available
    - bcs; boundary conditions for forward runs

Outputs:
    - N_aug; augmented cell data, with original data at the correct location

Contributors: Chase Christenson
%}


function N_aug = augmentCellMaps_2D_knownParams_wAC(N, t, h, dt, kp, d, alpha, tx_params, bcs)

    %How many gaps do we need to fill in?
    [sy,sx,nt] = size(N);
    gaps = nt-1;
    
    N_aug = [];
    
    %For each of the gaps we want to fill
    for i = 1:gaps
        %Get measured maps on the edge of the gap
        N0 = N(:,:,i);
        
        %Get time to run for current gap
        t_end = t(i+1) - t(i);
        t_out = 0:dt:t_end - dt; %We do not last image from simulation
        
        temp = RXDIF_2D_wAC_comb(N0,kp.*ones(sy,sx),d,alpha,tx_params,t_out,h,dt,bcs);
        
        N_aug = cat(3, N_aug, temp);
    end
    N_aug = cat(3, N_aug, N(:,:,end));
end

