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


function N_aug = augmentCellMaps_2D_knownParams(N, t, h, dt, kp, d, bcs)

    %How many gaps do we need to fill in?
    [sy,sx,nt] = size(N);
    gaps = nt-1;
    
    N_aug = N(:,:,1); %Initialized with original map
    
    %For each of the gaps we want to fill
    for i = 1:gaps
        %Get measured maps on the edge of the gap
        N1 = N(:,:,i);
        N2 = N(:,:,i+1);
        
        %Get time to run for current gap
        t_end = t(i+1) - t(i);
        t_out = dt:dt:t_end - dt; %We do not need first or last image from simulation
        
        %Estimate parameters for forward runs
        %Proliferation from exponential growth
        if(isempty(kp))
            kp = log(sum(N2,'all')/sum(N1,'all'))/t_end;
        end
        %Diffusivity from simplified fick's law
        %Change in radius squared over time
        if(isempty(d))
            d = ((h^2 * numel(N2(N2>0.10))) - (h^2 * numel(N1(N1>0.10))))^2/(4*t_end);
        end
        
        temp = RXDIF_2D(N1,kp.*ones(sy,sx),d,t_out,h,dt,bcs);
        
        N_aug = cat(3, N_aug, temp);
    end

end

