function [N_aug,kp] = Augment_wMC_wTRX_comb_3D(N, t, h, dz, dt, bcs, bounds, ntp_cal, tx_params, M, E, nu, matX, matY, matZ)
    thresh = 0.15;

    [sy,sx,sz,~] = size(N);
    N_aug = zeros(sy,sx,sz,t(end)/dt+1);
    N_aug(:,:,:,1) = N(:,:,:,1);
    
    %% Loop once for 2 time point calibrations, twice for 3
    for i = 1:ntp_cal
    %% Get parameters for forward run
        N1 = N(:,:,:,i);
        N2 = N(:,:,:,i+1);
        
        if(i==1)
            t1 = 0;
        else
            t1 = t(i-1);
        end
        t2 = t(i);
    
        %Proliferation - estimation of change in cells based on logisitc growth equation
        kp = -1*log(N2./N1)./(t2-t1); kp(isnan(kp)) = 0; kp(isinf(kp)) = 0;
        
        idx = find(kp);
        kp_change = normalize(kp(idx), 'range', [bounds.kp_bounds(1), bounds.kp_bounds(end)]);
        kp = zeros(size(kp));
        kp(idx) = kp_change;
        
        %Diffusivity - estimation of change in radius over time (area of circle assumption)
        V1 = numel(N1(N1>thresh))*h^2*dz;
        V2 = numel(N2(N2>thresh))*h^2*dz;
        
        %Alpha parameters
        alpha1 = (bounds.alpha_bounds(1) + bounds.alpha_bounds(end))/2;
        
        del_r = (3*V2/(4*pi))^(1/3) - (3*V1/(4*pi))^(1/3);
        d = del_r^2 / (6*(t2-t1));

        if(d<bounds.d_bounds(1))
            d = bounds.d_bounds(1);
        elseif(d>bounds.d_bounds(end))
            d = bounds.d_bounds(end);
        end
        
        if(i==1)
            kp_out = kp;
        end
        
        
    %% Run forward eval
        if(ntp_cal == 1) % two time point calibration, run twice with same kp from the two maps
            for j = 1:2
                if(j==1)
                    initial = N(:,:,:,1);
                    t_out = dt:dt:t(1);
                    t_idx = t_out./dt + 1;
                else
                    initial = N(:,:,:,2);
                    t_out = [t(1)+dt:dt:t(2)]-t(1);
                    t_idx = (t_out + t(1))./dt + 1;
                end
                N_aug(:,:,:,t_idx) = RXDIF_3D_wMC_wAC_comb(initial, kp, d, alpha1, tx_params, t_out, h, dz, dt, bcs, M, E, nu, matX, matY, matZ);
            end
        else % three time point calibration, run twice, once eith each kp map
            if(i==1)
                initial = N(:,:,:,1);
                t_out = dt:dt:t(1);
                t_idx = t_out./dt + 1;
            else
                initial = N(:,:,:,2);
                t_out = [t(1)+dt:dt:t(2)]-t(1);
                t_idx = (t_out + t(1))./dt + 1;
            end
            N_aug(:,:,:,t_idx) = RXDIF_3D_wMC_wAC_comb(initial, kp, d, alpha1, tx_params, t_out, h, dz, dt, bcs, M, E, nu, matX, matY, matZ);
        end
    end
    if(ntp_cal == 2) % Add 3rd time point to end if needed
        N_aug = cat(4,N_aug, N(:,:,:,end)); 
    end
    
    kp = kp_out;
end