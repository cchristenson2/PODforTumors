function [N_aug,kp] = Augment_wMC_comb_2D(N, t, h, dt, bcs, bounds, ntp_cal, tx_params, M, E, nu, matX, matY)
    thresh = 0.15;

    [sy,sx,~] = size(N);
    N_aug = zeros(sy,sx,t(end)/dt+1);
    N_aug(:,:,1) = N(:,:,1);
    
    %% Loop once for 2 time point calibrations, twice for 3
    for i = 1:ntp_cal
    %% Get parameters for forward run
        N1 = N(:,:,i);
        N2 = N(:,:,i+1);
        
        if(i==1)
            t1 = 0;
        else
            t1 = t(i-1);
        end
        t2 = t(i);
    
        %Proliferation - estimation of change in cells based on logisitc growth equation
        kp = -1*log(N2./N1)./(t2-t1);
        kp(kp<bounds.kp_bounds(1)) = bounds.kp_bounds(1);
        kp(kp>bounds.kp_bounds(end)) = bounds.kp_bounds(end);
        kp(isnan(kp)) = 0;

        %Diffusivity - estimation of change in radius over time (area of circle assumption)
        A1 = numel(N1(N1>thresh))*h^2;
        A2 = numel(N2(N2>thresh))*h^2;
        
        %Alpha parameters
        alpha1 = (bounds.alpha_bounds(1) + bounds.alpha_bounds(end))/2;
%         alpha2 = (bounds.alpha_bounds(1) + bounds.alpha_bounds(end))/2;
        
        del_r = sqrt(A2/pi) - sqrt(A1/pi);
        d = del_r^2 / (4*(t2-t1));

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
                    initial = N(:,:,1);
                    t_out = dt:dt:t(1);
                    t_idx = t_out./dt + 1;
                else
                    initial = N(:,:,2);
                    t_out = [t(1)+dt:dt:t(2)]-t(1);
                    t_idx = (t_out + t(1))./dt + 1;
                end
                N_aug(:,:,t_idx) = RXDIF_2D_wMC_comb(initial, kp, d, alpha1, tx_params, t_out, h, dt, bcs, M, E, nu, matX, matY);
            end
        else % three time point calibration, run twice, once eith each kp map
            if(i==1)
                initial = N(:,:,1);
                t_out = dt:dt:t(1);
                t_idx = t_out./dt + 1;
            else
                initial = N(:,:,2);
                t_out = [t(1)+dt:dt:t(2)]-t(1);
                t_idx = (t_out + t(1))./dt + 1;
            end
            N_aug(:,:,t_idx) = RXDIF_2D_wMC_comb(initial, kp, d, alpha1, tx_params, t_out, h, dt, bcs, M, E, nu, matX, matY);
        end
    end
    if(ntp_cal == 2) % Add 3rd time point to end if needed
        N_aug = cat(3,N_aug, N(:,:,end)); 
    end
    
    kp = kp_out;
end