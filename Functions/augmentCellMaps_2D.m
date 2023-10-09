%{
Averages maps in order to fill in gaps between time points, no parameter estimation required

Inputs:
    N - cell maps to fill gaps between
    t - time of cell maps
    depth - number of averages to add

Outputs:
    N_aug - full time course with true maps, and averaged maps inbetween
    t_spacing - time of each map in the augmented cell data

%}


function [N_aug, t_spacing] = augmentCellMaps_2D(N, t, depth)

    N_aug = N;
    t_spacing = [0, t];
    
    thresh = 0.1;
    

    %% Create new midpoint for every depth, 1 = one average, 2 = two averages ...
    for j = 1:depth
        
        nt = numel(t_spacing); %Number of images in current dataset
        
        curr = 0;
        for i = 1:nt-1
            
            N_mid = (N_aug(:,:,i+curr)+N_aug(:,:,i+1+curr))./2;
            
            N_mid(N_mid<thresh) = 0;
            if(j==1)
                N_mid = imgaussfilt(N_mid, 0.5);
            end
            
            t_mid = (t_spacing(i+curr)+t_spacing(i+1+curr))/2;
            
            N_aug = cat(3, N_aug(:,:,1:i+curr), N_mid, N_aug(:,:,i+1+curr:end));
            t_spacing = [t_spacing(1:i+curr), t_mid, t_spacing(i+1+curr:end)];
            
            curr = curr+1;
        end
    end
    
end