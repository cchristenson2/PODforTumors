%{ 
Builds ROI to use for assigning parameters

Inputs:
    - tumor struct
    - number of calibrated time points

Outputs:
    - ROI for proliferating region

Contributors: Casey Stowers
%}
function [ROI] = buildROI(tumor, ntp_cal)
    BreastMask = tumor.Mask;
    N = tumor.N;
    ROI = zeros(size(BreastMask));
    BigROI = zeros(size(BreastMask));
    
    [~,~,temp] = size(BreastMask);
    if(temp==1) %Tumor is 2D
        BigROI(:,:) = sum(N(:,:,1:1+ntp_cal), 3);
        BigROI(BigROI>0) = 1;
        BigROI(isnan(BigROI)) = 0;
        BigROIswell = imdilate(BigROI(:,:), strel('disk',ceil(3*4/tumor.CF)));
        BigROIswell = bwconvhull(BigROIswell);
%         ROI(:,:) = BigROIswell.*BW(:,:);
        ROI(:,:) = BigROIswell;
        
    else %Tumor is 3D
        sz = temp;
        for gg = 1:sz
            temp = sum(N(:,:,gg,1:1+ntp_cal), 4);
            temp(temp>0) = 1;
            temp(isnan(temp)) = 0;
            BigROIswell = imdilate(temp, strel('disk',ceil(3*4/tumor.CF)));
            BigROIswell = bwconvhull(BigROIswell);
%             ROI(:,:,gg) = BigROIswell.*BW(:,:,gg);
            ROI(:,:,gg) = BigROIswell;
        end
    end
end