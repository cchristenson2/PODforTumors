%{ 
Get projection matrix for a given set of snapshots

Inputs:
    - Snapshot matrix
    - desired rank

Outputs:
    - V; projection matrix 
    - U; full left handed vectors
    - S; singular values
    - V_; full right handed vectors

Contributors: Chase Christenson
%}

function [Vx, Vy, Vz] = getProjectionMatrix_MC(Ux, Uy, Uz,k)
    [~,~,temp1,temp2] = size(Ux);
    if(temp2==1) %2D
        Ux_vert = zeros(numel(Ux(:,:,1)), temp1);
        Uy_vert = zeros(numel(Uy(:,:,1)), temp1);
        for i = 1:temp1
            Ux_vert(:,i) = reshape(Ux(:,:,i),[],1);
            Uy_vert(:,i) = reshape(Uy(:,:,i),[],1);
        end
    else %3D
        Ux_vert = zeros(numel(Ux(:,:,:,1)), temp2);
        Uy_vert = zeros(numel(Uy(:,:,:,1)), temp2);
        Uz_vert = zeros(numel(Uz(:,:,:,1)), temp2);
        for i = 1:temp2
            Ux_vert(:,i) = reshape(Ux(:,:,:,i),[],1);
            Uy_vert(:,i) = reshape(Uy(:,:,:,i),[],1);
            Uz_vert(:,i) = reshape(Uz(:,:,:,i),[],1);
        end
    end
    
    if(temp2==1) %2D
        [V1,~,~] = svd(Ux_vert);
        [V2,~,~] = svd(Uy_vert);
        
        Vx = V1(:,1:k);
        Vy = V2(:,1:k);
        Vz = [];
    else
        [V1,~,~] = svd(Ux_vert);
        [V2,~,~] = svd(Uy_vert);
        [V3,~,~] = svd(Uz_vert);
        
        Vx = V1(:,1:k);
        Vy = V2(:,1:k);
        Vz = V3(:,1:k);
    end
end