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

function [V, k, V_, U, S] = getProjectionMatrix(N,k)
    [~,~,temp1,temp2] = size(N);
    if(temp2==1) %2D
        N_vert = zeros(numel(N(:,:,1)), temp1);
        for i = 1:temp1
            N_vert(:,i) = reshape(N(:,:,i),[],1);
        end
    else %3D
        N_vert = zeros(numel(N(:,:,:,1)), temp2);
        for i = 1:temp2
            N_vert(:,i) = reshape(N(:,:,:,i),[],1);
        end
    end

    %Compute SVD of snapshot matrix
    [U,S,V_] = svd(N_vert);
    
    %Compute cummulative energy
    sings = S(S>0);
    cummulative_e = cumsum(sings);
    max_e = max(cummulative_e);
    if(k == 0)
        k = numel(cummulative_e(cummulative_e < (max_e*0.995))) + 2;
    else
        k = k;
    end

    %Set projection basis
    V = U(:,1:k);
end