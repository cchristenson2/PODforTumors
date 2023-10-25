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

function [V_d,Vx,Vy] = getDisplacementProjection_2D(Ux, Uy, k)

    [Vx,kx] = getProjectionMatrix(Ux,k);
    [Vy,ky] = getProjectionMatrix(Uy,k);
    
    n = size(Vx,1);
    
    V_d = sparse([Vx, zeros(n,ky); zeros(n,kx), Vy]);

end