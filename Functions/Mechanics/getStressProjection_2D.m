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

function [V_e] = getStressProjection_2D(Sxx, Syy, Sxy, k)

    [Vxx,kxx] = getProjectionMatrix(Sxx,k);
    [Vyy,kyy] = getProjectionMatrix(Syy,k);
    [Vxy,kxy] = getProjectionMatrix(Sxy,k);
    
    n = size(Vxx,1);
    
    V_e = sparse([Vxx, zeros(n,kyy), zeros(n,kxy); zeros(n,kxx), Vyy, zeros(n,kxy); zeros(n,kxx), zeros(n,kyy), Vxy]);

end