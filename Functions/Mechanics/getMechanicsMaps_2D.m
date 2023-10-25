function [Damper, VM, Ux, Uy, Exx, Eyy, Exy, Sxx, Syy, Sxy] = getMechanicsMaps_2D(N, M, E, nu, matX, matY)
    [sy,sx,nt] = size(N);
    
    Damper = zeros(sy,sx,nt);
    VM = zeros(sy,sx,nt);
    Ux = zeros(sy,sx,nt);
    Uy = zeros(sy,sx,nt);
    Exx = zeros(sy,sx,nt);
    Eyy = zeros(sy,sx,nt);
    Exy = zeros(sy,sx,nt);
    Sxx = zeros(sy,sx,nt);
    Syy = zeros(sy,sx,nt);
    Sxy = zeros(sy,sx,nt);
    
    
    for i = 1:nt
        [temp_damper, temp_vm, temp_x, temp_y, temp_Exx, temp_Eyy, temp_Exy, temp_Sxx, temp_Syy, temp_Sxy] = ...
            get_damper(matX, matY, N(:,:,i), M, E, nu);
        
        Damper(:,:,i) = reshape(temp_damper, sy,sx);
        VM(:,:,i) = reshape(temp_vm, sy,sx);
        
        Ux(:,:,i) = reshape(temp_x,sy,sx);
        Uy(:,:,i) = reshape(temp_y,sy,sx);
        
        Exx(:,:,i) = reshape(temp_Exx,sy,sx);
        Eyy(:,:,i) = reshape(temp_Eyy,sy,sx);
        Exy(:,:,i) = reshape(temp_Exy,sy,sx);
        
        Sxx(:,:,i) = reshape(temp_Sxx,sy,sx);
        Syy(:,:,i) = reshape(temp_Syy,sy,sx);
        Sxy(:,:,i) = reshape(temp_Sxy,sy,sx);
        
    end
end