function [stress_maps] = getStressMaps_2D(N, M, E, nu, matX, matY)
    [sy,sx,nt] = size(N);
    
    stress_maps = zeros(sy,sx,nt);
    
    for i = 1:nt
        stress_maps(:,:,i) = reshape(get_damper(matX, matY, N(:,:,i), M, E, nu), sy, sx);
    end
end