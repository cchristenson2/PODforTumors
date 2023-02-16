function [stress_maps] = getStressMaps_3D(N, M, E, nu, matX, matY, matZ)
    [sy,sx,sz,nt] = size(N);
    
    stress_maps = zeros(sy,sx,sz,nt);
    
    for i = 1:nt
        stress_maps(:,:,:,i) = reshape(get_damper_3D(matX, matY, matZ, N(:,:,:,i), M, E, nu), sy, sx);
    end
end