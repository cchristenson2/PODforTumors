function [Ux, Uy, Uz] = getDisplacementMaps_3D(N, M, E, nu, matX, matY, matZ)
    [sy,sx,sz,nt] = size(N);
    
    Ux = zeros(sy,sx,sz,nt);
    Uy = zeros(sy,sx,sz,nt);
    for i = 1:nt
        [~, ~, temp_x, temp_y, temp_z] = get_damper_3D(matX, matY, matZ, N(:,:,:,i), M, E, nu);
        Ux(:,:,:,i) = reshape(temp_x,sy,sx,sz);
        Uy(:,:,:,i) = reshape(temp_y,sy,sx,sz);
        Uz(:,:,:,i) = reshape(temp_z,sy,sx,sz);
    end
end