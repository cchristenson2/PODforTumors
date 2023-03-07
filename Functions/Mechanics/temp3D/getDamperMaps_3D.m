function Dampers = getDamperMaps_3D(N, M, E, nu, matX, matY, matZ)
    [sy,sx,sz,nt] = size(N);
    
    Dampers = zeros(sy,sx,sz,nt);
    for i = 1:nt
        temp = get_damper_3D(matX, matY, matZ, N(:,:,:,i), M, E, nu);
        Dampers(:,:,:,i) = reshape(temp,sy,sx,sz);
    end
end