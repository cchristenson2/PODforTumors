function [Dampers] = getDamperMaps_2D(N, M, E, nu, matX, matY)
    [sy,sx,nt] = size(N);
    
    Dampers = zeros(sy,sx,nt);
    for i = 1:nt
        temp = get_damper(matX, matY, N(:,:,i), M, E, nu);
        Dampers(:,:,i) = reshape(temp,sy,sx);
    end
end