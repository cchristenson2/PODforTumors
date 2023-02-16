function [Ux, Uy] = getDisplacementMaps_2D(N, M, E, nu, matX, matY)
    [sy,sx,nt] = size(N);
    
    Ux = zeros(sy,sx,nt);
    Uy = zeros(sy,sx,nt);
    for i = 1:nt
        [~, ~, temp_x, temp_y] = get_damper(matX, matY, N(:,:,i), M, E, nu);
        Ux(:,:,i) = reshape(temp_x,sy,sx);
        Uy(:,:,i) = reshape(temp_y,sy,sx);
    end
end