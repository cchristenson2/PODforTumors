function [Lib] = buildLocalDiffuseLibrary(N, V, Vd, k, num, bounds, h, dz, bcs)
    for j = 1:k
        vec = linspace(bounds.low(j),bounds.up(j), num);
        eval(['Lib.Mode',num2str(j),'.vec = vec;']);
        for i = 1:num
            eval(['Lib.Mode',num2str(j),'.OP',num2str(i),' = V'' * assembleA(N, (Vd(:,j)*vec(i)), h, dz, bcs) * V;']);
        end
    end
end