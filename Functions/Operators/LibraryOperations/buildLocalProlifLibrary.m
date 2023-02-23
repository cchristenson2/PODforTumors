function [LibB, LibH] = buildLocalProlifLibrary(N, V, k, num, bounds)
    for j = 1:k
        vec = linspace(bounds.low(j),bounds.up(j), num);
        eval(['LibB.Mode',num2str(j),'.vec = vec;']);
        eval(['LibH.Mode',num2str(j),'.vec = vec;']);
        for i = 1:num
            tempB = zeros(k,k);
            tempH = zeros(k,k^2);
            
            temp_map = V(:,j) * vec(i);
            for l = 1:numel(N)
                tempB = tempB + V(l,:)' * temp_map(l) * V(l,:);
                tempH = tempH + V(l,:)' * temp_map(l) * kron(V(l,:),V(l,:));
            
            end
            
            eval(['LibB.Mode',num2str(j),'.OP',num2str(i),' = tempB;']);
            eval(['LibH.Mode',num2str(j),'.OP',num2str(i),' = tempH;']);
        end
    end
end