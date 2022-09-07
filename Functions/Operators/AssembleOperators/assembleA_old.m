%{
Build operator for diffusivity with square domain assumption

Input:
    - number of elements in domain
    - diffusivity
    - grid spacing

Output:
    - A operator

Contributors: Graham Pash
%}
function [A] = assembleA_old(n, d, h)
    L = assembleL(n);
    L = applyBC(L);
%     L = applyBC_v2(L, bcs);
    L = (d/h^2)*L;
    A = L;
end

function [L] = assembleL(n)
    n = sqrt(n); %temp added
    n = int16(n);
    I = eye(n, 'like', sparse(n,n));
    e = ones(n,1);
    T = spdiags([e -4*e e],[-1 0 1],n,n);
    S = spdiags([e e],[-1 1],n,n);
    L = (kron(I,T) + kron(S,I));
end

function [L] = applyBC_v2(L, bcs)
% Apply Neumann BC discrete laplacian

    n = sqrt(size(L,1));
    [sy,sx,~] = size(bcs)
    
    for x = 1:sx
        for y = 1:sy
            i = y+(x-1)*x;
            %Remove outside of mask contributions
            if(bcs(y,x,1)==2)
                L(i,:) = 0;
            else
                if(bcs(y,x,1)==-1) %Check if top wall
                    L(i,i-1) = 2;
                elseif(bcs(y,x,1)==1) %Check bottom wall
                    L(i,i+1) = 2;
                end
                if(bcs(y,x,2)==-1) %Check if left wall
                    L(i,i-n) = 2;
                elseif(bcs(y,x,2)==1) %Check if right wall
                    L(i,i+n) = 2;
                end
            end
        end
    end
end

function [L] = applyBC(L)
% Apply Neumann BC discrete laplacian

    n = sqrt(size(L,1));

    % Apply BC to left wall
    for i = 1:n:n^2
        L(i, i+1) = 2;
    end

    % Apply BC to right wall
    for i = 1:n:n^2
        L(i, i+1) = 2;
    end

    % Apply BC to bottom wall
    for i = 1:n
        L(i, i+n) = 2;
    end

    % Apply BC to top wall
    for i = 1:n
        L(end-i, end-i-n) = 2;
    end
end