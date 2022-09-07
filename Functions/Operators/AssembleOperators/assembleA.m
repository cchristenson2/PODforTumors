%{
Build operator for diffusivity with boundary conditions around the breast

Input:
    - number of elements in domain
    - diffusivity
    - grid spacing
    - boundary map

Output:
    - A operator

Contributors: Chase Christenson, Graham Pash
%}
function [A] = assembleA(n, d, h, bcs)
    L = assembleL(n);
    L = applyBC(L, bcs);
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

function [L] = applyBC(L, bcs)
% Apply Neumann BC discrete laplacian

    n = sqrt(size(L,1));
    [sy,sx,~] = size(bcs);

    for x = 1:sx
        for y = 1:sy
%             disp(['(y,x) = (',num2str(y),',',num2str(x),')']);
            i = y + (x-1)*sx;
            %Remove outside of mask contributions
            if(bcs(y,x,1)==2)
                L(i,:) = 0;
            else
                
                if(bcs(y,x,1)==-1) %Check if top wall
                    L(i,i+1) = 2;
                elseif(bcs(y,x,1)==1) %Check bottom wall
                    L(i,i-1) = 2;
                end
                
                if(bcs(y,x,2)==-1) %Check if left wall
                    L(i,i+n) = 2;
                elseif(bcs(y,x,2)==1) %Check if right wall
                    L(i,i-n) = 2;
                end
            end
        end
    end
end
