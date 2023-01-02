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
function [A] = assembleA(N, d, h, dz, bcs)
    L = assembleL(N, h, dz);
    L = applyBC(L, bcs, h, dz);
    L = d*L;
    A = L;
end

function [L] = assembleL(N, h, dz)
    [sy,sx,sz] = size(N);
    
    if(sz==1) %2D build
        n = sy*sx;
        e = ones(n,1);
        Y = spdiags([e -2*e e]./h, [-1 0 1],n,n);
        X = spdiags([e -2*e e]./h, [-sy 0 sy],n,n);
        L = sparse(X+Y);
        
    else %3D build
        n = sy*sx*sz;
        e = ones(n,1);
        Y = spdiags([e -2*e e]./h, [-1 0 1],n,n);
        X = spdiags([e -2*e e]./h, [-sy 0 sy],n,n);
        Z = spdiags([e -2*e e]./dz, [-1*(sy*sx) 0 (sy*sx)],n,n);
        L = sparse(X+Y+Z);
    end

end

function [L] = applyBC(L, bcs, h, dz)
% Apply Neumann BC discrete laplacian
    [sy,sx,sz,~] = size(bcs);
    n_3D = sy*sx*sz;
    n_2D = sy*sx;
    
    if(sz<=2) %2D operator
        for x = 1:sx
            for y = 1:sy
%                 disp(['(y,x) = (',num2str(y),',',num2str(x),')']);
                i = y + (x-1)*sy;
                %Remove outside of mask contributions
                if(bcs(y,x,1)==2)
                    L(i,:) = 0;
                else

                    if(bcs(y,x,1)==-1) %Check if top wall
                        L(i,i+1) = 2/h;
                        if(i-1>=1)
                            L(i,i-1) = 0;
                        end
                    elseif(bcs(y,x,1)==1) %Check bottom wall
                        L(i,i-1) = 2/h;
                        if(i+1<=n_2D)
                            L(i,i+1) = 0;
                        end
                    end

                    if(bcs(y,x,2)==-1) %Check if left wall
                        L(i,i+sy) = 2/h;
                        if(i-sy>=1)
                            L(i,i-sy) = 0;
                        end
                    elseif(bcs(y,x,2)==1) %Check if right wall
                        L(i,i-sy) = 2/h;
                        if(i+sy<=n_2D)
                            L(i,i+sy) = 0;
                        end
                    end
                end
            end
        end
    else %3D operator
        for z = 1:sz
            for x = 1:sx
                for y = 1:sy
        %             disp(['(y,x) = (',num2str(y),',',num2str(x),')']);
                    i = y + (x-1)*sy + (z-1)*sx*sy;
%                     disp([x, y, z]);
                    %Remove outside of mask contributions
                    if(bcs(y,x,z,1)==2)
                        L(i,:) = 0;
                    else
                        if(bcs(y,x,z,1)==-1) %Check if top wall
                            L(i,i+1) = 2/h;
                            if(i-1>=1)
                                L(i,i-1) = 0;
                            end
                        elseif(bcs(y,x,z,1)==1) %Check bottom wall
                            L(i,i-1) = 2/h;
                            if(i+1<=n_3D)
                                L(i,i+1) = 0;
                            end
                        end

                        if(bcs(y,x,z,2)==-1) %Check if left wall
                            L(i,i+sy) = 2/h;
                            if(i-sy>=1)
                                L(i,i-sy) = 0;
                            end
                        elseif(bcs(y,x,z,2)==1) %Check if right wall
                            L(i,i-sy) = 2/h;
                            if(i+sy<=n_3D)
                                L(i,i+sy) = 0;
                            end
                        end
                        
                        if(bcs(y,x,z,3)==-1) %Check if under wall
                            L(i,i+sy*sx) = 2/dz;
%                             disp(size(L));
                            if(i-(sy*sx)>=1)
                                L(i,i-(sy*sx)) = 0;
                            end
                        elseif(bcs(y,x,z,3)==1) %Check if above wall
                            L(i,i-sy*sx) = 2/dz;
                            if(i+(sy*sx)<=n_3D)
                                L(i,i+(sy*sx)) = 0;
                            end
                        end
                        
%                         disp(size(L));
                    end
                end
            end
        end
    end
end
