%{
Build operator for diffusivity with boundary conditions around the breast

Input:
    - number of elements in domain
    - diffusivity
    - grid spacing
    - boundary map

Output:
    - A2 operator - for gradient on diffusion

Contributors: Chase Christenson, Graham Pash
%}
function [A2] = assembleGrad(N, h, dz, bcs)

    div = assembleDiv(N, h, dz, bcs);

    
    A2 = div;
    
    A2 = applyBC(A2, bcs, h, dz);
    
%     A = L+div;
end

%Build divergence for div(d) * div(N) with boundaries included
%Boundaries have to be built now since d is included in this operator
function [div] = assembleDiv(N, h, dz, bcs)
    [sy,sx,sz] = size(N);
    
    if(sz==1) %2D build
        n = sy*sx;
        div = zeros(n,n);
        
        for x = 1:sx
            for y = 1:sy
                i = y + (x-1)*sy;
                % Y Direction
                if(bcs(y,x,1)==0) %Only on interior nodes
                    div(i,i+1) = 1/(2*h);
                    div(i,i-1) = -1/(2*h);
                end
                % X Direction
                if(bcs(y,x,2)==0) %Only on interior nodes
                    div(i,i+sy) = 1/(2*h);
                    div(i,i-sy) = -1/(2*h);
                end
            end
        end
        
        div = sparse(div);
        
        
    else %3D build
        n = sy*sx*sz;
        div = zeros(n,n);
        
        for z = 1:sz
            for x = 1:sx
                for y = 1:sy
                    i = y + (x-1)*sy + (z-1)*sx*sy;
                    
                    % Y Direction
                    if(bcs(y,x,z,1)==0) %Only on interior nodes
                        div(i,i+1) = 1/(2*h);
                        div(i,i-1) = -1/(2*h);
                    end
                    
                    % X Direction
                    if(bcs(y,x,z,2)==0) %Only on interior nodes
                        div(i,i+sy) = 1/(2*h);
                        div(i,i-sy) = -1/(2*h);
                    end

                    % Z Direction
                    if(bcs(y,x,z,3)==0) %Only on interior nodes
                        div(i,i+(sy*sx)) = 1/(2*dz);
                        div(i,i-(sy*sx)) = -1/(2*dz);
                    end
                    
                end
            end
        end
        
        div = sparse(div);
    end
end


function [L] = applyBC(L, bcs, h, dz, d)
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
                        L(i,i+1) = 2*1/h;
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
                        L(i,i+sy) = 2*1/h;
                        if(i-sy>=1)
                            L(i,i-sy) = 0;
                        end
                    elseif(bcs(y,x,2)==1) %Check if right wall
                        L(i,i-sy) = 2*1/h;
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
                            L(i,i+1) = 2*1/h;
                            if(i-1>=1)
                                L(i,i-1) = 0;
                            end
                        elseif(bcs(y,x,z,1)==1) %Check bottom wall
                            L(i,i-1) = 2*1/h;
                            if(i+1<=n_3D)
                                L(i,i+1) = 0;
                            end
                        end

                        if(bcs(y,x,z,2)==-1) %Check if left wall
                            L(i,i+sy) = 2*1/h;
                            if(i-sy>=1)
                                L(i,i-sy) = 0;
                            end
                        elseif(bcs(y,x,z,2)==1) %Check if right wall
                            L(i,i-sy) = 2*1/h;
                            if(i+sy<=n_3D)
                                L(i,i+sy) = 0;
                            end
                        end
                        
                        if(bcs(y,x,z,3)==-1) %Check if under wall
                            L(i,i+sy*sx) = 2*1/dz;
%                             disp(size(L));
                            if(i-(sy*sx)>=1)
                                L(i,i-(sy*sx)) = 0;
                            end
                        elseif(bcs(y,x,z,3)==1) %Check if above wall
                            L(i,i-sy*sx) = 2*1/dz;
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
