function [matX, matY, matZ] = gradN_matrix(sx, sy, sz, h, dz, bcf)
% Initialize matrices for each
matX = zeros(sy*sx*sz,sy*sx*sz);
matY = zeros(sy*sx*sz,sy*sx*sz);
matZ = zeros(sy*sx*sz,sy*sx*sz);

count = 0;
for z = 1:sz
    for x = 1:sx
        for y = 1:sy
            count = count + 1;
            itm(y,x,z) = count;
        end
    end
end


% Fill our matrices
val_h = 1 / 2 / h;
val_z = 1 / 2 / dz;
for z = 1:sz
    for y = 1:sy
        for x = 1:sx
            count = itm(y,x,z);
            % Now we go through all of our boundary checks
            % On the boundaries we are going to set the gradient to be zero, so
            % we just don't change those rows
            
            if(bcf(y,x,z,1)==0)
                matY(count, itm(y+1,x,z)) = val_h;
                matY(count, itm(y-1,x,z)) = -val_h;
            end
            
            if(bcf(y,x,z,2)==0)
                matX(count, itm(y,x+1,z)) = val_h;
                matX(count, itm(y,x-1,z)) = -val_h;
            end
            
            if(bcf(y,x,z,3)==0)
                matZ(count, itm(y,x,z+1)) = val_z;
                matZ(count, itm(y,x,z-1)) = -val_z;
            end
        end
    end
end

end