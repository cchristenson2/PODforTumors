function [matX, matY] = gradN_matrix(sx, sy, h, bcf)
% Initialize matrices for each
matX = zeros(sy*sx,sy*sx);
matY = zeros(sy*sx,sy*sx);


count = 0;
for x = 1:sx
    for y = 1:sy
        count = count + 1;
        itm(y,x) = count;
    end
end


% Fill our matrices
val = 1 / 2 / h;
for y = 1:sy
    for x = 1:sx
        count = itm(y,x);
        % Now we go through all of our boundary checks
        % On the boundaries we are going to set the gradient to be zero, so
        % we just don't change those rows

        
        % on no boundaries
        if and(bcf(y,x,1) == 0, bcf(y,x,2) == 0)
            matX(count, itm(y,x+1)) = val;
            matX(count, itm(y,x-1)) = -val;

            matY(count, itm(y+1,x)) = -val;
            matY(count, itm(y-1,x)) = val;

        % x is interior, y is not
        elseif and(bcf(y,x,1) ~= 0, bcf(y,x,2) == 0)
            matX(count, itm(y,x+1)) = val;
            matX(count, itm(y,x-1)) = -val;
        
        % y is interior, x is not
        elseif and(bcf(y,x,1) == 0, bcf(y,x,2) ~= 0)
            matY(count, itm(y+1,x)) = -val;
            matY(count, itm(y-1,x)) = val;
        % o/w we are out of the mask or on boundaries in both ways
        else
            continue
        end
    end
end

end