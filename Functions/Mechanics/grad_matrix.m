%{
Function to build gradient operators for current domain with boundaries
2D inputs only

inputs:
    - h: in-plane resolution
    - bcf: boundary condition function
outputs:
    - matX, matY: gradient operators for x and y directions
%}

function [matX, matY] = grad_matrix(h, bcf)
[sy,sx,~] = size(bcf);

% Initialize matrices for each
matX = zeros(sy*sx,sy*sx);
matY = zeros(sy*sx,sy*sx);

% Fill our matrices
val = 1/(2*h);
count = 1;
for x = 1:sx
    for y = 1:sy
        % Now we go through all of our boundary checks
        % On the boundaries we are going to set the gradient to be zero, so
        % we just don't change those rows
        
        if(bcf(y,x,2)==0)
            matX(count, count+sy) = val;
            matX(count, count-sy) = -val;
        elseif(bcf(y,x,2)==-1)
            matX(count, count+sy) = val*2;
            matX(count, count) = -val*2;
        elseif(bcf(y,x,2)==1)
            matX(count, count) = val*2;
            matX(count, count-sy) = -val*2;
        end
        
        if(bcf(y,x,1)==0)
            matY(count, count+1) = val;
            matY(count, count-1) = -val;
        elseif(bcf(y,x,1)==-1)
            matY(count, count+1) = val*2;
            matY(count, count) = -val*2;
        elseif(bcf(y,x,1)==1)
            matY(count, count) = val*2;
            matY(count, count-1) = -val*2;
        end
        
        count = count + 1;
    end
end



end