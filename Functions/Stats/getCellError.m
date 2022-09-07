%{ 
Calculate Cell % Error

Inputs:
    - Image a (typically measured data)
    - Image b (typically simulation)
    - threshold for significant voxels

Outputs:
    - Cell % Error of total region great than threshold

Contributors: Chase Christenson
%}
function err = getCellError(a,b,thresh)
    idx1 = find(a>thresh);
    idx2 = find(b>thresh);
    
    if(~isempty(idx1))
        cells_a = sum(a(idx1),'all');
    else
        cells_a = 0;
    end
    if(~isempty(idx2))
        cells_b = sum(b(idx2),'all');
    else
        cells_b = 0;
    end
    
    err = 100*(cells_b - cells_a)/cells_a;
    if(cells_a == 0 && cells_b == 0)
        err = 0;
    end
end