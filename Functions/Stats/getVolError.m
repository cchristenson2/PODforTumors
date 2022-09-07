%{ 
Calculate Volume % Error

Inputs:
    - Image a (typically measured data)
    - Image b (typically simulation)
    - threshold for significant voxels

Outputs:
    - Volume % Error of total region great than threshold

Contributors: Chase Christenson
%}
function err = getVolError(a,b,h,dz,thresh)
    idx1 = find(a>thresh);
    idx2 = find(b>thresh);
    
    if(~isempty(idx1))
        vol_a = numel(a(idx1))*h*h*dz;
    else
        vol_a = 0;
    end
    if(~isempty(idx2))
        vol_b = numel(b(idx2))*h*h*dz;
    else
        vol_b = 0;
    end
    
    err = 100*(vol_b - vol_a)/vol_a;
    if(vol_a == 0 && vol_b == 0)
        err = 0;
    end
end