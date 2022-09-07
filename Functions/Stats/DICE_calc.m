%{ 
Calculate DICE

Inputs:
    - Image a (typically measured data)
    - Image b (typically simulation)
    - threshold for significant voxels

Outputs:
    - DICE of image overlap

Contributors: Chase Christenson
%}
function dice = DICE_calc(a,b,thresh)
    idx1 = find(a>thresh);
    idx2 = find(b>thresh);
    idx3 = intersect(idx1, idx2);
    dice = 2*numel(idx3)/(numel(idx1)+numel(idx2));
end