%{ 
Calculate CCC for overlapping regions above a threshold

Inputs:
    - Image a (typically measured data)
    - Image b (typically simulation)
    - threshold for significant voxels

Outputs:
    - CCC of image overlap

Contributors: Chase Christenson
%}
function ccc = CCC_overlap(a,b,thresh)
    
    
    idx1 = find(a>thresh);
    idx2 = find(b>thresh);
    idx3 = intersect(idx1, idx2);
    
    a = a(idx3);
    b = b(idx3);
    
    n = numel(a);
    
    
    mu_a = mean(a);
    mu_b = mean(b);

    cov = 0;
    varx = 0;
    vary = 0;
    for i = 1:n
       cov = cov + ((a(i)-mu_a)*(b(i)-mu_b)); 
       varx = varx + (a(i)-mu_a)^2;
       vary = vary + (b(i)-mu_b)^2;
    end
    cov = cov/n;
    varx = varx/n;
    vary = vary/n;

    ccc = 2*cov/(varx+vary+(mu_a - mu_b)^2);
end