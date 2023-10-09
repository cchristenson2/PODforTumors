%{ 
Calculate CCC for whole image

Inputs:
    - Image a (typically measured data)
    - Image b (typically simulation)

Outputs:
    - CCC of images
%}
function ccc = CCC_calc(a,b)
    n = numel(a);

    a = a(:);
    b = b(:);
    
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