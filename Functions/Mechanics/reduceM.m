

function M_r = reduceM(M,V)
    [n,~] = size(M);
    n = n/2;
    
    XX = M(1:n,1:n);
    YY = M(n+1:2*n,n+1:2*n);
    XY = M(1:n,n+1:2*n);
    YX = M(n+1:2*n,1:n);
    
    M_r = [[(V' * XX * V),(V' * XY * V)];[(V' * YX * V),(V' * YY * V)]];

end