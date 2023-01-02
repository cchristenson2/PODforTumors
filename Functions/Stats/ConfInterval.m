function delta = ConfInterval(J,res)
    Js = diag(J'*J); % J = Jacobian
    [~,R] = qr(J,0);
    rankJ = size(J,2);
    Rinv = R \ eye(size(R));
    pinvJTJ = Rinv*Rinv';
    
    n=numel(res);
    mse = sum(abs(res).^2)/(n-rankJ);
    Sigma = mse*pinvJTJ;
    % Sigma = covariance matrix

%     S2 = sqrtm(Sigma);  

    se = sqrt(diag(Sigma));
    alpha = 0.05; 
    delta = se * tinv(1-alpha/2,n);
end