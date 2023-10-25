function [S_mat] = buildStrainMat(N,E,v)
% Factor we need to multiply by
fac = E(:)/(1+v)/(1-2*v);
% Matrix we need to multiply by

[sy,sx] = size(N);

S_mat = sparse(sy*sx*3, sy*sx*3);
for j = 1:sx
    for i = 1:sy
        ind = i + (j - 1)*sy;
        yy_jump = sy*sx;
        xy_jump = sy*sx*2;
        
        %XX-stress
        S_mat(ind,ind) = fac(ind)*v;
        S_mat(ind,ind+yy_jump) = fac(ind)*(1-v);
        
        %YY-stress
        S_mat(ind+yy_jump,ind) = fac(ind)*(1-v);
        S_mat(ind+yy_jump,ind+yy_jump) = fac(ind)*v;
        
        %XY-stress
        S_mat(ind+xy_jump,ind+xy_jump) = fac(ind)*(1-2*v);
    end
end

end