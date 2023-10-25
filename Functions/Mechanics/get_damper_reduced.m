%Improved cell gradient calculation

function [damper, s_vm] = get_damper_reduced(Vd_gradXY_V, Ve_gradXYY_Vd, Vs_SMAT_Ve, N, M_r, Vs)
    k = numel(N);

    lambda1 = 2.5e-3;
    lambda2 = 2.5e-3; 
%     lambda2 = 1000; 

    % Calculate displacement
    grad_N = Vd_gradXY_V * [N(:); N(:)];
    Uv_r = M_r\(lambda1.*grad_N);
    %Return to full size
%     Uv = Vd * Uv_r;
    
    
    %Strain calculation
    e_r = Ve_gradXYY_Vd * [Uv_r; Uv_r(1:k)];
    %Return to full size
%     e = Ve * e_r;
    
    
    %Stress calculation
    stress_r = Vs_SMAT_Ve * e_r;
    %Return to full size
    stress = Vs * stress_r;
    
    num = numel(stress)/3;
    
    s_vm = gsqrt(stress(1:num).^2-stress(1:num).*stress(num+1:num*2)+stress(num+1:num*2).^2+3*stress(num*2+1:num*3).^2);
    
    damper = exp(-1*s_vm * lambda2);
end