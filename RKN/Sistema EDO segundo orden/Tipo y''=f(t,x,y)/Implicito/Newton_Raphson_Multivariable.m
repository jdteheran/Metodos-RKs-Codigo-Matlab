function [K] = Newton_Raphson_Multivariable(f,x,t)
% f = [F1 ; F2 ; ... ; Fn]
% x = [X1 ; X2 ; ... ; Xn]

format long;
syms kp1x kp2x kp3x kp4x kp1y kp2y kp3y kp4y;

j = jacobian(f);
var = transpose(symvar(f));

for i=1:100
    xk = x - subs(j,var,x)\subs(f,var,x);
    
    xk = double(xk);    
    if abs(xk-x)<=t
        break;
    end
    
    x = xk;
end
K = xk(:,1);
end