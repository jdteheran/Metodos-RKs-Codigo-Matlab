function [K] = Newton_Raphson_Multivariable(f,x,t)
% f = [F1 ; F2 ; ... ; Fn]
% x = [X1 ; X2 ; ... ; Xn]

format long;
syms k1 k2 k3 k4;

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