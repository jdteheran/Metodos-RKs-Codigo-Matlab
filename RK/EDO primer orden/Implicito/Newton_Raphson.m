function [K] = Newton_Raphson(f,x,t)

format long;
syms k;

df = diff(f);

for i = 1:1000
    xk = x - subs(f,x)/subs(df,x);    
    
    xk = double(xk);
    
    if abs(xk-x)<t      
        break;
    end
    x = xk;
end

K = xk;
end

