function [y] = Gauss_RK1_I(a,b,n,y0)
%EJEMPLO: Gauss_RK1_I(0,7,200,pi^2);

syms k;
format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);

K1 = 0;

y(1,1) = y0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    K1 = Newton_Raphson(f(x(i,1) + (1/2)*h, y(i,1) + (1/2)*h*k) - k, K1, 10^-8);
    y(i+1,1) = y(i,1) + h*K1;
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end