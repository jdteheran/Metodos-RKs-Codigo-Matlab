function [y] = Gauss_RK2_I(a,b,n,y0)
%EJEMPLO: Gauss_RK2_I(0,7,200,pi^2);

syms k1 k2;
format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);

K1 = 0;
K2 = 0;

y(1,1) = y0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    
    F = [f(x(i,1) + ((3-sqrt(3))/6)*h, y(i,1) + (1/4)*h*k1 + ((3-2*sqrt(3))/12)*h*k2) - k1 ;
         f(x(i,1) + ((3+sqrt(3))/6)*h, y(i,1) + ((3+2*sqrt(3))/12)*h*k1 + (1/4)*h*k2) - k2];
    K = Newton_Raphson_Multivariable(F, [K1;K2], 10^-8);
    
    %TENER EN CUENTA COMO EL ALGORITMO DE NEWTON_RAPHSON TE DEVUELVE LAS
    %VARIABLES
    K1 = K(1,1);
    K2 = K(2,1);
    
    y(i+1,1) = y(i,1) + (h/2)*(K1+K2);
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end