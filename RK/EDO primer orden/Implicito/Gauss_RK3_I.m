function [y] = Gauss_RK3_I(a,b,n,y0)
%EJEMPLO: Gauss_RK3_I(0,7,100,pi^2); Se toma sus 10 min

syms k1 k2 k3;
format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);

K1 = 0;
K2 = 0;
K3 = 0;

y(1,1) = y0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    
    F = [f(x(i,1) + ((5-sqrt(15))/10)*h, y(i,1) + (5/36)*h*k1 + ((10-3*sqrt(15))/45)*h*k2 + ((25-6*sqrt(15))/180)*h*k3) - k1 ;
         f(x(i,1) + (1/2)*h, y(i,1) + ((10+3*sqrt(15))/72)*h*k1 + (2/9)*h*k2 + ((10-3*sqrt(15))/72)*h*k3) - k2 ;
         f(x(i,1) + ((5+sqrt(15))/10)*h, y(i,1) + ((25+6*sqrt(15))/180)*h*k1 + ((10+3*sqrt(15))/45)*h*k2 + (5/36)*h*k3) - k3];
    K = Newton_Raphson_Multivariable(F, [K1;K2;K3], 10^-8);
    
    %TENER EN CUENTA COMO EL ALGORITMO DE NEWTON_RAPHSON TE DEVUELVE LAS
    %VARIABLES
    K1 = K(1,1);
    K2 = K(2,1);
    K3 = K(3,1);
    
    y(i+1,1) = y(i,1) + (h/18)*(5*K1+8*K2+5*K3);
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end