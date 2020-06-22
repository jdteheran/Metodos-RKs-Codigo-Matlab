function [x,y] = Gauss_RK1_I(a,b,n,x0,y0)
%EJEMPLO: [x,y] = Gauss_RK1_I(0.01,1/2,100,0.0099999166667361106977490148082,0.0001105170921759550699799478876);

syms kx1 ky1;
format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);

KX1 = 1; KY1 = 1;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    
    F = [f(t(i,1) + (1/2)*h, x(i,1) + (1/2)*h*kx1, y(i,1) + (1/2)*h*ky1) - kx1 ;
         g(t(i,1) + (1/2)*h, x(i,1) + (1/2)*h*kx1, y(i,1) + (1/2)*h*ky1) - ky1];
    K = Newton_Raphson_Multivariable(F, [KX1;KY1], 10^-8);
    
    %TENER EN CUENTA COMO EL ALGORITMO DE NEWTON_RAPHSON TE DEVUELVE LAS
    %VARIABLES
    KX1 = K(1,1);
    KY1 = K(2,1);    
    
    x(i+1,1) = x(i,1) + h*KX1;
    y(i+1,1) = y(i,1) + h*KY1;
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end