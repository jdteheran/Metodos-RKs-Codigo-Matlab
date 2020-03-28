function [x,y] = Gauss_RK2_I(a,b,n,x0,y0)
%EJEMPLO: [x,y] = Gauss_RK2_I(0.01,1/2,100,0.0099999166667361106977490148082,0.000110517092175955069979947);

syms kx1 kx2 ky1 ky2;
format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);

KX1 = 1; KY1 = 1;
KX2 = 1; KY2 = 1;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    
    F = [f(t(i,1) + ((3-sqrt(3))/6)*h, x(i,1) + (1/4)*h*kx1 + ((3-2*sqrt(3))/12)*h*kx2, y(i,1) + (1/4)*h*ky1 + ((3-2*sqrt(3))/12)*h*ky2) - kx1 ;
         g(t(i,1) + ((3-sqrt(3))/6)*h, x(i,1) + (1/4)*h*kx1 + ((3-2*sqrt(3))/12)*h*kx2, y(i,1) + (1/4)*h*ky1 + ((3-2*sqrt(3))/12)*h*ky2) - ky1 ;
         f(t(i,1) + ((3+sqrt(3))/6)*h, x(i,1) + ((3+2*sqrt(3))/12)*h*kx1 + (1/4)*h*kx2, y(i,1) + ((3+2*sqrt(3))/12)*h*ky1 + (1/4)*h*ky2) - kx2 ;
         g(t(i,1) + ((3+sqrt(3))/6)*h, x(i,1) + ((3+2*sqrt(3))/12)*h*kx1 + (1/4)*h*kx2, y(i,1) + ((3+2*sqrt(3))/12)*h*ky1 + (1/4)*h*ky2) - ky2];
     
    K = Newton_Raphson_Multivariable(F, [KX1;KX2;KY1;KY2], 10^-8);
    
    KX1 = K(1,1);
    KX2 = K(2,1);
    KY1 = K(3,1);    
    KY2 = K(4,1);
    
    x(i+1,1) = x(i,1) + (h/2)*(KX1+KX2);
    y(i+1,1) = y(i,1) + (h/2)*(KY1+KY2);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end