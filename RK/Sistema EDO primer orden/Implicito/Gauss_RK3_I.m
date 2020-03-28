function [x,y] = Gauss_RK3_I(a,b,n,x0,y0)
%EJEMPLO: [x,y] = Gauss_RK3_I(0.01,1/2,100,0.0099999166667361106977490148082,0.000110517092175955069979947);

syms kx1 kx2 kx3 ky1 ky2 ky3;
format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);

KX1 = 1; KY1 = 1;
KX2 = 1; KY2 = 1;
KX3 = 1; KY3 = 1;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    
    F = [f(t(i,1) + ((5-sqrt(15))/10)*h, x(i,1) + (5/36)*h*kx1 + ((10-3*sqrt(15))/45)*h*kx2 + ((25-6*sqrt(15))/180)*h*kx3, y(i,1) + (5/36)*h*ky1 + ((10-3*sqrt(15))/45)*h*ky2 + ((25-6*sqrt(15))/180)*h*ky3) - kx1 ;
         g(t(i,1) + ((5-sqrt(15))/10)*h, x(i,1) + (5/36)*h*kx1 + ((10-3*sqrt(15))/45)*h*kx2 + ((25-6*sqrt(15))/180)*h*kx3, y(i,1) + (5/36)*h*ky1 + ((10-3*sqrt(15))/45)*h*ky2 + ((25-6*sqrt(15))/180)*h*ky3) - ky1 ;
         f(t(i,1) + (1/2)*h, x(i,1) + ((10+3*sqrt(15))/72)*h*kx1 + (2/9)*h*kx2 + ((10-3*sqrt(15))/72)*h*kx3, y(i,1) + ((10+3*sqrt(15))/72)*h*ky1 + (2/9)*h*ky2 + ((10-3*sqrt(15))/72)*h*ky3) - kx2 ;
         g(t(i,1) + (1/2)*h, x(i,1) + ((10+3*sqrt(15))/72)*h*kx1 + (2/9)*h*kx2 + ((10-3*sqrt(15))/72)*h*kx3, y(i,1) + ((10+3*sqrt(15))/72)*h*ky1 + (2/9)*h*ky2 + ((10-3*sqrt(15))/72)*h*ky3) - ky2 ;
         f(t(i,1) + ((5+sqrt(15))/10)*h, x(i,1) + ((25+6*sqrt(15))/180)*h*kx1 + ((10+3*sqrt(15))/45)*h*kx2 + (5/36)*h*kx3, y(i,1) + ((25+6*sqrt(15))/180)*h*ky1 + ((10+3*sqrt(15))/45)*h*ky2 + (5/36)*h*ky3) - kx3 ;
         g(t(i,1) + ((5+sqrt(15))/10)*h, x(i,1) + ((25+6*sqrt(15))/180)*h*kx1 + ((10+3*sqrt(15))/45)*h*kx2 + (5/36)*h*kx3, y(i,1) + ((25+6*sqrt(15))/180)*h*ky1 + ((10+3*sqrt(15))/45)*h*ky2 + (5/36)*h*ky3) - ky3];
     
    K = Newton_Raphson_Multivariable(F, [KX1;KX2;KX3;KY1;KY2;KY3], 10^-8);
    
    KX1 = K(1,1);
    KX2 = K(2,1);
    KX3 = K(3,1);
    KY1 = K(4,1);    
    KY2 = K(5,1);
    KY3 = K(6,1);
    
    x(i+1,1) = x(i,1) + (h/18)*(5*KX1+8*KX2+5*KX3);
    y(i+1,1) = y(i,1) + (h/18)*(5*KY1+8*KY2+5*KY3);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end