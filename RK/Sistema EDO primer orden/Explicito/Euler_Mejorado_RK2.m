function [x,y] = Euler_Mejorado_RK2(a,b,n,x0,y0)
%EJEMPLO: [x,y]=Euler_Mejorado_RK2(0.01,1/2,100,0.0099999166667361106977490148082,0.000110517092175955069979947);

format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);
xtemp = 0;
ytemp = 0;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    xtemp = x(i,1) + h*f(t(i,1),x(i,1),y(i,1));
    ytemp = y(i,1) + h*g(t(i,1),x(i,1),y(i,1));
    x(i+1,1) = x(i,1) + h/2*(f(t(i,1),x(i,1),y(i,1)) + f(a + i*h,xtemp,ytemp));
    y(i+1,1) = y(i,1) + h/2*(g(t(i,1),x(i,1),y(i,1)) + g(a + i*h,xtemp,ytemp));
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end