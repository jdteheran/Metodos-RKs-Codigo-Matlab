function [y] = Euler_Mejorado_RK2(a,b,n,y0)
%EJEMPLO: Euler_Mejorado_RK2(0,7,200,pi^2);

format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);
ytemp = 0;

y(1,1) = y0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    ytemp = y(i,1) + h*f(x(i,1),y(i,1));
    y(i+1,1) = y(i,1) + h/2*(f(x(i,1),y(i,1)) + f(a + i*h,ytemp));
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end

