function [y] = RK4(a,b,n,y0)
%EJEMPLO: RK4(0,7,200,pi^2);

format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);
K1 = 0;
K2 = 0;
K3 = 0;
K4 = 0;

y(1,1) = y0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    K1 = f(x(i,1),y(i,1));
    K2 = f(x(i,1) + (1/2)*h,y(i,1) + (1/2)*h*K1);
    K3 = f(x(i,1) + (1/2)*h,y(i,1) + (1/2)*h*K2);
    K4 = f(x(i,1) + h,y(i,1) + h*K3);
    y(i+1,1) = y(i,1) + h*((1/6)*K1 + (2/6)*K2 + (2/6)*K3 + (1/6)*K4);
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end

