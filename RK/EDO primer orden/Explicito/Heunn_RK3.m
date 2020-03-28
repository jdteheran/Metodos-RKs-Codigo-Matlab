function [y] = Heunn_RK3(a,b,n,y0)
%EJEMPLO: Heunn_RK3(0,7,200,pi^2);

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
    K1 = f(x(i,1),y(i,1));
    K2 = f(x(i,1) + (1/3)*h,y(i,1) + (1/3)*h*K1);
    K3 = f(x(i,1) + (2/3)*h,y(i,1) + (1/3)*h*K2);
    y(i+1,1) = y(i,1) + h*((1/4)*K1 + (3/4)*K3);
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end

