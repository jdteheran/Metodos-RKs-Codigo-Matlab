function [y yprima] = RKN_direct_orden_3(a,b,n,y0,yprima0)
%EJEMPLO: RKN_direct_orden_3();

format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);
yprima = zeros(n,1);
K1 = 0;
K2 = 0;
K3 = 0;

y(1,1) = y0;
yprima(1,1) = yprima0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    K1 = f(x(i,1),y(i,1),yprima(i,1));
    K2 = f(x(i,1) + (1/2)*h,y(i,1) + (1/2)*h*yprima(i,1) + (1/8)*h^2*K1,yprima(i,1) + (1/2)*h*K1);
    K3 = f(x(i,1) + (3/4)*h,y(i,1) + (3/4)*h*yprima(i,1) + (9/32)*h^2*K2,yprima(i,1) + (3/4)*h*K2);
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((2/9)*K1 + (1/6)*K2 + (1/9)*K3);
    yprima(i+1,1) = yprima(i,1) + h*((2/9)*K1 + (1/3)*K2 + (4/9)*K3);
end
x(n,1) = a + n*h;
plot(x,y);
%plot(x,yprima);
%plot(y,yprima);
grid;
end