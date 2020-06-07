function [y yprima] = RKN_orden_4(a,b,n,y0,yprima0)
%EJEMPLO: [y yprima] = RKN_orden_4(0,25,200,1,1);

format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);
yprima = zeros(n,1);
K1 = 0;
K2 = 0;
K3 = 0;
K4 = 0;

y(1,1) = y0;
yprima(1,1) = yprima0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    K1 = f(x(i,1),y(i,1),yprima(i,1));
    K2 = f(x(i,1) + (1/2)*h,y(i,1) + (1/2)*h*yprima(i,1) + (1/8)*h^2*K1,yprima(i,1) + (1/2)*h*K1);
    K3 = f(x(i,1) + (1/2)*h,y(i,1) + (1/2)*h*yprima(i,1) + (1/8)*h^2*K1,yprima(i,1) + (1/2)*h*K2);
    K4 = f(x(i,1) + h,y(i,1) + h*yprima(i,1) + (1/2)*h^2*K3,yprima(i,1) + h*K3);
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((1/6)*K1 + (1/6)*K2 + (1/6)*K3);
    yprima(i+1,1) = yprima(i,1) + h*((1/6)*K1 + (2/6)*K2 + (2/6)*K3 + (1/6)*K4);
end
x(n,1) = a + n*h;
%plot(x,y);
%plot(x,yprima);
plot(y,yprima);
grid;
end