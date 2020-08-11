function [y yprima] = RKN_simple_orden_5(a,b,n,y0,yprima0)
%EJEMPLO: [y yprima] = RKN_simple_orden_5(0,10,200,0,1);

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
    K1 = f(x(i,1),y(i,1));
    K2 = f(x(i,1) + (1/5)*h,y(i,1) + (1/5)*h*yprima(i,1) + (1/50)*h^2*K1);
    K3 = f(x(i,1) + (2/3)*h,y(i,1) + (2/3)*h*yprima(i,1) + h^2*(-(1/27)*K1 + (7/27)*K2));
    K4 = f(x(i,1) + h,y(i,1) + h*yprima(i,1) + h^2*((3/10)*K1 - (2/35)*K2 + (9/35)*K3));
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((14/336)*K1 + (100/336)*K2 + (54/336)*K3);
    yprima(i+1,1) = yprima(i,1) + h*((14/336)*K1 + (125/336)*K2 + (162/336)*K3 + (35/336)*K4);
end
x(n,1) = a + n*h;
plot(x,y);
%plot(x,yprima);
%plot(y,yprima);
grid;
end