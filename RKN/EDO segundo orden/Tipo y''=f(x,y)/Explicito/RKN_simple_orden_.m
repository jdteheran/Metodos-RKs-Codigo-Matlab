function [y yprima] = RKN_simple_orden_(a,b,n,y0,yprima0)
%EJEMPLO: [y yprima] = RKN_simple_orden_(0,10,200,0,1);

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
    K2 = f(x(i,1) + (2/5)*h,y(i,1) + (2/5)*h*yprima(i,1) + (2/25)*h^2*K1);
    K3 = f(x(i,1) + (2/3)*h,y(i,1) + (2/3)*h*yprima(i,1) + h^2*(2/9)*K1);
    K4 = f(x(i,1) + (4/5)*h,y(i,1) + (4/5)*h*yprima(i,1) + h^2*((4/25)*K1 + (4/25)*K2));
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((23/192)*K1 + (75/192)*K2 - (27/192)*K3 + (25/192)*K4);
    yprima(i+1,1) = yprima(i,1) + h*((23/192)*K1 + (125/192)*K2 - (81/192)*K3 + (125/192)*K4);
end
x(n,1) = a + n*h;
%plot(x,y);
%plot(x,yprima);
plot(y,yprima);
grid;
end