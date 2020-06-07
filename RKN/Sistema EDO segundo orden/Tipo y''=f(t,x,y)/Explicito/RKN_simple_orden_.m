function [x,y,xprima,yprima] = RKN_simple_orden_(a,b,n,x0,y0,xprima0,yprima0)
%EJEMPLO: [x,y,xprima,yprima] = RKN_simple_orden_(0,2.5,200,0,1,0,0);

format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);
xprima = zeros(n,1);
yprima = zeros(n,1);
KX1 = 0; KY1 = 0;
KX2 = 0; KY2 = 0;
KX3 = 0; KY3 = 0;
KX4 = 0; KY4 = 0;

x(1,1) = x0;
y(1,1) = y0;
xprima(1,1) = xprima0;
yprima(1,1) = yprima0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    KX1 = f(t(i,1),x(i,1),y(i,1));
    KY1 = g(t(i,1),x(i,1),y(i,1));
    KX2 = f(t(i,1) + (2/5)*h,x(i,1) + (2/5)*h*xprima(i,1) + (2/25)*h^2*KX1,y(i,1) + (2/5)*h*yprima(i,1) + (2/25)*h^2*KY1);
    KY2 = g(t(i,1) + (2/5)*h,x(i,1) + (2/5)*h*xprima(i,1) + (2/25)*h^2*KX1,y(i,1) + (2/5)*h*yprima(i,1) + (2/25)*h^2*KY1);
    KX3 = f(t(i,1) + (2/3)*h,x(i,1) + (2/3)*h*xprima(i,1) + h^2*(2/9)*KX1,y(i,1) + (2/3)*h*yprima(i,1) + h^2*(2/9)*KY1);
    KY3 = g(t(i,1) + (2/3)*h,x(i,1) + (2/3)*h*xprima(i,1) + h^2*(2/9)*KX1,y(i,1) + (2/3)*h*yprima(i,1) + h^2*(2/9)*KY1);
    KX4 = f(t(i,1) + (4/5)*h,x(i,1) + (4/5)*h*xprima(i,1) + h^2*((4/25)*KX1 + (4/25)*KX2),y(i,1) + (4/5)*h*yprima(i,1) + h^2*((4/25)*KY1 + (4/25)*KY2));
    KY4 = g(t(i,1) + (4/5)*h,x(i,1) + (4/5)*h*xprima(i,1) + h^2*((4/25)*KX1 + (4/25)*KX2),y(i,1) + (4/5)*h*yprima(i,1) + h^2*((4/25)*KY1 + (4/25)*KY2));
    x(i+1,1) = x(i,1) + h*xprima(i,1) + h^2*((23/192)*KX1 + (75/192)*KX2 - (27/192)*KX3 + (25/192)*KX4);
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((23/192)*KY1 + (75/192)*KY2 - (27/192)*KY3 + (25/192)*KY4);
    xprima(i+1,1) = xprima(i,1) + h*((23/192)*KX1 + (125/192)*KX2 - (81/192)*KX3 + (125/192)*KX4);
    yprima(i+1,1) = yprima(i,1) + h*((23/192)*KY1 + (125/192)*KY2 - (81/192)*KY3 + (125/192)*KY4);
t(n,1) = a + n*h;
%plot(t,x,t,y);
%plot(t,xprima,t,yprima);
plot(x,xprima,y,yprima);
grid;
end