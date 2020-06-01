function [x,y,xprima,yprima] = RKN_simple_orden_5(a,b,n,x0,y0,xprima0,yprima0)
%EJEMPLO: [x,y,xprima,yprima] = RKN_simple_orden_5(0,2.5,100,0,1,0,0);

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
    KX2 = f(t(i,1) + (1/5)*h,x(i,1) + (1/5)*h*xprima(i,1) + (1/50)*h^2*KX1,y(i,1) + (1/5)*h*yprima(i,1) + (1/50)*h^2*KY1);
    KY2 = g(t(i,1) + (1/5)*h,x(i,1) + (1/5)*h*xprima(i,1) + (1/50)*h^2*KX1,y(i,1) + (1/5)*h*yprima(i,1) + (1/50)*h^2*KY1);
    KX3 = f(t(i,1) + (2/3)*h,x(i,1) + (2/3)*h*xprima(i,1) + h^2*(-(1/27)*KX1 + (7/27)*KX2),y(i,1) + (2/3)*h*yprima(i,1) + h^2*(-(1/27)*KY1 + (7/27)*KY2));
    KY3 = g(t(i,1) + (2/3)*h,x(i,1) + (2/3)*h*xprima(i,1) + h^2*(-(1/27)*KX1 + (7/27)*KX2),y(i,1) + (2/3)*h*yprima(i,1) + h^2*(-(1/27)*KY1 + (7/27)*KY2));
    KX4 = f(t(i,1) + h,x(i,1) + h*xprima(i,1) + h^2*((3/10)*KX1 - (2/35)*KX2 + (9/35)*KX3),y(i,1) + h*yprima(i,1) + h^2*((3/10)*KY1 - (2/35)*KY2 + (9/35)*KY3));
    KY4 = g(t(i,1) + h,x(i,1) + h*xprima(i,1) + h^2*((3/10)*KX1 - (2/35)*KX2 + (9/35)*KX3),y(i,1) + h*yprima(i,1) + h^2*((3/10)*KY1 - (2/35)*KY2 + (9/35)*KY3));
    x(i+1,1) = x(i,1) + h*xprima(i,1) + h^2*((14/336)*KX1 + (100/336)*KX2 + (54/336)*KX3);
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((14/336)*KY1 + (100/336)*KY2 + (54/336)*KY3);
    xprima(i+1,1) = xprima(i,1) + h*((14/336)*KX1 + (125/336)*KX2 + (162/336)*KX3 + (35/336)*KX4);
    yprima(i+1,1) = yprima(i,1) + h*((14/336)*KY1 + (125/336)*KY2 + (162/336)*KY3 + (35/336)*KY4);
t(n,1) = a + n*h;
plot(t,x,t,y);
%plot(t,xprima,t,yprima);
grid;
end