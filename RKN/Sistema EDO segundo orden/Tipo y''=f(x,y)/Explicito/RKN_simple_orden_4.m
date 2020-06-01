function [x,y,xprima,yprima] = RKN_simple_orden_4(a,b,n,x0,y0,xprima0,yprima0)
%EJEMPLO: [x,y,xprima,yprima] = RKN_simple_orden_4(0,2.5,100,0,1,0,0);

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

x(1,1) = x0;
y(1,1) = y0;
xprima(1,1) = xprima0;
yprima(1,1) = yprima0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    KX1 = f(t(i,1),x(i,1),y(i,1));
    KY1 = g(t(i,1),x(i,1),y(i,1));
    KX2 = f(t(i,1) + (1/2)*h,x(i,1) + (1/2)*h*xprima(i,1) + (1/8)*h^2*KX1,y(i,1) + (1/2)*h*yprima(i,1) + (1/8)*h^2*KY1);
    KY2 = g(t(i,1) + (1/2)*h,x(i,1) + (1/2)*h*xprima(i,1) + (1/8)*h^2*KX1,y(i,1) + (1/2)*h*yprima(i,1) + (1/8)*h^2*KY1);
    KX3 = f(t(i,1) + h,x(i,1) + h*xprima(i,1) + (1/2)*h^2*KX2,y(i,1) + h*yprima(i,1) + (1/2)*h^2*KY2);
    KY3 = g(t(i,1) + h,x(i,1) + h*xprima(i,1) + (1/2)*h^2*KX2,y(i,1) + h*yprima(i,1) + (1/2)*h^2*KY2);
    x(i+1,1) = x(i,1) + h*xprima(i,1) + h^2*((1/6)*KX1 + (1/3)*KX2);
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((1/6)*KY1 + (1/3)*KY2);
    xprima(i+1,1) = xprima(i,1) + h*((1/6)*KX1 + (4/6)*KX2 + (1/6)*KX3);
    yprima(i+1,1) = yprima(i,1) + h*((1/6)*KY1 + (4/6)*KY2 + (1/6)*KY3);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
%plot(t,xprima,t,yprima);
grid;
end