function [x,y] = Heunn_RK3(a,b,n,x0,y0)
%EJEMPLO: [x,y] = Heunn_RK3(0.01,1/2,100,0.0099999166667361106977490148082,0.000110517092175955069979947);

format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);
KX1 = 0; KY1 = 0;
KX2 = 0; KY2 = 0;
KX3 = 0; KY3 = 0;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    KX1 = f(t(i,1),x(i,1),y(i,1));
    KY1 = g(t(i,1),x(i,1),y(i,1));
    KX2 = f(t(i,1) + (1/3)*h,x(i,1) + (1/3)*h*KX1,y(i,1) + (1/3)*h*KY1);
    KY2 = g(t(i,1) + (1/3)*h,x(i,1) + (1/3)*h*KX1,y(i,1) + (1/3)*h*KY1);
    KX3 = f(t(i,1) + (2/3)*h,x(i,1) + (2/3)*h*KX2,y(i,1) + (2/3)*h*KY2);
    KY3 = g(t(i,1) + (2/3)*h,x(i,1) + (2/3)*h*KX2,y(i,1) + (2/3)*h*KY2);
    x(i+1,1) = x(i,1) + h*((1/4)*KX1 + (3/4)*KX3);
    y(i+1,1) = y(i,1) + h*((1/4)*KY1 + (3/4)*KY3);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end