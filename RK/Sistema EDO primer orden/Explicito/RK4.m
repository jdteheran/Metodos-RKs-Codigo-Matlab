function [x,y] = RK4(a,b,n,x0,y0)
%EJEMPLO: [x,y] = RK4(0.01,1/2,100,0.0099999166667361106977490148082,0.0001105170921759550699799478876);

format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);
KX1 = 0; KY1 = 0;
KX2 = 0; KY2 = 0;
KX3 = 0; KY3 = 0;
KX4 = 0; KY4 = 0;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    KX1 = f(t(i,1),x(i,1),y(i,1));
    KY1 = g(t(i,1),x(i,1),y(i,1));
    KX2 = f(t(i,1) + (1/2)*h,x(i,1) + (1/2)*h*KX1,y(i,1) + (1/2)*h*KY1);
    KY2 = g(t(i,1) + (1/2)*h,x(i,1) + (1/2)*h*KX1,y(i,1) + (1/2)*h*KY1);
    KX3 = f(t(i,1) + (1/2)*h,x(i,1) + (1/2)*h*KX2,y(i,1) + (1/2)*h*KY2);
    KY3 = g(t(i,1) + (1/2)*h,x(i,1) + (1/2)*h*KX2,y(i,1) + (1/2)*h*KY2);
    KX4 = f(t(i,1) + h,x(i,1) + h*KX3,y(i,1) + h*KY3);
    KY4 = g(t(i,1) + h,x(i,1) + h*KX3,y(i,1) + h*KY3);
    x(i+1,1) = x(i,1) + h*((1/6)*KX1 + (2/6)*KX2 + (2/6)*KX3 + (1/6)*KX4);
    y(i+1,1) = y(i,1) + h*((1/6)*KY1 + (2/6)*KY2 + (2/6)*KY3 + (1/6)*KY4);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end