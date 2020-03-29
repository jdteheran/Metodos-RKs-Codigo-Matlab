function [x,y] = Euler_RK1(a,b,n,x0,y0)
%EJEMPLO: [x,y] = Euler_RK1(0.01,1/2,100,0.0099999166667361106977490148082,0.0001105170921759550699799478876);

format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    x(i+1,1) = x(i,1) + h*f(t(i,1),x(i,1),y(i,1));
    y(i+1,1) = y(i,1) + h*g(t(i,1),x(i,1),y(i,1));
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end

