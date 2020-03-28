function [y] = Kuntzmann_Butcher_RK4_I(a,b,n,y0)
%EJEMPLO: Kuntzmann_Butcher_RK4_I(0,7,50,pi^2)

syms k1 k2 k3 k4;
format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);

K1 = 0;
K2 = 0;
K3 = 0;
K4 = 0;

W1 = (1/8)-(sqrt(30)/144);              W1P = (1/8)+(sqrt(30)/144);
W2 = (1/2)*sqrt((15+2*sqrt(30))/35);    W2P = (1/2)*sqrt((15-2*sqrt(30))/35);
W3 = W2*((1/6)+(sqrt(30)/24));          W3P = W2P*((1/6)-(sqrt(30)/24));
W4 = W2*((1/21)+((5*sqrt(30))/168));    W4P = W2P*((1/21)-(5*sqrt(30)/168));
W5 = W2-2*W3;                           W5P = W2P - 2*W3P;

y(1,1) = y0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    
    F = [f(x(i,1) + ((1/2)-W2)*h, y(i,1) + W1*h*k1 + (W1P-W3+W4P)*h*k2 + (W1P-W3-W4P)*h*k3 + (W1-W5)*h*k4) - k1 ;
         f(x(i,1) + ((1/2)-W2P)*h, y(i,1) + (W1-W3P+W4)*h*k1 + (W1P)*h*k2 + (W1P-W5P)*h*k3 + (W1-W5P-W4)*h*k4) - k2 ;
         f(x(i,1) + ((1/2)+W2P)*h, y(i,1) + (W1+W3P+W4)*h*k1 + (W1P+W5P)*h*k2 + (W1P)*h*k3 + (W1+W5P-W4)*h*k4) - k3 ;
         f(x(i,1) + ((1/2)+W2)*h, y(i,1) + (W1+W5)*h*k1 + (W1P+W3+W4P)*h*k2 + (W1P+W3-W4P)*h*k3 + (W1)*h*k4) - k4];
    K = Newton_Raphson_Multivariable(F, [K1;K2;K3;K4], 10^-8);
    K1 = K(1,1);
    K2 = K(2,1);
    K3 = K(3,1);
    K4 = K(4,1);
    
    y(i+1,1) = y(i,1) + (2*h)*(W1*K1+W1P*K2+W1P*K3+W1*K4);
end
x(n,1) = a + n*h;
plot(x,y);
grid;
end