function [x,y] = Kuntzmann_Butcher_RK4_I(a,b,n,x0,y0)
%EJEMPLO: [x,y] = Kuntzmann_Butcher_RK4_I(0.01,1/2,20,0.0099999166667361106977490148082,0.000110517092175955069979947);

syms kx1 kx2 kx3 kx4 ky1 ky2 ky3 ky4;
format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);

KX1 = 1; KY1 = 1;
KX2 = 1; KY2 = 1;
KX3 = 1; KY3 = 1;
KX4 = 1; KY4 = 1;

W1 = (1/8)-(sqrt(30)/144);              W1P = (1/8)+(sqrt(30)/144);
W2 = (1/2)*sqrt((15+2*sqrt(30))/35);    W2P = (1/2)*sqrt((15-2*sqrt(30))/35);
W3 = W2*((1/6)+(sqrt(30)/24));          W3P = W2P*((1/6)-(sqrt(30)/24));
W4 = W2*((1/21)+((5*sqrt(30))/168));    W4P = W2P*((1/21)-(5*sqrt(30)/168));
W5 = W2-2*W3;                           W5P = W2P - 2*W3P;

x(1,1) = x0;
y(1,1) = y0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
    
    F = [f(t(i,1) + ((1/2)-W2)*h, x(i,1) + W1*h*kx1 + (W1P-W3+W4P)*h*kx2 + (W1P-W3-W4P)*h*kx3 + (W1-W5)*h*kx4, y(i,1) + W1*h*ky1 + (W1P-W3+W4P)*h*ky2 + (W1P-W3-W4P)*h*ky3 + (W1-W5)*h*ky4) - kx1 ;
         g(t(i,1) + ((1/2)-W2)*h, x(i,1) + W1*h*kx1 + (W1P-W3+W4P)*h*kx2 + (W1P-W3-W4P)*h*kx3 + (W1-W5)*h*kx4, y(i,1) + W1*h*ky1 + (W1P-W3+W4P)*h*ky2 + (W1P-W3-W4P)*h*ky3 + (W1-W5)*h*ky4) - ky1 ;
         f(t(i,1) + ((1/2)-W2P)*h, x(i,1) + (W1-W3P+W4)*h*kx1 + (W1P)*h*kx2 + (W1P-W5P)*h*kx3 + (W1-W5P-W4)*h*kx4, y(i,1) + (W1-W3P+W4)*h*ky1 + (W1P)*h*ky2 + (W1P-W5P)*h*ky3 + (W1-W5P-W4)*h*ky4) - kx2 ;
         g(t(i,1) + ((1/2)-W2P)*h, x(i,1) + (W1-W3P+W4)*h*kx1 + (W1P)*h*kx2 + (W1P-W5P)*h*kx3 + (W1-W5P-W4)*h*kx4, y(i,1) + (W1-W3P+W4)*h*ky1 + (W1P)*h*ky2 + (W1P-W5P)*h*ky3 + (W1-W5P-W4)*h*ky4) - ky2 ;
         f(t(i,1) + ((1/2)+W2P)*h, x(i,1) + (W1+W3P+W4)*h*kx1 + (W1P+W5P)*h*kx2 + (W1P)*h*kx3 + (W1+W5P-W4)*h*kx4, y(i,1) + (W1+W3P+W4)*h*ky1 + (W1P+W5P)*h*ky2 + (W1P)*h*ky3 + (W1+W5P-W4)*h*ky4) - kx3 ;
         g(t(i,1) + ((1/2)+W2P)*h, x(i,1) + (W1+W3P+W4)*h*kx1 + (W1P+W5P)*h*kx2 + (W1P)*h*kx3 + (W1+W5P-W4)*h*kx4, y(i,1) + (W1+W3P+W4)*h*ky1 + (W1P+W5P)*h*ky2 + (W1P)*h*ky3 + (W1+W5P-W4)*h*ky4) - ky3 ;
         f(t(i,1) + ((1/2)+W2)*h, x(i,1) + (W1+W5)*h*kx1 + (W1P+W3+W4P)*h*kx2 + (W1P+W3-W4P)*h*kx3 + (W1)*h*kx4, y(i,1) + (W1+W5)*h*ky1 + (W1P+W3+W4P)*h*ky2 + (W1P+W3-W4P)*h*ky3 + (W1)*h*ky4) - kx4 ;
         g(t(i,1) + ((1/2)+W2)*h, x(i,1) + (W1+W5)*h*kx1 + (W1P+W3+W4P)*h*kx2 + (W1P+W3-W4P)*h*kx3 + (W1)*h*kx4, y(i,1) + (W1+W5)*h*ky1 + (W1P+W3+W4P)*h*ky2 + (W1P+W3-W4P)*h*ky3 + (W1)*h*ky4) - ky4];
     
    K = Newton_Raphson_Multivariable(F, [KX1;KX2;KX3;KX4;KY1;KY2;KY3;KY4], 10^-8);
    
    KX1 = K(1,1);
    KX2 = K(2,1);
    KX3 = K(3,1);
    KX4 = K(4,1);
    KY1 = K(5,1);    
    KY2 = K(6,1);
    KY3 = K(7,1);
    KY4 = K(8,1);
    
    x(i+1,1) = x(i,1) + (2*h)*(W1*KX1+W1P*KX2+W1P*KX3+W1*KX4);
    y(i+1,1) = y(i,1) + (2*h)*(W1*KY1+W1P*KY2+W1P*KY3+W1*KY4);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
grid;
end