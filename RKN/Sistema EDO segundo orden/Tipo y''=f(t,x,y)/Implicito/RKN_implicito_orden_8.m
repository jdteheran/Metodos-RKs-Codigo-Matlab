function [x y xprima yprima] = RKN_implicito_orden_8(a,b,n,x0,y0,xprima0,yprima0)
%EJEMPLO: [x y xprima yprima] = RKN_implicito_orden_8(0,2.5,200,0,1,0,0);

syms kp1x kp2x kp3x kp4x kp1y kp2y kp3y kp4y;
format long;

h = (b-a)/n;
t = zeros(n,1);
x = zeros(n,1);
y = zeros(n,1);
xprima = zeros(n,1);
yprima = zeros(n,1);

KP1X = 0;    KP1Y = 0;
KP2X = 0;    KP2Y = 0;
KP3X = 0;    KP3Y = 0;
KP4X = 0;    KP4Y = 0;

W1 = (1/8)-(sqrt(30)/144);              W1P = (1/8)+(sqrt(30)/144);
W2 = (1/2)*sqrt((15+2*sqrt(30))/35);    W2P = (1/2)*sqrt((15-2*sqrt(30))/35);
W3 = W2*((1/6)+(sqrt(30)/24));          W3P = W2P*((1/6)-(sqrt(30)/24));
W4 = W2*((1/21)+((5*sqrt(30))/168));    W4P = W2P*((1/21)-(5*sqrt(30)/168));
W5 = W2-2*W3;                           W5P = W2P - 2*W3P;

AR11 = ((2*W2 + 1)^2*(16*W2^3 - 16*W2^2 - 80*W2*W2P^2 + 12*W2 + 20*W2P^2 - 3))/(960*W2*(W2^2 - W2P^2));   AR12 = -(W1P*(2*W2P + 1)^3*(50*W2*W2P - 21*W2P - 15*W2 + 14*W2P^2 + 6))/(1920*W1*W2*(W2^2 - W2P^2));             AR13 = (W1P*(2*W2P - 1)^3*(21*W2P - 15*W2 - 50*W2*W2P + 14*W2P^2 + 6))/(1920*W1*W2*(W2^2 - W2P^2));              AR14 = ((2*W2 - 1)^3*(2*W2^2 + 3*W2 - 20*W2P^2 + 3))/(960*W2*(W2^2 - W2P^2));
AR21 = (W1*(2*W2 + 1)^3*(50*W2*W2P - 15*W2P - 21*W2 + 14*W2^2 + 6))/(1920*W1P*W2P*(W2^2 - W2P^2));        AR22 = -((2*W2P + 1)^2*(- 80*W2^2*W2P + 20*W2^2 + 16*W2P^3 - 16*W2P^2 + 12*W2P - 3))/(960*W2P*(W2^2 - W2P^2));   AR23 = -((2*W2P - 1)^3*(- 20*W2^2 + 2*W2P^2 + 3*W2P + 3))/(960*W2P*(W2^2 - W2P^2));                              AR24 = -(W1*(2*W2 - 1)^3*(21*W2 - 15*W2P - 50*W2*W2P + 14*W2^2 + 6))/(1920*W1P*W2P*(W2^2 - W2P^2));
AR31 = -(W1*(2*W2 + 1)^3*(15*W2P - 21*W2 - 50*W2*W2P + 14*W2^2 + 6))/(1920*W1P*W2P*(W2^2 - W2P^2));       AR32 = ((2*W2P + 1)^3*(20*W2^2 - 2*W2P^2 + 3*W2P - 3))/(960*W2P*(W2^2 - W2P^2));                                 AR33 = -((2*W2P - 1)^2*(- 80*W2^2*W2P - 20*W2^2 + 16*W2P^3 + 16*W2P^2 + 12*W2P + 3))/(960*W2P*(W2^2 - W2P^2));   AR34 = (W1*(2*W2 - 1)^3*(21*W2 + 15*W2P + 50*W2*W2P + 14*W2^2 + 6))/(1920*W1P*W2P*(W2^2 - W2P^2));
AR41 = -((2*W2 + 1)^3*(- 2*W2^2 + 3*W2 + 20*W2P^2 - 3))/(960*W2*(W2^2 - W2P^2));                          AR42 = (W1P*(2*W2P + 1)^3*(15*W2 - 21*W2P - 50*W2*W2P + 14*W2P^2 + 6))/(1920*W1*W2*(W2^2 - W2P^2));              AR43 = -(W1P*(2*W2P - 1)^3*(15*W2 + 21*W2P + 50*W2*W2P + 14*W2P^2 + 6))/(1920*W1*W2*(W2^2 - W2P^2));             AR44 = ((2*W2 - 1)^2*(16*W2^3 + 16*W2^2 - 80*W2*W2P^2 + 12*W2 - 20*W2P^2 + 3))/(960*W2*(W2^2 - W2P^2));

A11 = W1;          A12 = W1P-W3+W4P;   A13 = W1P-W3-W4P;   A14 = W1-W5;
A21 = W1-W3P+W4;   A22 = W1P;          A23 = W1P-W5P;      A24 = W1-W5P-W4;
A31 = W1+W3P+W4;   A32 = W1P+W5P;      A33 = W1P;          A34 = W1+W5P-W4;
A41 = W1+W5;       A42 = W1P+W3+W4P;   A43 = W1P+W3-W4P;   A44 = W1;

C1 = (1/2)-W2; C2 = (1/2)-W2P; C3 = (1/2)+W2P; C4 = (1/2)+W2;

B1 = 2*W1; B2 = 2*W1P; B3 = 2*W1P; B4 = 2*W1;

BR1 = B1*(1-C1); BR2 = B2*(1-C2); BR3 = B3*(1-C3); BR4 = B4*(1-C4);

x(1,1) = x0;
y(1,1) = y0;
xprima(1,1) = xprima0;
yprima(1,1) = yprima0;
for i = 1:n-1
    t(i,1) = a + (i-1)*h;
   
    F = [f(t(i,1) + C1*h, x(i,1) + C1*h*xprima(i,1) + h^2*(AR11*kp1x + AR12*kp2x + AR13*kp3x + AR14*kp4x), y(i,1) + C1*h*yprima(i,1) + h^2*(AR11*kp1y + AR12*kp2y + AR13*kp3y + AR14*kp4y), xprima(i,1) + h*(A11*kp1x + A12*kp2x + A13*kp3x + A14*kp4x), yprima(i,1) + h*(A11*kp1y + A12*kp2y + A13*kp3y + A14*kp4y)) - kp1x;
         g(t(i,1) + C2*h, x(i,1) + C2*h*xprima(i,1) + h^2*(AR21*kp1x + AR22*kp2x + AR23*kp3x + AR24*kp4x), y(i,1) + C2*h*yprima(i,1) + h^2*(AR21*kp1y + AR22*kp2y + AR23*kp3y + AR24*kp4y), xprima(i,1) + h*(A21*kp1x + A22*kp2x + A23*kp3x + A24*kp4x), yprima(i,1) + h*(A21*kp1y + A22*kp2y + A23*kp3y + A24*kp4y)) - kp1y;
         f(t(i,1) + C3*h, x(i,1) + C3*h*xprima(i,1) + h^2*(AR31*kp1x + AR32*kp2x + AR33*kp3x + AR34*kp4x), y(i,1) + C3*h*yprima(i,1) + h^2*(AR31*kp1y + AR32*kp2y + AR33*kp3y + AR34*kp4y), xprima(i,1) + h*(A31*kp1x + A32*kp2x + A33*kp3x + A34*kp4x), yprima(i,1) + h*(A31*kp1y + A32*kp2y + A33*kp3y + A34*kp4y)) - kp2x;
         g(t(i,1) + C4*h, x(i,1) + C4*h*xprima(i,1) + h^2*(AR41*kp1x + AR42*kp2x + AR43*kp3x + AR44*kp4x), y(i,1) + C4*h*yprima(i,1) + h^2*(AR41*kp1y + AR42*kp2y + AR43*kp3y + AR44*kp4y), xprima(i,1) + h*(A41*kp1x + A42*kp2x + A43*kp3x + A44*kp4x), yprima(i,1) + h*(A41*kp1y + A42*kp2y + A43*kp3y + A44*kp4y)) - kp2y;
         f(t(i,1) + C1*h, x(i,1) + C1*h*xprima(i,1) + h^2*(AR11*kp1x + AR12*kp2x + AR13*kp3x + AR14*kp4x), y(i,1) + C1*h*yprima(i,1) + h^2*(AR11*kp1y + AR12*kp2y + AR13*kp3y + AR14*kp4y), xprima(i,1) + h*(A11*kp1x + A12*kp2x + A13*kp3x + A14*kp4x), yprima(i,1) + h*(A11*kp1y + A12*kp2y + A13*kp3y + A14*kp4y)) - kp3x;
         g(t(i,1) + C2*h, x(i,1) + C2*h*xprima(i,1) + h^2*(AR21*kp1x + AR22*kp2x + AR23*kp3x + AR24*kp4x), y(i,1) + C2*h*yprima(i,1) + h^2*(AR21*kp1y + AR22*kp2y + AR23*kp3y + AR24*kp4y), xprima(i,1) + h*(A21*kp1x + A22*kp2x + A23*kp3x + A24*kp4x), yprima(i,1) + h*(A21*kp1y + A22*kp2y + A23*kp3y + A24*kp4y)) - kp3y;
         f(t(i,1) + C3*h, x(i,1) + C3*h*xprima(i,1) + h^2*(AR31*kp1x + AR32*kp2x + AR33*kp3x + AR34*kp4x), y(i,1) + C3*h*yprima(i,1) + h^2*(AR31*kp1y + AR32*kp2y + AR33*kp3y + AR34*kp4y), xprima(i,1) + h*(A31*kp1x + A32*kp2x + A33*kp3x + A34*kp4x), yprima(i,1) + h*(A31*kp1y + A32*kp2y + A33*kp3y + A34*kp4y)) - kp4x;
         g(t(i,1) + C4*h, x(i,1) + C4*h*xprima(i,1) + h^2*(AR41*kp1x + AR42*kp2x + AR43*kp3x + AR44*kp4x), y(i,1) + C4*h*yprima(i,1) + h^2*(AR41*kp1y + AR42*kp2y + AR43*kp3y + AR44*kp4y), xprima(i,1) + h*(A41*kp1x + A42*kp2x + A43*kp3x + A44*kp4x), yprima(i,1) + h*(A41*kp1y + A42*kp2y + A43*kp3y + A44*kp4y)) - kp4y];
     
    K = Newton_Raphson_Multivariable(F, [KP1X;KP2X;KP3X;KP4X;KP1Y;KP2Y;KP3Y;KP4Y], 10^-8);
    
    %TENER EN CUENTA COMO EL ALGORITMO DE NEWTON_RAPHSON TE DEVUELVE LAS
    %VARIABLES
    KP1X = K(1,1);
    KP1Y = K(2,1);
    KP2X = K(3,1);
    KP2Y = K(4,1);
    KP3X = K(5,1);
    KP3Y = K(6,1);
    KP4X = K(7,1);
    KP4Y = K(8,1);
    
    x(i+1,1) = x(i,1) + h*xprima(i,1) + h^2*(BR1*KP1X + BR2*KP2X + BR3*KP3X + BR4*KP4X);
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*(BR1*KP1Y + BR2*KP2Y + BR3*KP3Y + BR4*KP4Y);
    xprima(i+1,1) = xprima(i,1) + h*(B1*KP1X + B2*KP2X + B3*KP3X + B4*KP4X);
    yprima(i+1,1) = yprima(i,1) + h*(B1*KP1Y + B2*KP2Y + B3*KP3Y + B4*KP4Y);
end
t(n,1) = a + n*h;
plot(t,x,t,y);
%plot(x,yprima);
%plot(y,yprima);
grid;
end