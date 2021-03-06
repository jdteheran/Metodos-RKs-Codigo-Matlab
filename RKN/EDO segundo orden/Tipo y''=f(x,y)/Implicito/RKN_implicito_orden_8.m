function [y yprima] = RKN_implicito_orden_8(a,b,n,y0,yprima0)
%EJEMPLO: [y yprima] = RKN_implicito_orden_8(0,10,200,0,1);

syms kp1 kp2 kp3 kp4;
format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);
yprima = zeros(n,1);

KP1 = 0;
KP2 = 0;
KP3 = 0;
KP4 = 0;

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


y(1,1) = y0;
yprima(1,1) = yprima0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
   
    F = [f(x(i,1) + C1*h, y(i,1) + C1*h*yprima(i,1) + h^2*(AR11*kp1 + AR12*kp2 + AR13*kp3 + AR14*kp4), yprima(i,1) + h*(A11*kp1 + A12*kp2 + A13*kp3 + A14*kp4)) - kp1;
         f(x(i,1) + C2*h, y(i,1) + C2*h*yprima(i,1) + h^2*(AR21*kp1 + AR22*kp2 + AR23*kp3 + AR24*kp4), yprima(i,1) + h*(A21*kp1 + A22*kp2 + A23*kp3 + A24*kp4)) - kp2;
         f(x(i,1) + C3*h, y(i,1) + C3*h*yprima(i,1) + h^2*(AR31*kp1 + AR32*kp2 + AR33*kp3 + AR34*kp4), yprima(i,1) + h*(A31*kp1 + A32*kp2 + A33*kp3 + A34*kp4)) - kp3;
         f(x(i,1) + C4*h, y(i,1) + C4*h*yprima(i,1) + h^2*(AR41*kp1 + AR42*kp2 + AR43*kp3 + AR44*kp4), yprima(i,1) + h*(A41*kp1 + A42*kp2 + A43*kp3 + A44*kp4)) - kp4];
    K = Newton_Raphson_Multivariable(F, [KP1;KP2;KP3;KP4], 10^-8);
    
    %TENER EN CUENTA COMO EL ALGORITMO DE NEWTON_RAPHSON TE DEVUELVE LAS
    %VARIABLES
    KP1 = K(1,1);
    KP2 = K(2,1);
    KP3 = K(3,1);
    KP4 = K(4,1);
    
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*(BR1*KP1 + BR2*KP2 + BR3*KP3 + BR4*KP4);
    yprima(i+1,1) = yprima(i,1) + h*(B1*KP1 + B2*KP2 + B3*KP3 + B4*KP4);
end
x(n,1) = a + n*h;
plot(x,y);
%plot(x,yprima);
%plot(y,yprima);
grid;
end