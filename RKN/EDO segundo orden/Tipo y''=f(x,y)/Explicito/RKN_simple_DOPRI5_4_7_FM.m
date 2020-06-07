function [y yprima] = RKN_simple_DOPRI5_4_7_FM(a,b,n,y0,yprima0)
%EJEMPLO: [y yprima] = RKN_simple_DOPRI5_4_7_FM(0,10,200,0,1);

format long;

h = (b-a)/n;
x = zeros(n,1);
y = zeros(n,1);
yprima = zeros(n,1);
K1 = 0;
K2 = 0;
K3 = 0;
K4 = 0;
K5 = 0;
K6 = 0;
K7 = 0;

y(1,1) = y0;
yprima(1,1) = yprima0;
for i = 1:n-1
    x(i,1) = a + (i-1)*h;
    K1 = f(x(i,1),y(i,1));
    K2 = f(x(i,1) + (1/5)*h,y(i,1) + (1/5)*h*yprima(i,1) + (1/5)*h^2*K1);
    K3 = f(x(i,1) + (3/10)*h,y(i,1) + (3/10)*h*yprima(i,1) + h^2*((3/40)*K1 + (9/40)*K2));
    K4 = f(x(i,1) + (4/5)*h,y(i,1) + (4/5)*h*yprima(i,1) + h^2*((44/45)*K1 - (56/15)*K2 + (32/9)*K3));
    K5 = f(x(i,1) + (8/9)*h,y(i,1) + (8/9)*h*yprima(i,1) + h^2*((19372/6561)*K1 - (25360/2187)*K2 + (64448/6561)*K3 - (212/729)*K4));
    K6 = f(x(i,1) + h,y(i,1) + h*yprima(i,1) + h^2*((9017/3168)*K1 - (355/33)*K2 + (46732/5247)*K3 + (49/176)*K4 - (5103/18656)*K5));
    K7 = f(x(i,1) + h,y(i,1) + h*yprima(i,1) + h^2*((35/384)*K1 + (500/1113)*K3 + (125/192)*K4 - (2187/6784)*K5 + (11/84)*K6));
    y(i+1,1) = y(i,1) + h*yprima(i,1) + h^2*((35/384)*K1 + (500/1113)*K3 + (125/192)*K4 - (2187/6784)*K5 + (11/84)*K6);
    yprima(i+1,1) = yprima(i,1) + h*((5179/57600)*K1 + (7571/16695)*K3 + (393/640)*K4 - (92097/339200)*K5 + (187/2100)*K6 + (1/40)*K7);
end
x(n,1) = a + n*h;
%plot(x,y);
%plot(x,yprima);
plot(y,yprima);
grid;
end