function [gval] = g(t,x,y,xprima,yprima)
gval = exp(xprima + 2*x)*((sec(t))^2+log(y^2));
end

