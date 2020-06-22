function [gval] = g(t,x,y,xprima,yprima)
gval = - 2*x - 4*t^2*y/(x^2+y^2);
end

