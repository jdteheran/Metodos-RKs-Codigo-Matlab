function [fval] = f(t,x,y)
fval = 2*y/(x^2+y^2) - 4*t^2*x;
end