function [gval] = g(t,x,y)
gval = 2*t*y^2/exp(sqrt(t))+y*sqrt(sin(t))/(2*x)+2*csc(t)*exp(sqrt(t))*x^2;
end

