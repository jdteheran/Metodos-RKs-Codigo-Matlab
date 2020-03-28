function [fval] = f(t,x,y)
fval = cot(t)*x*y/(2*(exp(sqrt(t)))*tan(t^2))+(1/(2*t))*x;
end