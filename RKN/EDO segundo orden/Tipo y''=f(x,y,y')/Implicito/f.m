function [fval] = f(x,y,yprima)
fval = yprima/(x-y^2);
end