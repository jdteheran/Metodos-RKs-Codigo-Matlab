syms t;
f = sqrt(t*sin(t));
g = exp(sqrt(t))*tan(t^2);
fplot(f,[0 0.5]);
hold on;
fplot(g,[0 0.5]);
grid;