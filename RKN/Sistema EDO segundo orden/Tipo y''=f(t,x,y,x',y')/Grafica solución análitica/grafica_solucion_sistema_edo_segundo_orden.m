syms t;
f = log(sec(t));
g = exp(tan(t));
fplot(f,[0 1]);
hold on;
fplot(g,[0 1]);
grid;