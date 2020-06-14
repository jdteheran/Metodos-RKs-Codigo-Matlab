syms t;
f = sin(t^2);
g = cos(t^2);
fplot(f,[0 2.5]);
hold on;
fplot(g,[0 2.5]);
grid;