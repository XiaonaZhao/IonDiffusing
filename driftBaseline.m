function deBaseline = driftBaseline(independentVariable, dependentVariable)

plot(independentVariable, dependentVariable)
hold on

[x, y] = ginput;
% ginput gathers an unlimited number of points until you press the RETURN
% key.
p = polyfit(x, y, 1);
f = polyval(p, independentVariable);

deBaseline = dependentVariable - f;

plot(independentVariable, deBaseline);
legend('data', 'linear fit')
hold off
