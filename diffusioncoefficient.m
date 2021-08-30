% diffusion coefficient
distance = 71;
% t = (50-19)/100;
t = (375-350)/100;
c = 308.8;
cmax = 710.1;

L = distance*250e-9; % Unit: m

D = (L^2)/(-4*t*log(c/cmax))