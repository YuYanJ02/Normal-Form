% -------------------------------------------------------------------------
% 三体动力学方程
%
% Copyright(C) 2015/07/25 by Chen Zhang,
% School of Astronautics, Beihang University
% chenzhang.buaa@gmail.com
% -------------------------------------------------------------------------
function xdot = fun_crtbp_(t , x , auxdata)
mu = auxdata.mu;

d2 = (x(1) + mu)^2 + x(2)^2 + x(3)^2;
r2 = (x(1) - (1 - mu))^2 + x(2)^2 + x(3)^2;
d3 = d2^1.5;
r3 = r2^1.5;

Ux = x(1)  - (1 - mu) * (x(1) + mu) / d3 - mu * (x(1) - (1 - mu)) / r3;
Uy = x(2) - (1 - mu) * x(2) / d3 - mu*x(2) / r3;
Uz = -(1 - mu) * x(3) / d3 - mu*x(3) / r3;

xdot = [x(4);
    x(5);
    x(6);
    Ux + 2*x(5);
    Uy - 2*x(4);
    Uz];
end
% -------------------------------------------------------------------------