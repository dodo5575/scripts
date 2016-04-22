% computeHardcore.m

function v = computeHardcore(d, eps0, eps1, rad0, rad1)

eps = sqrt(eps0*eps1);
radius = rad0 + rad1;

radius6 = radius*radius*radius*radius*radius*radius;
d6 = d.*d.*d.*d.*d.*d;

v = eps*(d6<radius6).*((radius6*radius6./(d6.*d6) - 2*radius6./d6) + 1);


