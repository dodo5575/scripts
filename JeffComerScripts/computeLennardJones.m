% computeLennardJones.m

function v = computeLennardJones(r, eps0, eps1, rad0, rad1)

eps = sqrt(eps0*eps1);
radius = rad0 + rad1;

radius6 = radius*radius*radius*radius*radius*radius;
r6 = r.*r.*r.*r.*r.*r;

v = eps*(radius6*radius6./(r6.*r6) - 2*radius6./r6);


