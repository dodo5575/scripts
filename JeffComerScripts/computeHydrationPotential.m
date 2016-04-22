% computeHydrationPotential.m

function v = computeHydrationPotential(r, c0, c1, c2, c3, c4)

d = c1./r;
d6 = d.*d.*d.*d.*d.*d;

v = c0*exp((c1-r)/c2).*cos(c3*(c1-r)*pi) + c4*d6;

