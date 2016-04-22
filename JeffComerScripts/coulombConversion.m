clear all;

k = 1.3806504e-23; % J/K
temp = 295; % K
e = 1.602176487e-19; % C
eps0 = 8.854187817e-12; % F/m

coulomb = 1e10/(k*temp)*e^2/(4*pi*eps0) % A kT/e^2
