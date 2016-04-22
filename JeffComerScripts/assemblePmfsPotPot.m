clear all;

% Parameters:
run0 = 'reprise';
run1 = 'three';
name = 'pot-pot';
type = [1 1];
coulombConst = 566.443209/92;
charge = [1 -1];
radius = [1.76375 2.27];
eps = [0.0870 0.15];
dir = '.';
levelZ0 = 12;
levelZ1 = 15;
levelX0 = 4;
levelX1 = 6;
longDist = 30;

s1 = type(1);
s2 = type(2);

inName0 = sprintf('%s/pmf_%s_%s.dx.z.dat', dir, run0, 'chl-chl');
data0 = dlmread(inName0, ' ');
r0 = data0(:,1);
v0 = data0(:,2);
mv0 = meanOnInterval(r0, levelX0, levelX1, v0);
dx = r0(2)-r0(1);

inName1 = sprintf('%s/pmf_%s_%s.dx.z.dat', dir, run1, name);
data1 = dlmread(inName1, ' ');
r1 = data1(:,1);
v1 = data1(:,2);
mv1 = meanOnInterval(r1, levelX0, levelX1, v1);

% Fit the PMFs to one another.
v0 = v0 - mv0;
v1 = v1 - mv1;

% Shift horizontally.
shift0 = meanOnInterval(v0, 3, 6, r0);
shift1 = meanOnInterval(v1, 3, 6, r1);
r0 = r0 - shift0 + shift1;

r = (r0(1):dx:longDist)';
u = coulombConst*charge(s1)*charge(s2)./r + computeHardcore(r, eps(s1), eps(s2), radius(s1), radius(s2));
uc = coulombConst*charge(s1)*charge(s2)./r;
mu = meanOnInterval(r, levelX0, levelX1, u);

% Shift veritcally.
mv0 = meanOnInterval(r0, levelZ0, levelZ1, v0);
mv1 = meanOnInterval(r1, levelZ0, levelZ1, v1);
mu = meanOnInterval(r, levelZ0, levelZ1, u);
v0 = v0 - mv0 + mu;
v1 = v1 - mv1 + mu;

% Merge them together.
nearR = (0:dx:(r0(1)-dx))';
nearV = v0(1)*ones(size(nearR));
slope = (v0(3)-v0(1))/(2*dx);
for j=1:length(nearR)
    nearV(j) = slope*(nearR(j)-r0(1)) + v0(1);
end

farR = r;
farV = zeros(size(r));
for j=1:length(farV)
    if r(j) < levelZ1
        farV(j) = v0(j);
    else
        farV(j) = u(j);
    end
end

% Make the tabulated potentials.
tableR = [nearR; farR];
tableV = [nearV; farV];

% Plot stuff.
figure(1)
plot(r0, v0, 'r--', r1, v1, 'b--', r, u, 'k--', tableR, tableV, 'go')

% Write the results.
outName = sprintf('%s/helmholtz_%s.dat', dir, name);
dlmwrite(outName, [tableR tableV], ' ');
fprintf('Wrote %s.\n', outName);
