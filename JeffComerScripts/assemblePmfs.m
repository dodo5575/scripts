clear all;


% Input:
run0 = 'reprise';
run1 = 'ions';
type = [1 2];
typeName = {'pot' 'chl'};
dir = '.';
% Output
outName = 'new_code';
% Parameters:
charge = [1 -1];
radius = [1.76375 2.27];
eps = [0.0870 0.15];
coulombConst = 566.443209/92;
levelZ0 = 10;
levelZ1 = 12;
longDist = 30;

% Begin code.
s1 = type(1);
s2 = type(2);
name = sprintf('%s-%s', typeName{s1}, typeName{s2});

inName0 = sprintf('%s/pmf_%s_%s_%s.dx.z.dat', dir, 'radial', name, run0);
data0 = dlmread(inName0, ' ');
r0 = data0(:,1);
v0 = data0(:,2);
[mv0 ev0] = meanOnInterval(r0, levelZ0, levelZ1, v0);
dx = r0(2)-r0(1);

inName1 = sprintf('%s/pmf_%s_%s_%s.dx.z.dat', dir, 'radial', name, run1);
data1 = dlmread(inName1, ' ');
r1 = data1(:,1);
v1 = data1(:,2);
[mv1 ev1] = meanOnInterval(r1, levelZ0, levelZ1, v1);

r = (r0(1):dx:longDist)';
u = coulombConst*charge(s1)*charge(s2)./r + computeLennardJones(r, eps(s1), eps(s2), radius(s1), radius(s2));
uc = coulombConst*charge(s1)*charge(s2)./r;
mu = meanOnInterval(r, levelZ0, levelZ1, u);

% Fit the PMFs to one another.
v0 = v0 - mv0 + mu;
v1 = v1 - mv1 + mu;
% Average the pmfs.

% Merge them together.
nr = floor(r0(1)/dx)-1;
nearR0 = r0(1)-nr*dx;
nearR1 = r0(1)-dx;
nearR = (nearR0:dx:nearR1)';
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
tableFR = tableR(2:end)-0.5*dx;
tableF = -diff(tableV)/dx;

% Plot stuff.
figure(1)
plot(r0, v0, 'r--', r1, v1, 'b--', r, u, 'k--', tableR, tableV, 'go')
axis([0 inf -3 6])

% Write the results.
outForceName = sprintf('%s/force_%s.dat', dir, name);
dlmwrite(outForceName, [tableFR tableF], ' ');

outName = sprintf('%s/%s_%s.dat', dir, outName, name);
dlmwrite(outName, [tableR tableV], ' ');
fprintf('Wrote %s.\n', outName);
