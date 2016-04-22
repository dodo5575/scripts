% Gaussian random data.
% Author: Jeff Comer <jcomer2@illinois.edu>
clear all;

inName = 'matlab_window_fake.txt';
outPrefix = 'log/';
count = 10000;

data = dlmread(inName, ' ');
kBT = 0.5862292; % in kcal/mol at 295 K

ind = data(:,1);
x = data(:,2);
y = data(:,3);
z = data(:,4);
kx = data(:,5)/kBT;
ky = data(:,6)/kBT;
kz = data(:,7)/kBT;

t = 100*((1:count)'-1);
for j=1:length(x)
    % Find the standard deviations.
    dx = sqrt(1/kx(j));
    dy = sqrt(1/ky(j));
    dz = sqrt(1/kz(j));
    
    rx = x(j) + dx*randn(count,1);
    ry = y(j) + dy*randn(count,1);
    rz = z(j) + dz*randn(count,1);
    
    m = [t rx ry rz];
    dlmwrite(sprintf('%s%d.log.pos', outPrefix, ind(j)), m, ' ');
end
