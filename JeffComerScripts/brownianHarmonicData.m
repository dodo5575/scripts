clear all;

dt = 16e-6;
diffuse = 227; % in A^2/ns
spring = 4.0; % in kcal/mol/A^2
kT = 0.5862292; % in kcal/mol at 295 K
subSteps = 600;
subTrajNum = 10000;
refPos = 20;
pos0 = 20.5;
outFile = 'data/brownian.dat';

% Do Brownian Dynamics.
pos = zeros(subTrajNum*subSteps,1);
c1 = diffuse*dt/kT;
c2 = sqrt(2*diffuse*dt);

pos(1) = pos0;
n = subSteps*subTrajNum;
for j=2:n
    force = -spring*(pos(j-1) - refPos);
    pos(j) = pos(j-1) + c1*force + c2*randn;
end

subTraj = reshape(pos, subSteps, subTrajNum)';

dlmwrite(outFile, subTraj, ' ');
fprintf('Wrote %d subtrajectories.\n', subTrajNum);