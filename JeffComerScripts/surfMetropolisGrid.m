clear all;

% Parameters:
steps = 200000;
simArea = 2*25*25;
simNum = 200;
simRep = 16;
%diameter = 2.6; % in A
diameter = 9;
cutoff = 5;
surf = 2;
displayPeriod = 2000;
cutStep = floor(0.3*steps);
trajDir = 'traj';
outDir = 'output';

nameList = {'anneal' 'middling' 'raw1' 'phant' 'test'};
surfList = [12.91 13.25 8.574 16.66 13.25];
surfList1 = [67.87 68.16 63.34 71.57 68.16];
fracList = [0.96624 0.31193 0.18600 0.28 0.31193];
topFracList = [0.46035 0.14502 0.06585 0.14 0.14502];

name = nameList{surf};
surfZ = surfList(surf);
surfZ1 = surfList1(surf);

outFileText = sprintf('%s/monte_%s_d%d.txt', outDir, name, diameter);
outFile = sprintf('%s/surf_%s_d%d.dat', outDir, name, diameter);
outFileFreeEnergy = sprintf('%s/free_energy_%s_d%d.dat',outDir, name, diameter);
inFile = sprintf('shift_crop_pmf_%s_super.dx', name);
[data grid delta origin] = readdx(inFile);

%cutF = -1.5;
%data = (data<cutF).*(cutF*ones(size(data))) + (data>=cutF).*data;

% Add extra length to the bulk.
extraLength = 0.5*(surfZ + surfZ1) - origin(3);
extraNodesZ = floor(extraLength/delta(3,3));
gridNodes = length(data);
g = [grid(1) grid(2) extraNodesZ];
stateNodes = g(1)*g(2)*g(3);
extraNodes = stateNodes - gridNodes;
area = norm(cross(g(1)*delta(:,1),g(2)*delta(:,2)));
num = ceil(simNum*area/simArea/simRep);
disp(sprintf('Using %d particles.', num));

% Copy the pmf.
% Find the surface nodes.
pot = zeros(stateNodes,1);
isSurf = zeros(stateNodes,1);
surfNodes = 0;
surfW = 0;
bulkW = 0;
for ix=0:g(1)-1
    for iy=0:g(2)-1
        for iz=0:grid(3)-1
            node = 1 + iz + iy*g(3) + ix*g(3)*g(2);
            node0 = 1 + iz + iy*grid(3) + ix*grid(3)*grid(2);
            pot(node) = data(node0);
            
            r = delta*[ix; iy; iz] + origin;
            if r(3) < surfZ + cutoff
                surfNodes = surfNodes + 1;
                surfW = surfW + exp(-pot(node));
                isSurf(node) = 1;
            else
                bulkW = bulkW + exp(-pot(node));
            end
        end
    end
end
bulkW = bulkW + extraNodes;
indepFrac = surfW/(bulkW + surfW)
indepFreeEnergy = sum(exp(-pot).*pot)/sum(exp(-pot))

% Make an initial configuration.
partIndex = zeros(num,1);
surfCount0 = 0;
freeEnergy0 = 0;
for part=1:num
    % Check for collisions.
    collide = 1;
    while collide
        partIndex(part) = floor(rand*stateNodes);
        r1 = delta*getIndices(partIndex(part), g) + origin;
        
        collide = 0;
        for p=1:part-1
            r = delta*getIndices(partIndex(p), g) + origin;
            if norm(r-r1) < diameter
                collide = 1;
                break
            end
        end
    end
    
    freeEnergy0 = freeEnergy0 + pot(partIndex(part)+1);
    surfCount0 = surfCount0 + isSurf(partIndex(part)+1);
end

% Begin the Metropolis Monte Carlo.
disp(sprintf('Running %d Metropolis Monte Carlo steps.', steps));
surfTraj = zeros(steps,1);
freeEnergyTraj = zeros(steps,1);
acceptCount = 0;
for s=2:steps
    % Pick a particle.
    %part = ceil(rand*num);
    for part=1:num
        index0 = partIndex(part);
        
        % Get the initial free energy.
        energy0 = pot(index0+1);
        
        % Choose a trial index.
        index1 = floor(rand*stateNodes);
        r1 = delta*getIndices(index1, g) + origin;
        
        % Get the trial free energy.
        energy1 = pot(index1+1);
        
        collide = 0;
        % Check that this particle doesn't collide with any other
        % particle.
        for p=1:num
            if p==part, continue, end
            
            r = delta*getIndices(partIndex(p), g) + origin;
            if norm(r-r1) < diameter
                collide = 1;
                break
            end
        end
        
        % A collision is infinite energy--such moves are not accepted.
        % Accept other moves according to the metropolis criterion.
        if collide == 0 && (energy1 < energy0 || exp(energy0-energy1) > rand)
            % Accept the move.
            partIndex(part) = index1;
            %surfCount = surfCount + isSurf(index1+1) - isSurf(index0+1)
            acceptCount = acceptCount + 1;
        end
    end
    
    % Compute state values.
    surfTraj(s) = sum(isSurf(partIndex+1));
    freeEnergyTraj(s) = sum(pot(partIndex+1));
    
    if mod(s,displayPeriod) == 0
        if s > cutStep
            meanSurf = mean(surfTraj(cutStep:s))/num;
            meanFreeEnergy = mean(freeEnergyTraj(cutStep:s))/num;
        else
            meanSurf = surfTraj(s)/num;
            meanFreeEnergy = freeEnergyTraj(s)/num;
        end
        disp(sprintf('step %d: %f %f', s, meanSurf, meanFreeEnergy));
        
        % Write the trajectory file.
        outTraj = sprintf('%s/traj_%s_d%d.%d.xyz', trajDir, name, diameter, s);
        out = fopen(outTraj, 'w');
        fprintf(out, '%d\n', num);
        fprintf(out, 'DMMP system\n');
        for p=1:num
            r = delta*getIndices(partIndex(p), g) + origin;
            fprintf(out, 'P %f %f %f\n', r(1), r(2), r(3));
        end
        fclose(out);
        
        pause(0.05);
    end
end

acceptRatio = acceptCount/steps/num
meanSurf = mean(surfTraj(cutStep:end));
meanFrac = meanSurf/num
meanFreeEnergy = mean(freeEnergyTraj(cutStep:s))/num
indepFrac
indepFreeEnergy

out = fopen(outFileText, 'w');
fprintf(out, 'name: %s\n', name);
fprintf(out, 'particles %d\n', num);
fprintf(out, 'independent fraction: %f\n', indepFrac);
fprintf(out, 'mean fraction: %f\n', meanFrac);
fprintf(out, 'steps %d\n', steps);
fprintf(out, 'acceptance ratio %f\n', acceptRatio);
fprintf(out, 'mean free energy per particle %f\n', meanFreeEnergy);
fprintf(out, 'independent free energy per particle %f\n', indepFreeEnergy);
fclose(out);

dlmwrite(outFile, surfTraj/num, ' ');
dlmwrite(outFileFreeEnergy, freeEnergyTraj/num, ' ');
plot(freeEnergyTraj/num)
