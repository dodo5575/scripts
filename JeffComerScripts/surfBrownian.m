clear all;

% Parameters:
steps = 500000;
cutoff = 5;
surf = 4;
diffusion = 91; % in A^2/ns
dt = 0.0001; % ns
kT = 1; % kT (assuming input grids are in kT/A)
% Output:
displayPeriod = 1000;
trajDir = 'traj';
runName = 'size';

% Deprecated:
simArea = 2*25*25;
simNum = 200;
simRep = 16;

nameList = {'anneal' 'middling' 'raw1' 'phant'};
surfList = [12.91 13.25 8.574 16.66];
surfList1 = [67.87 68.16 63.34 71.57];
fracList = [0.96624 0.31193 0.18600 0.28];
topFracList = [0.46035 0.14502 0.06585 0.14];

name = nameList{surf};
surfZ = surfList(surf);
surfZ1 = surfList1(surf);
flowFrac = fracList(surf);

outFileText = sprintf('brownian_%s.txt', name);
outFile = sprintf('surf_%s.dat', name);
outFileFreeEnergy = sprintf('free_energy_%s.dat', name);
inFile = sprintf('shift_crop_pmf_%s_super.dx', name);
[data grid delta origin] = readdx(inFile);

%cutF = -1.5;
%data = (data<cutF).*(cutF*ones(size(data))) + (data>=cutF).*data;

% Add extra length to the bulk.
finalPosZ = 0.5*(surfZ + surfZ1);
lengthZ = finalPosZ - origin(3);
nodesZ = floor(lengthZ/delta(3,3));
gridNodes = length(data);
g = [grid(1); grid(2); nodesZ];
stateNodes = g(1)*g(2)*g(3);
extraNodes = stateNodes - gridNodes;
area = norm(cross(g(1)*delta(:,1),g(2)*delta(:,2)));
num = 1;
extraPos = origin + delta*(grid + [0; 0; 5]);
dim = delta*g;

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

% Open the output files.
outFile = sprintf('%s/traj_%s_%s.pdb', trajDir, runName, name);
outFileForce = sprintf('%s/fz_%s_%s.dat', trajDir, runName, name);
outForce = fopen(outFileForce, 'w');
outFileEvent = sprintf('%s/event_%s_%s.dat', trajDir, runName, name);
outEvent = fopen(outFileEvent, 'w');

% Make an initial configuration.
r = [0; 0; 24];
traj(1,:) = r';
trajCount = 1;

% Begin the Brownian Dynamics.
fprintf('Running %d Brownian dynamics steps.\n', steps);
eventTime = 0;
zc = surfZ + cutoff;
event = 0;
boundTime = 0;
boundSteps = 0;
r0 = r;
for s=1:steps
    % Get the time.
    t = s*dt;
    
    % Compute the force.
    if r(3) > extraPos(3)
        f = [0; 0; 0];
    else
        f = interpolateForce(pot, g, delta, origin, r);
    end
        
    % Perform a Brownian Dynamics step.
    rando = randn(3,1);
    r = r + diffusion*dt/kT*f + sqrt(2*diffusion*dt)*rando;
    
    % Wrap on x and y.
    if r(1) < origin(1)
        r(1) = r(1) + dim(1);
    end
    if r(1) >= origin(1) + dim(1)
        r(1) = r(1) - dim(1);
    end
    if r(2) < origin(2)
        r(2) = r(2) + dim(2);
    end
    if r(2) >= origin(2) + dim(2)
        r(2) = r(2) - dim(2);
    end
    
    % Reflect on z.
    if r(3) < origin(3)
        dz = r(3) - origin(3);
        r(3) = r(3) - 2*dz;
    end
    if r(3) >= finalPosZ;
        dz = r(3) - finalPosZ;
        r(3) = r(3) - 2*dz;
    end
    
    % Check for binding events.
    if event
        if r(3) > zc
            % A complete binding event.
            dur = t - eventTime;
            fprintf(outEvent, '%f %f\n', eventTime, dur);
            event = 0;
            boundTime = boundTime + dur;
        end
    else
        if r(3) < zc
            % New binding event.
            eventTime = t;
            event = 1;
        end
    end
    
    if mod(s,displayPeriod) == 0
        % Write the trajectory file.
        trajCount = trajCount + 1;
        traj(trajCount,:) = r0';
        
        fprintf('step: %d\n', s);
        fprintf(outForce, '%f %f\n', r0(3), f(3));
        
        %outTraj = sprintf('%s/traj_%s.%d.xyz', trajDir, name, s);
        %out = fopen(outTraj, 'w');
        %fprintf(out, '%d\n', num);
        %fprintf(out, 'DMMP system\n');
        %fprintf(out, 'P %f %f %f\n', r(1), r(2), r(3));
        %fclose(out);
    end
    r0 = r;
end
fclose(outForce);
fclose(outEvent);

timeFrac = boundTime/(dt*steps);
outFileFrac = sprintf('%s/frac_%s.txt', trajDir, name);
outFrac = fopen(outFileFrac, 'w');
fprintf(outFrac, '%f\n', timeFrac);
fprintf(outFrac, '%f\n', flowFrac);
fclose(outFrac);

fprintf('time fraction: %f\n', timeFrac);
fprintf('flow fraction: %f\n', flowFrac);

writePdbTraj(outFile, traj, dim);
