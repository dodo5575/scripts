clear all;

% Parameters:
steps = 200000;
num = 100;
diameter = 2.6; % in A
cutoff = 5;
surf = 1;
displayPeriod = 5000;
cutStep = floor(steps/10);

nameList = {'anneal' 'middling' 'raw1' 'phant'};
surfList = [12.91 13.25 8.574 16.66];
surfList1 = [67.87 68.16 63.34 71.57];
fracList = [0.96624 0.31193 0.18600 0.28];
topFracList = [0.46035 0.14502 0.06585 0.14];

name = nameList{surf};
surfZ = surfList(surf);
surfZ1 = surfList1(surf);

inFile = sprintf('shift_crop_pmf_%s_super.dx', name);
[data grid delta origin] = readdx(inFile);

% Add extra length to the bulk.
extraLength = 0.5*(surfZ + surfZ1) - origin(3);
extraNodes = grid(1)*grid(2)*floor(extraLength/delta(3,3));
gridNodes = length(data);
stateNodes = gridNodes + extraNodes;

% Find the surface nodes.
isSurf = zeros(stateNodes,1);
surfNodes = 0;
surfW = 0;
bulkW = 0;
for ix=0:grid(1)-1
    for iy=0:grid(2)-1
        for iz=0:grid(3)-1
            r = delta*[ix; iy; iz] + origin;
            if r(3) < surfZ + cutoff
                surfNodes = surfNodes + 1;
                node = 1 + iz + iy*grid(3) + ix*grid(3)*grid(2);
                surfW = surfW + exp(-data(node));
                isSurf(node) = 1;
            else
                bulkW = bulkW + exp(-data(node));
            end
        end
    end
end
bulkW = bulkW + extraNodes;
indepFrac = surfW/(bulkW + surfW);

% Initially put all particles in the extra space.
partIndex = gridNodes*ones(num,1);
surfTraj = zeros(steps,1);
surfCount = 0;
acceptCount = 0;
for s=1:steps
    % Pick a particle.
    part = ceil(rand*num);
    index0 = partIndex(part);
    
    % Get the initial free energy.
    if index0 >= gridNodes
        % in the extra space
        energy0 = 0;
    else
        % on the grid
        energy0 = data(index0+1);
    end
    
    % Choose a trial index.
    index1 = floor(rand*stateNodes);
    r1 = delta*getIndices(index1, grid) + origin;
    
    % Get the trial free energy.
    collide = 0;
    if index1 >= gridNodes
        % We picked a node in the extra space.
        energy1 = 0;
    else
        % We picked a node in the grid.
        energy1 = data(index1+1);
        
        % Check that this particle doesn't collide with any other
        % particle.
        for p=1:num
            if p==part, continue, end
            if partIndex(p) >= gridNodes, continue, end
            
            r = delta*getIndices(partIndex(p), grid) + origin;
            if norm(r-r1) < diameter
                collide = 1;
                break
            end
        end
    end

    % A collision is infinite energy--such moves are not accepted.
    % Accept other moves according to the metropolis criterion.
    if collide == 0 && (energy1 < energy0 || exp(energy0-energy1) > rand)
        % Accept the move.
        partIndex(part) = index1;
        surfCount = surfCount + isSurf(index1+1) - isSurf(index0+1);
        acceptCount = acceptCount + 1;
    end
  
    if mod(s,displayPeriod) == 0
        disp(sprintf('step %d: %f', s, surfCount));
    end
    surfTraj(s) = surfCount;
end

acceptRatio = acceptCount/steps
meanSurf = mean(surfTraj(cutStep:end));
meanFrac = meanSurf/num
indepFrac
plot(surfTraj)
