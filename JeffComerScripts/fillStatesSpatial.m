clear all;

% Parameters:
steps = 50000;
num = 100;
diameter = 2; % in A
cutoff = 5;
surf = 2;
displayPeriod = 1000;

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
extraNodes = grid(1)*grid(2)*(extraLength/delta(3,3));
gridNodes = length(data);
stateNodes = gridNodes + extraNodes;

surfNodes = 0;
surfW = 0;
bulkW = 0;
for ix=0:grid(1)-1
    for iy=0:grid(2)-1
        for iz=0:grid(3)-1
            r = delta*[0; 0; iz] + origin;
            if r(3) < surfZ + cutoff
                surfNodes = surfNodes + 1;
                node = 1 + iz + iy*grid(3) + ix*grid(3)*grid(2);
                surfW = surfW + exp(-data(node));
            else
                bulkW = bulkW + exp(-data(node));
            end
        end
    end
end
bulkW = bulkW + extraNodes;
fracIndep = surfW/(bulkW + surfW)

gridWeight = exp(-data);
weightSum = 0;
probSum = 0;
for s=1:steps
    % Clear the excluded spots.
    excludedPos = zeros(3,num);
    excluded = 0;
    surfCount = 0;
    freeEnergy = 0;
    
    for p=1:num
        picking = 1;
        while picking
            % Pick a state to fill.
            ind = floor(rand*stateNodes);
            
            if ind >= gridNodes
                % The free energy for the extra nodes is zero.
                % Extra nodes also exclude no space.
                picking = 0;
            else
                % Get the position of this point.
                l = getIndices(ind, grid);
                r = delta*l + origin;
                
                % Check that no other particles exclude this state.
                occupied = 0;
                if 0
                for k=1:excluded
                    d = norm(r-excludedPos(:,k));
                    if d < diameter
                        occupied = 1;
                        break
                    end
                end
                end
                
                if ~occupied
                    % We have picked a valid state.
                    picking = 0;
                    excluded = excluded + 1;
                    excludedPos(:,excluded) = r;
                                        
                    % Check if it is a surface state.
                    if r(3) < surfZ + cutoff
                        freeEnergy = freeEnergy + data(ind+1);
                        surfCount = surfCount + 1;
                    end
                end
            end
        end % picking loop    
    end % a single step
                
    % We have a complete state.
    % Add it to statistical sums.
    weight = exp(-freeEnergy);
    weightSum = weightSum + weight;
    probSum = probSum + weight*surfCount/num;
    
    if mod(s,displayPeriod) == 0
        disp(sprintf('step %d: %f', s, probSum/weightSum));
    end
    probTrack(s) = probSum;
end

%surf = num*probSum/weightSum
frac = probSum/weightSum

plot(probTrack,'-')
%hist3([trackX' trackY'])