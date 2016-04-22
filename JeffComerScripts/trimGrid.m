% Trim a grid.
clear all;

inFile = 'pore2.0_nostick2.dx';
outFile = 'pore2.0_nonstick2.dx';
%radius = 2/3.*66.857;
dx = 36.;
dy = 36.;

% Read the grid.
[data grid delta origin] = readdx(inFile);

lx0 = ceil((-dx-origin(1))/delta(1,1));
ly0 = ceil((-dy-origin(2))/delta(2,2));
lx1 = floor((dx-origin(1))/delta(1,1));
ly1 = floor((dy-origin(2))/delta(2,2));
lz0 = 1;
lz1 = grid(3);

gridNew = [lx1-lx0+1 ly1-ly0+1 lz1-lz0+1];
originNew = [-0.5*delta(1,1)*gridNew(1) -0.5*delta(2,2)*gridNew(2) ...
    -0.5*delta(3,3)*gridNew(3)];

n = prod(gridNew);
dataNew = zeros(n,1);
disp(sprintf('\nnew grid size: %d %d %d', gridNew(1), gridNew(2), gridNew(3)));

index = 1;
for ix=lx0:lx1
    for iy=ly0:ly1
        for iz=lz0:lz1
            l = iz + grid(3)*(iy-1) + grid(3)*grid(2)*(ix-1);
            dataNew(index) = data(l);
            
            index = index + 1;  
        end
    end
end

% Write the grid.
writedx(dataNew, gridNew, delta, originNew, outFile);


