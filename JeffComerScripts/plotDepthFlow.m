clear all;

% Parameters:
membrane = 4;
nameList = {'anneal' 'middling' 'raw1' 'phant'};
surfZList = [12.91 13.25 8.574 16.66];
run = 'super';
rangeV = [-12.5 0];
rangeZ = [-0.25 0.15];
%rangeZ = [0.15 0.5];

lx = 2.5;
ly = 2.5;
shift = 0;
%tick = [-1.25 -1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 1.25];
tick = [-1.0 -0.5 0 0.5 1.0];
plotType = 'eps';

% Input:
name = nameList{membrane};
surfZ = surfZList(membrane);
gridFile = sprintf('shift_crop_pmf_%s_%s.dx', name, run);
depthFile = sprintf('shift_crop_pmf_%s_%s.dx.depth', name, run);
topoFile = sprintf('../topography/topo_high_probe2_%s.dat', name);
flowFile = sprintf('../rogan/flow_combine_%s.dat', name);
profileFile = sprintf('shift_crop_pmf_%s_%s.dx.depth.profile', name, run);
% Output:
histFile = sprintf('depth_hist_%s.dat', name);
medianFile = sprintf('depth_median_%s.dat', name);
maxFile = sprintf('depth_max_%s.dat', name);
minFile = sprintf('depth_min_%s.dat', name);
min1File = sprintf('depth_min1_%s.dat', name);
depthPlot = sprintf('good_depth_%s.%s', name, plotType);
topoPlot = sprintf('good_topo_%s.%s', name, plotType);
flowPlot = sprintf('good_topo_%s.%s', name, plotType);
topoData = sprintf('topo_map_%s.dat', name);

levV = linspace(rangeV(1),rangeV(2),14);
levZ = linspace(rangeZ(1),rangeZ(2),14);
mapV = makeColorMap([0 0 0; 0.15 0.15 1; 0.6 1 0.6; 0.9 1 0.9; 1 1 1], length(levV)+1);
mapZ = makeColorMap([0 0 0; 0.4 0.34 0.28; 0.8 0.68 0.55; 1 1 0.85; 1 1 1], length(levZ)+1);

% Read the grid attributes.
[grid delta origin] = readDxAttr(gridFile);
nx = grid(1);
ny = grid(2);
nz = grid(3);

% Read the data.
data = dlmread(depthFile, ' ');
x = data(:,1);
y = data(:,2);
z0 = data(:,3);
v0 = data(:,4);
n = length(v0);

% Load the profile data.
prof = dlmread(profileFile, ' ');
profZ = prof(1,:)';

% Load the flow data.
data = dlmread(flowFile, ' ');
f = data(:,3);
df = std(f);
dv = std(v0);
%f = (f - mean(f))*dv/df + mean(v0);
f = (f - mean(f)) + mean(v0);

% Match the topographical data.
data = dlmread(topoFile, ' ');
z = zeros(size(v0));
%dx = x(nx+1)-x(1);
%dy = y(2)-y(1);
for j=1:n
    k = 1;
    dx = data(k,1)-x(j);
    dy = data(k,2)-y(j);
    closeK = k;
    closeDist = dx*dx + dy*dy;
    
    for k=2:length(data(:,1))
        dx = data(k,1)-x(j);
        dy = data(k,2)-y(j);
        dist = dx*dx + dy*dy;
        
        if (dist < closeDist)
            closeK = k;
            closeDist = dist;
        end
    end
    
    %disp(sprintf('%f %f and %f %f', data(closeK,1), data(closeK,2), x(j), y(j)));
    z(j) = data(closeK,3);
    disp(closeK)
end
surfZ = mean(z);

% Convert to nanometers from angstroms.
z0 = 0.1*z0;
z = 0.1*z;
x = 0.1*x;
y = 0.1*y;
profZ = 0.1*profZ;
surfZ = 0.1*surfZ;

% Find the surface level.
z = z - surfZ;
z0 = z0 - surfZ;
profZ = profZ - surfZ;
dlmwrite(topoData, [x y z], ' ');

% Make the histogram.
[freq binV] = hist(v0, 12);
binSize = binV(2)-binV(1);
%freq = freq/(length(v0)*binSize);
dlmwrite(histFile, [binV' freq'], ' ');

% Sort the profiles.
[v0Sort ind] = sort(v0);
% Write the median, maximum, and minimum profiles.
indMedian = ind(floor(n/2));
indMin = ind(1);
indMin1 = ind(2);
indMax = ind(end);
dlmwrite(medianFile, [profZ prof(indMedian+1,:)'], ' ');
dlmwrite(min1File, [profZ prof(indMin1+1,:)'], ' ');
dlmwrite(minFile, [profZ prof(indMin+1,:)'], ' ');
dlmwrite(maxFile, [profZ prof(indMax+1,:)'], ' ');

% Make the cention.
secX = zeros(nx,ny);
secY = zeros(nx,ny);
secZ0 = zeros(nx,ny);
secV0 = zeros(nx,ny);depth_median_phant.dat
secZ = zeros(nx,ny);
secF = zeros(nx,ny);
for ix=0:nx-1
  for iy=0:ny-1
    j = 1 + iy + ix*ny;
    secX(ix+1, iy+1) = x(j);
    secY(ix+1, iy+1) = y(j);
    secZ0(ix+1, iy+1) = z0(j);
    secV0(ix+1, iy+1) = v0(j);
    secZ(ix+1, iy+1) = z(j);
    secF(ix+1, iy+1) = f(j);
  end
end

% Make the section.
cenX = zeros(nx,ny);
cenY = zeros(nx,ny);
cenZ0 = zeros(nx,ny);
cenV0 = zeros(nx,ny);
cenZ = zeros(nx,ny);
cenF = zeros(nx,ny);
for ix=1:nx
  for iy=1:ny
    jx = ix + shift;
    jy = iy + shift;
    dx = 0.0;
    dy = 0.0;
    if (jx > nx)
       jx = jx - nx;
       dx = lx;
    end
    if (jy > ny)
        jy = jy - ny;
        dy = ly;
    end
    
    cenX(ix,iy) = secX(jx,jy) + dx;
    cenY(ix,iy) = secY(jx,jy) + dy; 
    cenZ0(ix,iy) = secZ0(jx,jy);
    cenV0(ix,iy) = secV0(jx,jy); 
    cenZ(ix,iy) = secZ(jx,jy);
    cenF(ix,iy) = secF(jx,jy);
  end
end


figure(1)
%pcolor(cenX, cenY, cenV0)
contourf(cenX, cenY, cenV0, levV)
colormap(mapV)
xlabel('x (nm)', 'FontSize', 16)
ylabel('y (nm)', 'FontSize', 16)
%axis vis3d
%shading('interp')
axis square
caxis(rangeV)
colorbar
set(gcf, 'color', 'white'); % sets the color to white
set(gca, 'TickDir', 'in', 'XTick', tick, 'YTick', tick);
set(gca, 'FontSize', 16);
tickSize = 0.03*(cenX(end)-cenX(1));
tickSize1 = 0.035*(cenX(end)-cenX(1));
for j=1:length(tick)
    h = line([tick(j) tick(j)], [cenY(1) cenY(1)+tickSize1]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [cenY(1) cenY(1)+tickSize]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([tick(j) tick(j)], [cenY(end)-tickSize1 cenY(end)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [cenY(end)-tickSize cenY(end)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([cenX(1) cenX(1)+tickSize1], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([cenX(1) cenX(1)+tickSize], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([cenX(end)-tickSize1 cenX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([cenX(end)-tickSize cenX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
end
saveas(gcf, depthPlot, plotType);

figure(2)
%pcolor(cenX, cenY, cenZ)
contourf(cenX, cenY, cenZ, levZ)
colormap(mapZ)
xlabel('x (nm)', 'FontSize', 16)
ylabel('y (nm)', 'FontSize', 16)
axis square
caxis(rangeZ)
colorbar
set(gcf, 'color', 'white'); % sets the color to white
set(gca, 'TickDir', 'in', 'XTick', tick, 'YTick', tick);
set(gca, 'FontSize', 16);
tickSize = 0.03*(cenX(end)-cenX(1));
tickSize1 = 0.035*(cenX(end)-cenX(1));
for j=1:length(tick)
    h = line([tick(j) tick(j)], [cenY(1) cenY(1)+tickSize1]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [cenY(1) cenY(1)+tickSize]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([tick(j) tick(j)], [cenY(end)-tickSize1 cenY(end)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [cenY(end)-tickSize cenY(end)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([cenX(1) cenX(1)+tickSize1], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([cenX(1) cenX(1)+tickSize], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([cenX(end)-tickSize1 cenX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([cenX(end)-tickSize cenX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
end
saveas(gcf, topoPlot, plotType);

figure(3)
%pcolor(cenX, cenY, cenV0)
contourf(cenX, cenY, cenF, levV)
colormap(mapV)
xlabel('x (nm)', 'FontSize', 16)
ylabel('y (nm)', 'FontSize', 16)
%axis vis3d
%shading('interp')
axis square
caxis(rangeV)
colorbar
set(gcf, 'color', 'white'); % sets the color to white
set(gca, 'TickDir', 'in', 'XTick', tick, 'YTick', tick);
set(gca, 'FontSize', 16);
tickSize = 0.03*(cenX(end)-cenX(1));
tickSize1 = 0.035*(cenX(end)-cenX(1));
for j=1:length(tick)
    h = line([tick(j) tick(j)], [cenY(1) cenY(1)+tickSize1]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [cenY(1) cenY(1)+tickSize]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([tick(j) tick(j)], [cenY(end)-tickSize1 cenY(end)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([tick(j) tick(j)], [cenY(end)-tickSize cenY(end)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([cenX(1) cenX(1)+tickSize1], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([cenX(1) cenX(1)+tickSize], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
    
    h = line([cenX(end)-tickSize1 cenX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 3.5, 'Color', [1 1 1]);
    h = line([cenX(end)-tickSize cenX(end)], [tick(j) tick(j)]);
    set(h,'LineWidth', 2, 'Color', [0 0 0]);
end
saveas(gcf, flowPlot, plotType);

figure(4)
plot(binV, freq, 'k-')
xlabel('pmf (k_B T)')
%ylabel('probability/(k_B T)')
ylabel('counts')

