clear all;


rangeV = [-1.2 0];
levV = linspace(rangeV(1),rangeV(2),12);

name = 'atom_empty';
secS = dlmread(strcat(name,'.y.sec'), ' ');
secZ = dlmread(strcat(name,'.z.sec'), ' ');
secV = dlmread(strcat(name,'.v.sec'), ' ');

mapV = makeColorMap([0 0 0; 0.15 0.15 1; 0.6 1 0.6; 0.9 1 0.9; 1 1 1], length(levV)+1);

% Plot things.
figure(1)
%contourf(secS, secZ, secEz, [0 1.2.^(-8:0)] )
%contourf(secS, secZ, secEz, logspace(-2,0,20))
[c h] = contourf(secS, secZ, secV, levV);
%[c h] = contourf(secS, secZ, secV, linspace(0,1.0,10))
colormap(mapV)
xlabel('y (nm)')
ylabel('z (nm)')
axis square
caxis(rangeV)
colorbar
%set(h, 'LineStyle', 'None');
%set(h, 'LineWidth', 2);
