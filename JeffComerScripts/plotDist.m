% Author: Jeff Comer <jcomer2@illinois.edu>
clear all;

fig = 1;
probMax = 3e-4;
%probMax = 0.6;
%name = 'dens_conform_ss24_str3';
%name = 'dens_conform_ss48_str4';
name = 'dens_conform_ss12_str5';
secS = dlmread(sprintf('%s_cyl.s',name), ' ');
secZ = dlmread(sprintf('%s_cyl.z',name), ' ');
secV = dlmread(sprintf('%s_cyl.pot',name), ' ');

secS = 0.1*secS;
secZ = 0.1*secZ;
%total = trapz(trapz(secS.*secV));
%secV = secV/total;

% Plot the contours.
figure(fig)
contourf(secS, secZ, secV, linspace(0,probMax,10))
colormap(jet)
xlabel('s (nm)')
ylabel('z (nm)')
axis([0 14 -14 14]);
%axis equal
colorbar
set(gcf, 'PaperOrientation', 'portrait');

if 0
figure(3)
contourf(secZ, secS, secV, linspace(0,probMax,10))
colormap(jet)
xlabel('z (nm)')
ylabel('s (nm)')
axis([-14 14 0 14]);
%axis equal
colorbar
set(gcf, 'PaperOrientation', 'landscape');
end



