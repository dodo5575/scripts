clear all;

name = {'coil_p1.2_6.5Va0-1' 'coil_p1.2_4Va0-5' 'coilx_p1.6_6.5Va0-2'...
    'coilx_p1.6_6.5V0-3' 'coilx_p1.6_4V0-14' 'coil_p2.0_2V0-7' 'coil_p2.0_4V0-2'};
constrict = {1.32 1.32 1.60 1.60 1.60 2.20 2.20};

angle = 10.0;
phi = pi/180.0*angle;

for j=1:length(name)
    data = dlmread(sprintf('last1_pos_%s.dat', name{j}), ' ');
    t = data(:,1);
    z = data(:,2);

    d = constrict{j} + 2.0*tan(phi)*abs(z);
    fileName = sprintf('diam_%s.dat', name{j});
    dlmwrite(fileName, [t d], ' ');
    disp(sprintf('Created diameter file "%s" with constriction %f.',...
        fileName, constrict{j}));
end



