%resistivity of 1MKCl and Mg (m/S)
resistivity=0.060796827;

%box size in angstrom
boxX = 106.1;
boxY = 50.5987;
boxZ = 92.481;


%read total current data (nA)
V100_10ns = load('square2plate-1MKCl-100mV/curr_square2plate_1MKCl_100mV.dat');
V100_50ns = load('square2plate-1MKCl-100mV/curr_square2plate_1MKCl_100mV_50ns.dat');
V250_10ns = load('square2plate-1MKCl-250mV/curr_square2plate_1MKCl_250mV.dat');
V250_50ns = load('square2plate-1MKCl-250mV/curr_square2plate_1MKCl_250mV_50ns.dat');
V500_10ns = load('square2plate-1MKCl-500mV/curr_square2plate_1MKCl_500mV.dat');
V500_50ns = load('square2plate-1MKCl-500mV/curr_square2plate_1MKCl_500mV_50ns.dat');

%combine current data (nA)
V100_60ns = [V100_10ns;V100_50ns];
V100_60ns = [mean(V100_60ns);V100_60ns;mean(V100_60ns)];
V250_60ns = [V250_10ns;V250_50ns];
V250_60ns = [mean(V250_60ns);V250_60ns;mean(V250_60ns)];
V500_60ns = [V500_10ns;V500_50ns];
V500_60ns = [mean(V500_60ns);V500_60ns;mean(V500_60ns)];

%calculate totatl resistance (ohm)
V100_resistance = 0.1 ./ (V100_60ns(:,2) * 10^-9);
V250_resistance = 0.25 ./ (V250_60ns(:,2) * 10^-9);
V500_resistance = 0.5 ./ (V500_60ns(:,2) * 10^-9);


%read origami thickness data (angstrom)
V100_DNA_z_10ns = load('square2plate-1MKCl-100mV/square2plate-1MKCl-100mV-OrigamiMinMax.dat');
V100_DNA_z_50ns = load('square2plate-1MKCl-100mV/square2plate-1MKCl-100mV-50ns-OrigamiMinMax.dat');
V250_DNA_z_10ns = load('square2plate-1MKCl-250mV/square2plate-1MKCl-250mV-OrigamiMinMax.dat');
V250_DNA_z_50ns = load('square2plate-1MKCl-250mV/square2plate-1MKCl-250mV-50ns-OrigamiMinMax.dat');
V500_DNA_z_10ns = load('square2plate-1MKCl-500mV/square2plate-1MKCl-500mV-OrigamiMinMax.dat');
V500_DNA_z_50ns = load('square2plate-1MKCl-500mV/square2plate-1MKCl-500mV-50ns-OrigamiMinMax.dat');

%combine origami thickness data (angstrom)
V100_DNA_z_60ns = [V100_DNA_z_10ns;V100_DNA_z_50ns];
V250_DNA_z_60ns = [V250_DNA_z_10ns;V250_DNA_z_50ns];
V500_DNA_z_60ns = [V500_DNA_z_10ns;V500_DNA_z_50ns];


%change the unit of origami thickness data (m)
V100_DNA_z = V100_DNA_z_60ns(:,2) * 10^-10;
V250_DNA_z = V250_DNA_z_60ns(:,2) * 10^-10;
V500_DNA_z = V500_DNA_z_60ns(:,2) * 10^-10;



%calculate the length of solution in z (m)
V100_Sol_z = (boxZ - V100_DNA_z_60ns(:,2)) * 10^-10;
V250_Sol_z = (boxZ - V250_DNA_z_60ns(:,2)) * 10^-10;
V500_Sol_z = (boxZ - V500_DNA_z_60ns(:,2)) * 10^-10;


%calculate the resistance of solution (ohm)
V100_Sol_resistance = (resistivity * V100_Sol_z) / (boxX * 10^-10 * boxY * 10^-10);
V250_Sol_resistance = (resistivity * V250_Sol_z) / (boxX * 10^-10 * boxY * 10^-10);
V500_Sol_resistance = (resistivity * V500_Sol_z) / (boxX * 10^-10 * boxY * 10^-10);


%calculate the resistance of origami (ohm)
V100_DNA_resistance = V100_resistance - V100_Sol_resistance;
V250_DNA_resistance = V250_resistance - V250_Sol_resistance;
V500_DNA_resistance = V500_resistance - V500_Sol_resistance;


%calculate the resistivity of origami (m/S)
V100_DNA_resistivity = (V100_DNA_resistance * (boxX * 10^-10 * boxY * 10^-10)) ./ V100_DNA_z; 
V250_DNA_resistivity = (V250_DNA_resistance * (boxX * 10^-10 * boxY * 10^-10)) ./ V250_DNA_z;
V500_DNA_resistivity = (V500_DNA_resistance * (boxX * 10^-10 * boxY * 10^-10)) ./ V500_DNA_z;


%calculate the conductivity of origami (S/m)
V100_DNA_conductivity = 1 ./ V100_DNA_resistivity;
V250_DNA_conductivity = 1 ./ V250_DNA_resistivity;
V500_DNA_conductivity = 1 ./ V500_DNA_resistivity;

mean(V100_DNA_conductivity)
std(V100_DNA_conductivity)
mean(V250_DNA_conductivity)
std(V250_DNA_conductivity)
mean(V500_DNA_conductivity)
std(V500_DNA_conductivity)

