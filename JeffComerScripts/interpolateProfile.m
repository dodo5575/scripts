clear all;

% Parameters:
surfList = {'flip' 'raw' 'middling' 'anneal'};
profileFile = 'cube_middling_hard.zavg.dat';
barrierZ = -8;

profile = dlmread(profileFile, ' ');
profileZ = profile(:,1);
profileV = profile(:,2);

for j=1:length(surfList)
    surf = surfList{j};
    inFile = sprintf('1d_%s_cen.dat', surf);
    
    data = dlmread(inFile, ' ');
    z = data(:,1);
    v = data(:,2);
    
    bz = barrierZ;
    z0 = z(1);
    a = v(1);
    b = (v(2)-v(1))/(z(2)-z(1));
    c = -(3*a/bz^2 + 2*b/bz);
    d = (2*a/bz^3 + b/bz^2);
    
    profileSurf = pchip(z, v, profileZ);
    n = length(profileZ);
    for k=1:n
        if profileZ(k) > z(end)
            profileSurf(k) = 0;
        elseif profileZ(k) < z(1)
            if profileZ(k) > barrierZ
                x = profileZ(k) - z0;
                profileSurf(k) = a + b*x + c*x^2 + d*x^3;
            else
                profileSurf(k) = 0;
            end
        end
    end
    
    outFile = sprintf('1d_%s_profile.dat', surf);
    dlmwrite(outFile, [profileZ profileSurf], ' ');
end

