clear all;
inGrid = 'diffGrid.dx';

outPrefix = 'diffusion_map/diffusion';
%outName = 'no_pot';
%diffuse = [237 237 181];
%outName = 'no_chl';
%diffuse = [244 244 187];
%outName = 'no_1M_pot';
%diffuse = [201 201 162];
outName = 'no_1M_chl';
diffuse = [197 197 157];

%outPrefix = 'diffusion_map/diffDiff';
%outName = 'no_1M_pot';
%diffuse = [201 189 162];
%outName = 'no_1M_chl';
%diffuse = [197 194 157];
%outName = 'no_pot';
%diffuse = [237 234 181];
%outName = 'no_chl';
%diffuse = [244 237 187];

% Region parameters.
sIn = 3;
zIn = 15;
zTop = 21;
sTop = 11;
zIntercept = 13;

outFile = sprintf('%s_%s.dx', outPrefix, outName);
[data grid delta origin] = readdx(inGrid);

nz = grid(3);
ny = grid(2);
nx = grid(1);
nynz = ny*nz;
for ix=0:nx-1
    for iy=0:ny-1
        for iz=0:nz-1
            r = delta*[ix; iy; iz] + origin;
            
            z = abs(r(3));
            s = sqrt(r(1)*r(1) + r(2)*r(2));
                        
            % Find which region we are in.
            where = 0;
            if z <= zIn && s <= sIn
                where = 2;
            elseif z <= zIn && s > sIn
                where = 3;
            elseif z <= zTop && s > sTop
                where = 3;
            elseif z > zTop
                where = 1;
            elseif z > zTop && s <= sIn
                where = 1;
            end
            
            % Handle the region for the corner of the pore.
            if where == 0
                if z > 1.0*s + 13
                    where = 1;
                else
                    where = 3;
                end
            end
            
            j = 1 + iz + iy*nz + ix*nynz;
            data(j) = diffuse(where);
        end
    end
end

writedx(data, grid, delta, origin, outFile);

