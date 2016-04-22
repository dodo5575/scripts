function [data grid delta origin] = readdx(fileName)
% Read a .dx file.
% [DATA GRID DELTA ORIGIN] = READDX(FILENAME) returns the potential
% value at each point in DATA, the number of grid points along each
% cell axis in GRID (3-vector), the lattice vectors in DELTA
% (3x3 matrix), and the cell origin in ORIGIN (3-vector).
% Data order is z fast, y medium, and x slow.
% Author: Jeff Comer <jcomer2@illinois.edu>

in = fopen(fileName, 'r');
if in == -1
    error('readdx:fileNotFound', 'The .dx file was not found.')
end
disp(sprintf('Reading file %s.', fileName))

% Read the file contents.
n = 0;
items = 0;
data = [];
grid = zeros(1,3);
delta = zeros(3,3);
origin = zeros(1,3);
nDelta = 1;
while 1
    lin = fgets(in);
    
    if ~ischar(lin), break, end
    if isempty(lin)
        msg = ferror(in);
        error('readdx:ioError', msg)
    end

    % Cut off the newline and check for an empty line.
    line = strtrim(lin);
    if length(line) == 0, continue, end
        
    % Ignore comments.
    if line(1) == '#', continue, end
    
    % Read a line of numeric data.
    if n ~= items
        d = str2num(line);
                
        for j=1:length(d)
            n = n + 1;
            data(n) = d(j);
        end

    elseif all(line(1:6) == 'object')
        remain = line;
        while ~isempty(remain)
            [str, remain] = strtok(remain, ' ');
            if strcmp(str, 'class'), break, end
        end
        
        [str, remain] = strtok(remain, ' ');
        if strcmp(str, 'gridpositions')
            % Get the grid size.
            while ~isempty(remain)
                [str, remain] = strtok(remain, ' ');
                if strcmp(str, 'counts'), break, end
            end
            grid = str2num(strtrim(remain));
            disp(sprintf('grid size: %s', num2str(grid)))
        
        elseif strcmp(str, 'array')
            % Initialize the data array.
            while ~isempty(remain)
                [str, remain] = strtok(remain, ' ');
                if strcmp(str, 'items'), break, end
                n = 0;
            end
            [str, remain] = strtok(remain, ' ');
            items = str2double(str);
            data = zeros(items,1);
        end
    
    elseif all(line(1:6) == 'origin')
        origin = str2num(strtrim(line(7:end)));
    
    elseif all(line(1:5) == 'delta')
        delta(nDelta,:) = str2num(strtrim(line(6:end)));
        nDelta = nDelta + 1;
    else
    end
end
origin = origin';
grid = grid';
delta = delta';
disp(sprintf('Successfully read %i values.', n))

fclose(in);
    





