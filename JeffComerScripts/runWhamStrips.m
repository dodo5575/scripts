clear all;

% Parameters:
name = 'anneal';
smdFreq = 100;
cutTime0 = 200000;
kt = 0.5862292; % in kcal/mol at 295 K
levelDist = 1.5;
inPrefix = '../';
indexFile = sprintf('window_%s.txt', name);
outName = sprintf('%s/strip_%s', name, name);
makePmf = 1;

% Begin
cutIndex = round(cutTime0/smdFreq);

% Load the index file.
in = fopen(indexFile, 'r');
nameList = {};
dataList = [];
indexN = 0;
while 1
  lin = fgets(in);
  if ~ischar(lin), break, end
  if isempty(lin)
        msg = ferror(in);
        error('ioError %s', msg)
  end

  % Cut off the newline and check for an empty line.
  line = strtrim(lin);
  if isempty(line), continue, end
 
  indexN = indexN + 1;
  [name,data] = strtok(line);
  nameList{indexN} = name;
  dataList(indexN,:) = str2num(data);
end
fclose(in);

% Find the unique strip xy positions.
stripPos(1,:) = dataList(1,1:2);
stripCount(1) = 1;
stripIndex(1,1) = 1;
stripN = 1;
for k=2:indexN
    xy = dataList(k,1:2);
    
    % Search to see if this already exists.
    found = 0;
    for s=1:length(stripPos(:,1))
        if all(stripPos(s,:) == xy)
            found = 1;
            break;
        end
    end
    
    if ~found
        stripN = stripN + 1;
        stripPos(stripN,:) = xy;1
        stripCount(stripN) = 1;
        stripIndex(stripN,stripCount(stripN)) = k;
    else
        stripCount(s) = stripCount(s) + 1;
        stripIndex(s,stripCount(s)) = k;
    end
    
end
stripCount = stripCount';

% Compute the WHAM over all the strip xy positions.
for s=3:stripN
    % Skip lone data sets.
    if stripCount(s) < 2, continue, end
    winN = stripCount(s);
    
    % Pull out the window center and spring constant.
    center = dataList(stripIndex(s,1:winN),3);
    spring = dataList(stripIndex(s,1:winN),6);
    % Convert everything to kT.
    spring = spring/kt;
    
    mat = regexp(nameList{stripIndex(s,1)}, '(pos)\d+', 'match');
    stripName = mat{1}(4:end);
    disp(sprintf('\nComputing WHAM for strip %s at %f %f.', stripName, stripPos(s,1), stripPos(s,2)));
    
    % Write an index file for this strip.
    indexName = sprintf('%s_index%s.txt', outName, stripName);
    out = fopen(indexName, 'w');
    fprintf(out, 'name %s\n', stripName);
    fprintf(out, 'position %f %f\n', stripPos(s,1), stripPos(s,2));
    fprintf(out, 'files\n');
    for f=1:winN
       fprintf(out, '%s\n', nameList{stripIndex(s,f)}); 
    end
    fclose(out);
    
    % Figure out the most data we can use.
    dataNMin = 1e100;
    for w=1:winN
        fileName = sprintf('%s%s', inPrefix, nameList{stripIndex(s,w)});
        data = dlmread(fileName, ' ');
        
        dataN = length(data(:,1));
        if dataN < dataNMin
            dataNMin = dataN;
        end
        disp(sprintf('Read %d points from %s.', dataN, fileName));
    end
    
    % Use the most data we can.
    steps = dataNMin-cutIndex+1;
    X = zeros(steps,winN);
    disp(sprintf('Using %d points.', steps));
    for w=1:winN
        fileName = sprintf('%s%s', inPrefix, nameList{stripIndex(s,w)});
        data = dlmread(fileName, ' ');
        
        X(:,w) = data(end-steps+1:end,4);
        %disp(sprintf('Read %d points from %s.', dataN, fileName));
    end
    
    % Run WHAM.
    disp('Beginning WHAM iterations');
    F = WHAM_F(X,center,spring,1);
    disp('WHAM iterations complete.');
   
    % Write the free energy constants.
    outFileF = sprintf('%s_F%s.dat', outName, stripName);
    dlmwrite(outFileF, F, ' ');
    
    % Set the zero of the pmf.
    if makePmf
        % Compute the PMF.
        [u,z,p] = WHAMDIST(F,X,X,center,spring,ones(size(X)),1,0.1);
        
        farZ1 = max(z) - levelDist;
        farZ0 = farZ1 - levelDist;
        levelSum = 0;
        levelN = 0;
        for k=1:length(z)
            if z(k) >= farZ0 && z(k) < farZ1
                levelSum = levelSum + u(k);
                levelN = levelN + 1;
            end
        end
        levelU = levelSum/levelN;
        u = u - levelU;
        disp(sprintf('Shifted the PMF by %f kT.', levelU));
        disp(sprintf('The minimum of the PMF is %f kT.', min(u)));
        
        % Write the profile.
        outFile = sprintf('%s_pmf%s.dat', outName, stripName);
        dlmwrite(outFile, [z' u], ' ');
    end
end


