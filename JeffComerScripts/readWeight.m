function weight = readWeight(logFile) 
% Author: Jeff Comer <jcomer2@illinois.edu>

in = fopen(logFile, 'r');
if in == -1
    error('readWeight:fileNotFound', 'The log file was not found.')
end

weight = 0;
while 1
    lin = fgets(in);
    
    if ~ischar(lin), break, end
    if isempty(lin)
        msg = ferror(in);
        error('readWeight:ioError', msg)
    end
    
    if (strcmp(lin(1:7), 'WEIGHT:')) 
       [prefix file weight] = strread(lin, '%s %s %d');
       break
    end
end

fclose(in);
