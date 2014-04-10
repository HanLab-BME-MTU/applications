function [ iterInfo ] = readIterationInfo(filename)
%READITERATIONINFO Read information from a txt file.

% Open text file
fileIter = fopen(filename, 'rt');
if fileIter == -1
    error('Could not open file "%s".\n', filename);
end

% Read lines (skip first one)
fgetl(fileIter);
no = 1;
while ~feof(fileIter)
    % Read line
    line = fgetl(fileIter);
    
    % Scan numeric values
    V = sscanf(line, '%f');
    iterInfo(no,:) = V(:); %#ok<AGROW>
    no = no+1;
end
fclose(fileIter);

end
