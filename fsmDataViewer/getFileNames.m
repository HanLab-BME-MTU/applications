function [path, fileNames, status] = getFileNames(fileName)
% 'getFileNames' gets the list of files that have the same body name
% and are located in the same directory that 'fileName'.
% TODO: use this function everywhere in qfsm.

if isempty(fileName)
    path = '';
    fileNames = {};
    status = 1;
    return;
end

[path, body, no, ext] = getFilenameBody(fileName);
fileNames = dir([path, filesep, body, '*', ext]);

% Check if there is any file.
if isempty(fileNames)
    status = 0;
    error(['No file name containing ' body ' can be found in ' path '.']);
end

% Rearrange files according to their number
filesOrder = zeros(length(fileNames), 1);
for i = 1:length(fileNames)
    [dummy1, dummy2, no] = getFilenameBody(fileNames(i).name);
    filesOrder(i) = str2double(no);
end
minNumber = min(filesOrder(:));
filesOrder = filesOrder - minNumber + 1;
fileNames = {fileNames(filesOrder).name};

status = 1;

end