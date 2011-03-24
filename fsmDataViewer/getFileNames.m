function [path, fileNames] = getFileNames(fileName)
% 'getFileNames' gets the list of files that have the same body name
% and are located in the same directory that 'fileName'.

if isempty(fileName)
    path = '';
    fileNames = {};
    return;
end

[path, body, ~, ext] = getFilenameBody(fileName);
fileNames = dir([path, filesep, body, '*', ext]);

% Check if there is any file.
if isempty(fileNames)
    error(['No file name containing ' body ' can be found in ' path '.']);
end

fileNames = {fileNames(:).name};

end