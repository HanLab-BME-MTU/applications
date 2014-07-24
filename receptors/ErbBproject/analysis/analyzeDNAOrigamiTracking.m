function success = analyzeDNAOrigamiTracking(file, name)
% analyzeDNAOrigamiTrackinh is a function that loads a mat file containing
% two lists of localizations and runs khoulouds tracking software on them.
%
%
%   inputs:
%           file, filename that contains point lists (point lists must be
%           named atto and cy3
%
%           name, string with the name of the file for saving the results
%
%   outputs: 
%            success

load(file);

fprintf(1,'Working on Atto');
analyzeDNAOrigamiPointListTracking(atto,[name '_atto_tracking.mat'])

analyzeDNAOrigamiPointListTracking(cy3,[name '_cy3_tracking.mat'])

end

