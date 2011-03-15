%[dname dpath] = getCellDirectory(data) returns the name of the cell directory corresponding to 'data'.

% Francois Aguet, March 14 2011

function [dname dpath] = getCellDirectory(data)

nCh = length(data.channels);

masterChannel = find(strcmp(data.channels, data.source));
slaveChannels = setdiff(1:nCh, masterChannel);

dpath = getParentDir(data.channels{slaveChannels(1)});
dname = getDirFromPath(dpath);
dpath = getParentDir(dpath);