%[dname dpath] = getCellDirectory(data) returns the name of the cell directory corresponding to 'data'.

% Francois Aguet, March 14 2011

function [dname dpath] = getCellDirectory(data, selector)

if nargin<2
    selector = 'cell';
end

nCh = length(data.channels);

masterChannel = find(strcmp(data.channels, data.source));
slaveChannels = setdiff(1:nCh, masterChannel);

if ~isempty(slaveChannels)
    dpath = getParentDir(data.channels{slaveChannels(1)});
    dname = getDirFromPath(dpath);
    dpath = getParentDir(dpath);
else
    dpath = data.source;
    dname = getDirFromPath(dpath);
    % verify that directory contains selector
    if ~isempty(strfind(dname, selector))
        error([data.source ' is not a valid movie directory.']);
    end
end