%spath = getShortPath(data) returns the 3 directories above cell level

% Francois Aguet, 05/13/2011

function spath = getShortPath(data)

mCh = strcmp(data.channels, data.source);
sCh = setdiff(1:length(data.channels),mCh);
spath = data.channels{sCh(1)};
idx = regexp(spath, filesep);
spath = dpath(idx(end-4)+1:idx(end-1));