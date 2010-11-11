function [data] = runDetection(data, overwrite)
% runDetection automatically runs Henry's detection software - with iclean
% 1 on all the movies specified in the experiment structure
% SYNOPSIS [data] = loadAndSaveDetection(data);
%
% INPUT     data      : experiment structure, containing a '.source' field (path to data location)
%           overwrite : optional, 1 to overwrite previous detection results (default 0).
%
% OUTPUT    none

% Francois Aguet, April 2010

if nargin<2
    overwrite = 0;
end

% loop over all entries in the structure to enter the image data necessary for the detection input
nExp = length(data);

parfor i = 1:nExp
    tifFiles = dir([data(i).source '*.tif*']);
    if isempty(tifFiles)
        error(['No *.TIF frames found in ' data(i).source]);
    else
        if ~(exist([data(i).source 'Detection'], 'dir') == 7)
            spotDetection(data(i).source, 1);
        elseif (overwrite)
            fprintf('Overwriting detection results for movie %d.\n', i);
            spotDetection(data(i).source, 1);
        end
    end
end