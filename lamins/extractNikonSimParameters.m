function [ s, T ] = extractNikonSimParameters( MD )
% extractNikonSimParameters Pulls the Nikon SIM reconstruction parameters
% from a MovieData structure metadata. This will probably work best on
% MovieData objects that use .nd2 files as a base

% if no argument given, scan current direcftory for .nd2 files
if(nargin < 1)
    D = dir('*.nd2');
    MD(length(D)) = MovieData;
    for i=1:length(D);
        MD(i) = MovieData.load([pwd filesep D(i).name]);
    end;
end
% if MD is a an array of MovieData objects, convert to a cell
if(~isscalar(MD))
    MD = num2cell(MD);
end
% if MD is a MovieList, convert to a cell
if(isa(MD,'MovieList'))
    if(isempty(MD.movies_))
        MD.sanityCheck;
    end
    MD = MD.movies_;
end
% if MD is a cell, then call function on each instance
% optionally, also output a table
if(iscell(MD))
    s = cellfun(@extractNikonSimParameters,MD,'UniformOutput',false);
    if(nargout > 1)
        s = s(~cellfun('isempty',s));
        T = struct2table([s{:}]);
    end
    return;
end


R = MD.getReader();
M = R.formatReader.getGlobalMetadata;
keys = M.keySet.toArray;
keysStr = arrayfun(@(x) x,keys,'Unif',false);
keyIdx = ~cellfun('isempty',regexp(keysStr,'UseChannel'));
s = [];
if(any(keyIdx))
    s = struct('movieDataPath_',MD.movieDataPath_,'movieDataFileName_',MD.movieDataFileName_);
    mandatoryFields = {'UseChannel','Objective','GratingPitch','ChannelWave','ContrastLow','IllumConstrst','NoiseSuppression','BlurSuppression'};
    for f = 1:length(mandatoryFields)
        s.(mandatoryFields{f}) = [];
    end
    reconstructionString = keysStr{keyIdx};
    parts = strsplit(reconstructionString, {' ',':'});
    for p = 1:2:length(parts)-1
%         if(isfield(s,parts{p}))
%             parts{p} = [parts{p} '_' num2str(p)];
%         end
        s.(parts{p}) = parts{p+1};
    end
else
    warning([MD.movieDataPath_ ' does not contain reconstruction parameters']);
end

if(nargout > 1)
    T = struct2table(s);
end

end

