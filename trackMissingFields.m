function [data] = trackMissingFieldsNEW(data,overwrite)
%
%
% Francois Aguet, Jan 2010

if nargin < 2 || isempty(overwrite)
    overwrite = 0;
end

for k = 1:length(data)
               
    if ~(exist([data(k).source 'TrackInfoMatrices'], 'dir')==7) || overwrite
        
        costMatrices    = data(k).tracksettings.costMat;
        gapCloseParam   = data(k).tracksettings.gapClosePar;
        iterParam       = data(k).tracksettings.iterPar;
        
        % now we're missing the variable movieInfo, which is the detection
        % data. If a valid detection structure is a field in data, read it,
        % else load the structure from the appropriate detection file
        if isfield(data(k), 'detection') && ~isempty(data(k).detection)
            movieInfo = data(k).detection;
        else
            loadfile = load([data(k).source 'Detection' filesep 'detectionResults.mat']);
            if isfield(loadfile, 'frameInfo')
                movieInfo = loadfile.frameInfo;
            else
                error('No detection data file of specified format found');
            end
        end
        if ~(exist([data(k).source 'TrackInfoMatrices'], 'dir')==7)
            mkdir([data(k).source 'TrackInfoMatrices']);
        end;
        fprintf('Tracking movie no. %d\n', k);
        saveResults.dir = [data.source 'TrackInfoMatrices' filesep] ;
        [trackNum,trackInfo] = trackWithGapClosing(movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam,saveResults);
        trackInfo(isnan(trackInfo))=0;
        trackInfo = sparse(trackInfo);
        save([data(k).source 'TrackInfoMatrices' filesep 'trackInfo.mat'], 'trackInfo');
    else
        fprintf('Movie no. %d was skipped because it has already been tracked\n', k);
    end 
end