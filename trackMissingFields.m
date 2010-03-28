function [exp] = trackMissingFields(exp,overwrite)
%
%
% Francois Aguet, Jan 2010

if nargin < 2 || isempty(overwrite)
    overwrite = 0;
end

for i = 1:length(exp)
               
    if ~(exist([exp(i).source 'TrackInfoMatrices' filesep 'trackInfo.mat'], 'file')==2) || overwrite
        
        costMatrices    = exp(i).tracksettings.costMat;
        gapCloseParam   = exp(i).tracksettings.gapClosePar;
        iterParam       = exp(i).tracksettings.iterPar;
        
        % now we're missing the variable movieInfo, which is the detection
        % data. If a valid detection structure is a field in exp, read it,
        % else load the structure from the appropriate detection file
        if isfield(exp(i), 'detection') && ~isempty(exp(i).detection)
            movieInfo = exp(i).detection;
        else
            loadfile = load([exp(i).source 'DetectionStructures' filesep 'detection.mat']);
            if isfield(loadfile, 'detection')
                movieInfo = loadfile.detection;
            elseif isfield(loadfile, 'cdet')
                movieInfo = loadfile.cdet;
            else
                error('No detection data file of specified format found');
            end
        end
        if ~(exist([exp(i).source 'TrackInfoMatrices'], 'dir')==7)
            mkdir([exp(i).source 'TrackInfoMatrices']);
        end;
        fprintf('Tracking movie no. %d\n', i);
        [trackNum,trackInfo] = trackWithGapClosing(movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam);
        trackInfo(isnan(trackInfo))=0;
        trackInfo = sparse(trackInfo);
        save([exp(i).source 'TrackInfoMatrices' filesep 'trackInfo.mat'], 'trackInfo');
    else
        fprintf('Movie no. %d was skipped because it has already been tracked\n', i);
    end 
end