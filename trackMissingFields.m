function [ex] = trackMissingFields(ex,force)
% track the missing fields in the experiment


if nargin < 2 || isempty(force)
    force = 0;
end

for i=1:length(ex)
    
    % if data needs to be tracked, then proceed here:
    cd(ex(i).source)
    if exist('TrackInfoMatrices','dir') == 0 || force
        
        costMatrices    = ex(i).tracksettings.costMat;
        gapCloseParam   = ex(i).tracksettings.gapClosePar;
        iterParam       = ex(i).tracksettings.iterPar;
        
        % now we're missing the variable movieInfo, which is the detection
        % data. If a valid detection structure is a field in ex, read it,
        % else set the load_det parameter to 1 so that the structure is
        % loaded from the appropriate detection file from the source path
        load_det = 0;
        if isfield(ex,'detection')
            detection = ex(i).detection;
            if ~isempty(detection)
                movieInfo = detection;
            else
                load_det = 1;
            end
        else
            load_det = 1;
        end
        
        % set appropriate saving and reading file names
        saveResultsDir  = ex(i).source;
        saveResultsFilename = 'trackInfo.mat';
        readDetectionFilename = 'detection.mat';
        
        
        % now load the detection data if necessary, using the appropriate
        % name
        if load_det == 1
            od = cd;
            cd(saveResultsDir);
            cd('DetectionStructures');
            loadfile = load(readDetectionFilename);
            if isfield(loadfile,'detection')
                movieInfo = loadfile.detection;
            elseif isfield(loadfile,'cdet')
                movieInfo = loadfile.cdet;
            else
                error('no detection data file of specified format found');
            end
            
            cd(od);
        end
            
        od = cd;
        cd(saveResultsDir);
        
        mkdir('TrackInfoMatrices');
        cd('TrackInfoMatrices');
        
        disp(['tracking movie #',num2str(i)]);
        [trackNum,trackInfo,errFlag] = trackWithGapClosing(movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam);
        trackInfo(isnan(trackInfo))=0;
        trackInfo = sparse(trackInfo);
        save(saveResultsFilename,'trackInfo');
        
        cd(od);
        
    else
        disp(['movie ' num2str(i) ' was skipped because it has already been tracked'])
    end
    
end % of for-loop

end % of function