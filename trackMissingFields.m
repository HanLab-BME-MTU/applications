function [ex] = trackMissingFields(ex)
% track the missing fields in the experiment

for i=1:length(ex)
    
    % BEFORE doing anything else, check if the data needs to be tracked -
    % the decision is made based on whether tracking data  - the lftInfo
    % field already exits or not
    
    trackVar = 0;
    
    if ~isfield(ex,'lftInfo')
        trackVar = 1;
    else
        if isempty(ex(i).lftInfo)
            trackVar = 1;
        end
    end
    
   
    % if data needs to be tracked, then proceed here:
    
    if trackVar == 1
        
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
        
        track_det = 0;
        % if a directory called TrackInfoMatrices already exists, go
        % there, if not create one
        if exist('TrackInfoMatrices')==7
            cd('TrackInfoMatrices');
            % if a trackInfo file already exists in this directory, read
            % it, else track
            if exist(saveResultsFilename)==2
                readTrackInfo = open(saveResultsFilename);
                if isfield(readTrackInfo,'trackInfo')==1
                    trackInfo = readTrackInfo.trackInfo;
                elseif isfield(readTrackInfo,'trackInfoMat')==1
                    trackInfo = readTrackInfo.trackInfoMat;
                end
            else
                track_det = 1;
            end
        else
            mkdir('TrackInfoMatrices');
            cd('TrackInfoMatrices');
            track_det = 1;
        end

        if track_det==1
            disp(['tracking movie #',num2str(i)]);
            [trackNum,trackInfo,errFlag] = trackWithGapClosing(movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam);
            trackInfo(isnan(trackInfo))=0;
            trackInfo = sparse(trackInfo);
            save(saveResultsFilename,'trackInfo');
        end
        
                
        cd(od);
%         % fill in trackInfo for the lifetime determination
%         ex(i).trackInfo = trackInfo;
%         ex(i).movieLength = length(movieInfo);
%         
%         if length(ex)>1
%             ex(i) = fillStructLifetimeInfo(ex(i));
%         else
%             ex = fillStructLifetimeInfo(ex);
%         end
%         
%         %... then delete it again to save space
%         ex(i).trackInfo = [];
    end
    
end % of for-loop

end % of function