function [tubIntensity,errFlag] = corrMotionTypeWithMTs(trackData,...
    fullImagePathName,diffAnalysisData,minTrackLen,imageRange)


%% input

errFlag = 0;

if nargin < 3
    disp('--corrMotionTypeWithMTs: Missing input arguments!');
    return
end

if nargin < 4 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 5
    imageRange = [];
end

%% calculation

%get number of input movies
numMovies = length(trackData);

%go over all movies in trackData ...
for iMovie = 1 : numMovies
    
    %% preamble
    
    %get this movie's tracks and their diffusion analysis results
    tracks = trackData(iMovie).tracks;
    diffAnalysisRes = diffAnalysisData(iMovie).diffAnalysisRes;
    
    %keep only tracks with length >= minTrackLen
    criteria.lifeTime.min = minTrackLen;
    indx = chooseTracks(tracks,criteria);
    clear criteria
    tracks = tracks(indx);
    diffAnalysisRes = diffAnalysisRes(indx);
    
    %also retain only tracks in image region of interest
    if ~isempty(imageRange)
        criteria.initialXCoord.min = imageRange(1,1,iMovie);
        criteria.initialXCoord.max = imageRange(1,2,iMovie);
        criteria.initialYCoord.min = imageRange(2,1,iMovie);
        criteria.initialYCoord.max = imageRange(2,2,iMovie);
        indx = chooseTracks(tracks,criteria);
        clear criteria
        tracks = tracks(indx);
        diffAnalysisRes = diffAnalysisRes(indx);
    end
    
    %put tracks in matrix format
    tracksMat = convStruct2MatIgnoreMS(tracks);
    xCoordMat = tracksMat(:,1:8:end);
    yCoordMat = tracksMat(:,2:8:end);
    
    %load image
    cellImage = imread(fullImagePathName{iMovie});
    cellImage = double(cellImage);
    cellImage = (cellImage - min(cellImage(:))) / (max(cellImage(:)) - min(cellImage(:)));
    
    %% motion types
    
    %get track segment classification from diffusion analysis
    trackSegmentClass = vertcat(diffAnalysisRes.classification);
    
    %get track segment length
    trackSegmentLength = getTrackSEL(tracksMat);
    trackSegmentLength = trackSegmentLength(:,3);
    
    %get indices of linear, Brownian, confined Brownian and undetermined tracks
    indxLin    = find( trackSegmentClass(:,1) == 1 | trackSegmentClass(:,2) == 3 );
    indxBrown  = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 2 );
    indxConf   = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 1 );
    indxUndet1 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
        & trackSegmentLength >= 5);
    indxUndet2 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
        & trackSegmentLength < 5);
    
    %go over the tracks of each category and get the distribution of
    %tubulin intensities at their positions
    motionType = {'Lin','Brown','Conf','Undet1','Undet2'};
    for iType = 1 : 5
        
        %get track indices in this motion category
        eval(['indxType = indx' motionType{iType} ';']);
        
        %get coordinates, rounded up to integer pixel values
        xCoordType = ceil(xCoordMat(indxType,:));
        yCoordType = ceil(yCoordMat(indxType,:));
        
        %convert from 2D to 1D indices
        xCoordType = xCoordType(~isnan(xCoordType));
        yCoordType = yCoordType(~isnan(yCoordType));
        linCoordType = sub2ind(size(cellImage),yCoordType,xCoordType);
        
        %read out tubulin intensities at these pixels
        tubIntType = double(cellImage(linCoordType));
        
        %store for output
        eval(['tubIntensity(iMovie).' motionType{iType} ...
            '.values = tubIntType;'])
        eval(['tubIntensity(iMovie).' motionType{iType} ...
            '.meanStd = [mean(tubIntType) std(tubIntType)];'])
        
    end
    
    %plot tracks, for debugging purposes
    plotTracksDiffAnalysis2D(tracks,diffAnalysisRes,[],1,cellImage,0,1);
    
end

