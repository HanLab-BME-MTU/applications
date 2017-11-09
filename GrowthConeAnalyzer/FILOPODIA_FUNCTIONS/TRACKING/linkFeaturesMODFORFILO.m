function [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
    nnDistFeatures,prevCost,errFlag] = linkFeatureMODFORFILO(movieInfo,...
    costMatParams,analInfo,kalmanFunctions,probDim,filterInfoPrev,...
    prevCost,verbose,saveDir)
%LINKFEATURESKALMAN links features between consecutive frames using LAP and possibly motion propagation using the Kalman filter
%
%SYNOPSIS [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
%    nnDistFeatures,prevCost,errFlag] = linkFeaturesKalman(movieInfo,...
%    costMatName,costMatParam,kalmanFunctions,probDim,filterInfoPrev,...
%    prevCost,verbose)
%
%INPUT  movieInfo      : Array of size equal to the number of frames
%                        in a movie, containing the fields:
%             .xCoord      : x-coordinates of detected features. 
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .zCoord      : z-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                            Optional. Skipped if problem is 2D. Default: zeros.
%             .amp         : "Intensities" of detected features.
%                            1st column: values (ones if not available),
%                            2nd column: standard deviation (zeros if not
%                            available).
%             .num         : Number of features. 
%                            Optional. Calculated if not supplied.  
%             .allCoord    : x,dx,y,dy,[z,dz] of features collected in one
%                            matrix. Optional. Calculated if not supplied.
%             .nnDist      : Distance from each feature to its nearest
%                            neighbor. Optional. Calculated if not supplied.
%       costMatName    : Name of cost matrix function used for linking.
%       costMatParam   : Parameters needed for cost matrix calculation. 
%                        Structure with fields specified by particular
%                        cost matrix used (costMatName).
%       kalmanFunctions: Names of Kalman filter functions for self-adaptive
%                        tracking. Structure with fields:
%             .reserveMem   : Reserves memory for kalmanFilterInfo.
%             .initialize   : Initializes the Kalman filter for an appearing
%                             feature.
%             .calcGain     : Calculates the Kalman gain after linking.
%                        For non-self-adaptive tracking, enter [].
%                        Optional. Default: [].
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%                        Optional. If not input, dimensionality will be
%                        derived from movieInfo.
%       filterInfoPrev : Structure array with number of entries equal to
%                        number of frames in movie. Contains at least the
%                        fields:
%             .stateVec    : Kalman filter state vector for each feature in frame.
%             .stateCov    : Kalman filter state covariance matrix for each feature in frame.
%             .noiseVar    : Variance of state noise for each feature in frame.
%                            Optional. Enter [] or nothing if not to be used.
%       prevCost       : Matrix of costs of actual assignments in previous
%                        round of linking. Optional. Default: Empty.
%       verbose        : 1 to show calculation progress, 0 otherwise.
%                        Optional. Default: 1.
%
%OUTPUT trackedFeatureIndx: Connectivity matrix of features between frames.
%                           Rows indicate continuous tracks, while columns 
%                           indicate frames. A track that ends before the
%                           last frame is followed by zeros, and a track
%                           that starts in a frame after the first frame
%                           is preceded by zeros. 
%       trackedFeatureInfo: The positions and "intensities" of the tracked
%                           features, in the same units as input. 
%                           Number of rows = number of tracks.
%                           Number of columns = 8*number of frames.
%                           Each row consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           NaN is used to indicate frames where the tracks
%                           do not exist.
%       kalmanFilterInfo  : Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields
%                           defined in kalmanFunctions.reserveMem (at
%                           least stateVec, stateCov and noiseVar).
%       nnDistFeatures    : Matrix indicating the nearest neighbor
%                           distances of features linked together within
%                           tracks.
%       prevCost          : Matrix of costs of actual assignments.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must not be empty.
%
%Khuloud Jaqaman, March 2007
% MODIFIED SLIGHTLY for QUICK FIX FOR consistency with filo object output
% extra input 
% analInfo: analInfo(iFrame).filoInfo(ifilo).xDescriptor (include all xy
% coords of filo)
%% Output

trackedFeatureIndx = [];
trackedFeatureInfo = [];
kalmanFilterInfo = [];
nnDistFeatures = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 5
    disp('--linkFeaturesKalman: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether tracking is self-adaptive 
if nargin < 4 || isempty(kalmanFunctions)
    kalmanFunctions = [];
    selfAdaptive = 0;
else
    selfAdaptive = 1;
end

%check whether z-coordinates were input, making problem potentially 3D
if isfield(movieInfo,'zCoord')
    probDimT = 3;
else    
    probDimT = 2;
end

%assign problem dimensionality if not input
if nargin < 5 || isempty(probDim)
    probDim = probDimT;
else
    if probDim == 3 && probDimT == 2
        disp('--linkFeaturesKalman: Inconsistency in input. Problem 3D but no z-coordinates.');
        errFlag = 1;
    end
end

%check whether a priori Kalman filter information is given
if nargin < 6 || isempty(filterInfoPrev)
    filterInfoPrev = [];
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%check whether previous costs have been input
if nargin < 7 || isempty(prevCost)
    prevCost = [];
end

%check whether verbose
if nargin < 8 || isempty(verbose)
    verbose = 1;
end

%exit if there are problems with input
if errFlag
    disp('--linkFeaturesKalman: Please fix input parameters.');
    return
end



%% preamble

%get number of frames in movie
numFrames = length(movieInfo);

%get number of features in each frame
if ~isfield(movieInfo,'num')
    for iFrame = 1 : numFrames
        movieInfo(iFrame).num = size(movieInfo(iFrame).xCoord,1);
    end
end

%collect coordinates and their std in one matrix in each frame
if ~isfield(movieInfo,'allCoord')
    if probDim == 2
        for iFrame = 1 : numFrames
            movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                movieInfo(iFrame).yCoord];
        end
    else
        for iFrame = 1 : numFrames
            movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                movieInfo(iFrame).yCoord movieInfo(iFrame).zCoord];
        end
    end
end

%calculate nearest neighbor distance for each feature in each frame
if ~isfield(movieInfo,'nnDist')

    for iFrame = 1 : numFrames
        
        switch movieInfo(iFrame).num

            case 0 %if there are no features

                %there are no nearest neighbor distances
                nnDist = zeros(0,1);

            case 1 %if there is only 1 feature

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000;

            otherwise %if there is more than 1 feature

                %compute distance matrix
                nnDist = createDistanceMatrix(movieInfo(iFrame).allCoord(:,1:2:end),...
                    movieInfo(iFrame).allCoord(:,1:2:end));

                %sort distance matrix and find nearest neighbor distance
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,2);

        end

        %store nearest neighbor distance
        movieInfo(iFrame).nnDist = nnDist;

    end

end

%% Linking

%make an array of the number of features per frame
numFeatures = zeros(numFrames,1);
for iFrame = 1 : numFrames
    numFeatures(iFrame) = movieInfo(iFrame).num;
end

%reserve memory for kalmanFilterInfo
if selfAdaptive

    % -- USER DEFINED FUNCTION -- %
    eval(['kalmanFilterInfo = ' kalmanFunctions.reserveMem ...
        '(numFrames,numFeatures,probDim);']);
    
else
    
    kalmanFilterInfo = zeros(numFrames,1);

end

%fill the feature indices in 1st frame in the connectivity matrix
trackedFeatureIndx = (1:movieInfo(1).num)';

%fill the nearest neighbor distances of features in first frame
nnDistFeatures = movieInfo(1).nnDist;

%initialize Kalman filter for features in 1st frame
if selfAdaptive

    % -- USER DEFINED FUNCTION -- %
    if usePriorInfo %use a priori information if available
        kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
        kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
        kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
    else
        eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
            '(movieInfo(1),probDim,costMatParam);'])
        kalmanFilterInfo(1).stateVec = filterInit.stateVec;
        kalmanFilterInfo(1).stateCov = filterInit.stateCov;
        kalmanFilterInfo(1).noiseVar = filterInit.noiseVar;
    end
    
end

%store the costs of previous links for features in first frame
%if no previous costs have been input
%in this case, store NaN since there are no previous links
if isempty(prevCost)
    prevCost = NaN(movieInfo(1).num,1);
else
    prevCost = max(prevCost(:))*ones(movieInfo(1).num,1);
end

%assign the lifetime of features in first frame
featLifetime = ones(movieInfo(1).num,1);

%initialize progress display
if verbose
    progressText(0,'Linking frame-to-frame');
end

% % % %for paper - get number of potential link per feature
% % % numPotLinksPerFeature = [];

%go over all frames
for iFrame = 1 : numFrames-1

    %get number of features
    numFeaturesFrame1 = movieInfo(iFrame).num; %in 1st frame
    numFeaturesFrame2 = movieInfo(iFrame+1).num; %in 2nd frame

    if numFeaturesFrame1 ~= 0 %if there are features in 1st frame

        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
            
            %% quick fix- MB % calculate cost Matrix for iframe --> iframe +1
           
            % extract info for the two frames to be linked 
            currentFrameInfo = analInfo(iFrame:iFrame+1); 
     
                   troubleshoot  = 0; 
             %[costMat,nonlinkMarker]
             %=costMatFilo(currentFrameInfo,neuriteEdgeT,costMatParams,troubleshoot,iFrame,saveDir);%
             %20140915 look up old costMat 
           
             [costMat,nonlinkMarker] = costMatFilo(currentFrameInfo,costMatParams,iFrame,saveDir); 
             
                 
%             calculate cost matrix
%             % -- USER DEFINED FUNCTION -- %
%             eval(['[costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker]'...
%                 ' = ' costMatName '(movieInfo,kalmanFilterInfo(iFrame),'...
%                 'costMatParam,nnDistFeatures(1:numFeaturesFrame1,:),'...
%                 'probDim,prevCost,featLifetime,trackedFeatureIndx,iFrame);'])

            % % %             %for paper - get number of potential links per feature
            % % %             numPotLinksPerFeature = [numPotLinksPerFeature; sum(...
            % % %                 costMat(1:numFeaturesFrame1,1:numFeaturesFrame2)...
            % % %                 ~=nonlinkMarker,2)];

            if any(costMat(:)~=nonlinkMarker) %if there are potential links

                %link features based on cost matrix, allowing for birth and death
                [link12,link21] = lap(costMat,nonlinkMarker,0);

                %get indices of features in 2nd frame that are connected to features in 1st frame
                indx2C = find(link21(1:numFeaturesFrame2)<=numFeaturesFrame1);

                %get indices of corresponding features in 1st frame
                indx1C = link21(indx2C);

                %find existing tracks that are not connected to features in 2nd frame
                indx1U = ones(size(trackedFeatureIndx,1),1);
                indx1U(indx1C) = 0;
                indx1U = find(indx1U);

                %assign space for new connectivity matrix,
                %for matrix of nearest neighbors of connected features,
                %and for matrix of linking costs
                tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2-length(indx2C),iFrame+1);
                tmpNN = NaN(size(tmp));
                tmpCost = NaN(size(tmp));

                %fill in the feature numbers in 2nd frame and their nearest
                %neighbor distances
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
                tmpNN(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

                %fill in the linking costs for features in 2nd frame that
                %got linked to features in 1st frame
                for i = 1 : length(indx2C)
                    tmpCost(indx2C(i),iFrame+1) = costMat(indx1C(i),indx2C(i));
                end
                
                %shuffle existing tracks to get the correct connectivity with 2nd frame
                %also shuffle nearest neighbor distances and previous linking costs
                tmp(indx2C,1:iFrame) = trackedFeatureIndx(indx1C,:);
                tmpNN(indx2C,1:iFrame) = nnDistFeatures(indx1C,:);
                tmpCost(indx2C,1:iFrame) = prevCost(indx1C,:);

                %add rows of tracks that are not connected to points in 2nd frame
                %also add nearest neighbor distances and previous linking costs
                tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx(indx1U,:);
                tmpNN(numFeaturesFrame2+1:end,1:iFrame) = nnDistFeatures(indx1U,:);
                tmpCost(numFeaturesFrame2+1:end,1:iFrame) = prevCost(indx1U,:);

                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;
                
                %update the nearest neighbor distances of connected features
                nnDistFeatures = tmpNN;
                
                %update the matrix of previous linking costs
                prevCost = tmpCost;
                
                %get track lifetimes for features in 2nd frame
                featLifetime = ones(numFeaturesFrame2,1);
                for iFeat = 1 : numFeaturesFrame2
                    featLifetime(iFeat) = length(find(trackedFeatureIndx(iFeat,:)~=0));
                end

                %use the Kalman gain from linking to get better estimates
                %of the state vector and its covariance matrix in 2nd frame
                %as well as state noise and its variance
                if selfAdaptive

                    % -- USER DEFINED FUNCTION -- %
                    if usePriorInfo %if prior information is supplied
                        eval(['[kalmanFilterInfo,errFlag] = ' kalmanFunctions.calcGain ...
                            '(trackedFeatureIndx(1:numFeaturesFrame2,:),'...
                            'movieInfo(iFrame+1),kalmanFilterInfoTmp,'...
                            'propagationScheme,kalmanFilterInfo,probDim,'...
                            'filterInfoPrev(iFrame+1),costMatParam,kalmanFunctions.initialize);'])
                    else %if no prior information is supplied
                        eval(['[kalmanFilterInfo,errFlag] = ' kalmanFunctions.calcGain ...
                            '(trackedFeatureIndx(1:numFeaturesFrame2,:),'...
                            'movieInfo(iFrame+1),kalmanFilterInfoTmp,'...
                            'propagationScheme,kalmanFilterInfo,probDim,'...
                            '[],costMatParam,kalmanFunctions.initialize);'])
                    end
                    
                end
                
            else %if there are no potential links

                %assign space for new connectivity matrix,
                %for matrix of nearest neighbors of connected features,
                %and for matrix of linking costs
                tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2,iFrame+1);
                tmpNN = NaN(size(tmp));
                tmpCost = NaN(size(tmp));

                %fill in the feature numbers in 2nd frame and their nearest
                %neighbor distances
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
                tmpNN(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

                %fill in the tracks upto 1st frame, their nearest
                %neighbor distances and their linking costs
                tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx;
                tmpNN(numFeaturesFrame2+1:end,1:iFrame) = nnDistFeatures;
                tmpCost(numFeaturesFrame2+1:end,1:iFrame) = prevCost;

                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;

                %update the nearest neighbor distances of connected features
                nnDistFeatures = tmpNN;

                %update the matrix of previous linking costs
                prevCost = tmpCost;

                %assign track lifetimes for features in 2nd frame
                featLifetime = ones(numFeaturesFrame2,1);

                %initialize Kalman filter for features in 2nd frame
                if selfAdaptive

                    % -- USER DEFINED FUNCTION -- %
                    if usePriorInfo %use a priori information if available
                        kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                        kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                        kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                    else
                        eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
                            '(movieInfo(iFrame+1),probDim,costMatParam);'])
                        kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                        kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                        kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                    end
                    
                end

            end %(if any(costMat(:)~=nonlinkMarker))

        else %if there are no features in 2nd frame

            %add a column of zeros for the 2nd frame in the track
            %connectivity matrix
            trackedFeatureIndx = [trackedFeatureIndx zeros(size(trackedFeatureIndx,1),1)];
            
            %add a column of NaNs for the 2nd frame in the nearest neighbor
            %matrix
            nnDistFeatures = [nnDistFeatures NaN(size(nnDistFeatures,1),1)];
            
            %add a column of NaNs for the 2nd frame in the matrix of
            %previous linking costs
            prevCost = [prevCost NaN(size(trackedFeatureIndx,1),1)];

            %assign track lifetimes for features in 2nd frame
            featLifetime = [];

        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    else %if there are no feature in 1st frame

        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame

            %assign space for new connectivity matrix,
            %for matrix of nearest neighbors of connected features,
            %and for matrix of linking costs
            tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2,iFrame+1);
            tmpNN = NaN(size(tmp));
            tmpCost = NaN(size(tmp));

            %fill in the feature numbers in 2nd frame and their nearest
            %neighbor distances
            tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
            tmpNN(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

            %fill in the tracks upto 1st frame, their nearest
            %neighbor distances and their linking costs
            tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx;
            tmpNN(numFeaturesFrame2+1:end,1:iFrame) = nnDistFeatures;
            tmpCost(numFeaturesFrame2+1:end,1:iFrame) = prevCost;

            %update the connectivity matrix "trackedFeatureIndx"
            trackedFeatureIndx = tmp;

            %update the nearest neighbor distances of connected features
            nnDistFeatures = tmpNN;
            
            %update the matrix of previous linking costs
            prevCost = tmpCost;

            %assign track lifetimes for features in 2nd frame
            featLifetime = ones(numFeaturesFrame2,1);

            %initialize Kalman filter for features in 2nd frame
            if selfAdaptive

                % -- USER DEFINED FUNCTION -- %
                if usePriorInfo %use a priori information if available
                    kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                    kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                else
                    eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
                        '(movieInfo(iFrame+1),probDim,costMatParam);'])
                    kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                    kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                end
                
            end

        else %if there are no features in 2nd frame

            %add a column of zeros for the 2nd frame in the track
            %connectivity matrix
            trackedFeatureIndx = [trackedFeatureIndx zeros(size(trackedFeatureIndx,1),1)];
            
            %add a column of NaNs for the 2nd frame in the nearest neighbor
            %matrix
            nnDistFeatures = [nnDistFeatures NaN(size(nnDistFeatures,1),1)];

            %add a column of NaNs for the 2nd frame in the matrix of
            %previous linking costs
            prevCost = [prevCost NaN(size(trackedFeatureIndx,1),1)];

            %assign track lifetimes for features in 2nd frame
            featLifetime = [];

        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    end %(if numFeaturesFrame1 ~= 0 ... else ...)

    %display progress
    if verbose
        progressText(iFrame/(numFrames-1),'Linking frame-to-frame');
    end

end %(for iFrame=1:numFrames-1)

%get total number of tracks
numTracks = size(trackedFeatureIndx,1);

%find the frame where each track begins and then sort the vector
frameStart = zeros(numTracks,1);
for i=1:numTracks
    frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
end
[frameStart,indx] = sort(frameStart);

%rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the 
%same time in descending order from longest to shortest.
trackedFeatureIndx = trackedFeatureIndx(indx,:);

%also re-arrange the matrix indicating nearest neighbor distances
nnDistFeatures = nnDistFeatures(indx,:);

%clear some memory
clear costMat tmp tmpNN tmpCost

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system

%reserve space for matrix
trackedFeatureInfo = NaN*ones(size(trackedFeatureIndx,1),8*numFrames);

%for now, the matrix always has space for z, hence I have to treat 2D and
%3D differently.
%at some point I can change the postprocessing functions to handle 2D and
%3D differently, and so I can save them differently, in which case I won't
%need the if-else here.
if probDim ==2

    %go over all frames
    for iFrame = 1 : numFrames

        %find rows that have a feature index
        indx1 = find(trackedFeatureIndx(:,iFrame)~=0);

        %if there are detected features in this frame
        if ~isempty(indx1)

            %find the feature index
            indx2 = trackedFeatureIndx(indx1,iFrame);

            %store its information
%             trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
%                 [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
%                 zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,1) ...
%                 movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
%                 zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,2)];

trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
             [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
             zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,1) ...
                 movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
                 zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,2)];


        end

    end

else

    %go over all frames
    for iFrame = 1 : numFrames

        %find rows that have a feature index
        indx1 = find(trackedFeatureIndx(:,iFrame)~=0);

        %if there are detected features in this frame
        if ~isempty(indx1)

            %find the feature index
            indx2 = trackedFeatureIndx(indx1,iFrame);

            %store its information
            trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
                [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
                movieInfo(iFrame).amp(indx2,1) ...
                movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
                movieInfo(iFrame).amp(indx2,2)];

        end

    end

end

%% %%%%% ~~ the end ~~ %%%%%


%old stuff
%
% nnDistTracks = movieInfo(1).nnDist;
%                 tmpNN2(indx2C) = nnDistTracks(indx1C);
%                 tmpNN2(numFeaturesFrame2+1:end) = nnDistTracks(indx1U);
%                 nnDistTracks = tmpNN2;
%                 tmpNN = max(1,iFrame+1-nnWindow);
%                 nnDistTracks(1:numFeaturesFrame2) = min(nnDistFeatures...
%                     (1:numFeaturesFrame2,tmpNN:end),[],2);
%                 tmpNN2(numFeaturesFrame2+1:end) = nnDistTracks;
%                 nnDistTracks = tmpNN2;
%                 nnDistTracks(1:numFeaturesFrame2) = movieInfo(iFrame+1).nnDist;
%             tmpNN2(numFeaturesFrame2+1:end) = nnDistTracks;
%             nnDistTracks = tmpNN2;
%             nnDistTracks(1:numFeaturesFrame2) = movieInfo(iFrame+1).nnDist;
