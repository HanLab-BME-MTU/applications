function [trackedFeatureIndx,trackedFeatureInfo,errFlag] = linkFeaturesSR(...
    movieInfo,searchRadius,probDim,verbose)
%LINKFEATURESSR links features between consecutive frames using LAP for super-resolution applications
%
%SYNOPSIS [trackedFeatureIndx,trackedFeatureInfo,nnDistFeatures,prevCost,...
%    errFlag] = linkFeaturesSR(movieInfo,costMatParam,probDim,...
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
%       searchRadius   : Search radius for linking features between
%                        cosecutive frames.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%                        Optional. If not input, dimensionality will be
%                        derived from movieInfo.
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
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must not be empty.
%
%Khuloud Jaqaman, August 2011

%% Output

trackedFeatureIndx = [];
trackedFeatureInfo = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--linkFeaturesSR: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether z-coordinates were input, making problem potentially 3D
if isfield(movieInfo,'zCoord')
    probDimT = 3;
else    
    probDimT = 2;
end

%assign problem dimensionality if not input
if nargin < 3 || isempty(probDim)
    probDim = probDimT;
else
    if probDim == 3 && probDimT == 2
        disp('--linkFeaturesSR: Inconsistency in input. Problem 3D but no z-coordinates.');
        errFlag = 1;
    end
end

%check whether verbose
if nargin < 4 || isempty(verbose)
    verbose = 1;
end

%exit if there are problems with input
if errFlag
    disp('--linkFeaturesSR: Please fix input parameters.');
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

%% Linking

%make an array of the number of features per frame
numFeatures = zeros(numFrames,1);
for iFrame = 1 : numFrames
    numFeatures(iFrame) = movieInfo(iFrame).num;
end

%fill the feature indices in 1st frame in the connectivity matrix
trackedFeatureIndx = (1:movieInfo(1).num)';

%get number of particles in whole movie and calculate a worst-case scenario
%number of tracks
%it can be that the final number of tracks is even larger than this worst
%case scenario. Every time the auxiliary matrices (defined below) run out
%of rows, another "numTracksWorstCase" rows are added to them.
numTracksWorstCase = round(sum(numFeatures)/10);

%initialize auxiliary matrices for storing information related to tracks
%that end in the middle of the movie
trackedFeatureIndxAux = zeros(numTracksWorstCase,numFrames);
rowEnd = numTracksWorstCase;

%initialize progress display
if verbose
    progressText(0,'Linking frame-to-frame');
end

%go over all frames
for iFrame = 1 : numFrames-1

    %get number of features
    numFeaturesFrame1 = movieInfo(iFrame).num; %in 1st frame
    numFeaturesFrame2 = movieInfo(iFrame+1).num; %in 2nd frame

    if numFeaturesFrame1 ~= 0 %if there are features in 1st frame
        
        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
            
            %calculate modified positions of features in 1st frame
            %basically the position is taken as the average from the
            %multiple appearances upto the first frame
            movieInfoTmp = movieInfo(iFrame:iFrame+1);
            %             coordInfoTmp = coordAmpMatFromIndicesSparse(trackedFeatureIndx,...
            %                 movieInfo(1:iFrame),iFrame,probDim);
            %             xCoordTmp = full(coordInfoTmp(:,1:8:end));
            %             xCoordTmp(xCoordTmp==0) = NaN;
            %             xCoordTmp = nanmean(xCoordTmp,2);
            %             movieInfoTmp(1).allCoord(:,1) = xCoordTmp;
            %             yCoordTmp = full(coordInfoTmp(:,2:8:end));
            %             yCoordTmp(yCoordTmp==0) = NaN;
            %             yCoordTmp = nanmean(yCoordTmp,2);
            %             movieInfoTmp(1).allCoord(:,3) = yCoordTmp;
            %             if probDim == 3
            %                 zCoordTmp = full(coordInfoTmp(:,3:8:end));
            %                 zCoordTmp(zCoordTmp==0) = NaN;
            %                 zCoordTmp = nanmean(zCoordTmp,2);
            %                 movieInfoTmp(1).allCoord(:,5) = zCoordTmp;
            %             end
            
            %calculate cost matrix
            [costMat,nonlinkMarker,errFlag] = costMatLinkSR(movieInfoTmp,...
                searchRadius);
            
            if any(costMat(:)~=nonlinkMarker) %if there are potential links

                %link features based on cost matrix, allowing for birth and death
                [link12,link21] = lap(costMat,nonlinkMarker,0);

                %get indices of features in 2nd frame that are connected to features in 1st frame
                indx2C = find(link21(1:numFeaturesFrame2)<=numFeaturesFrame1);

                %get indices of corresponding features in 1st frame
                indx1C = link21(indx2C);

                %find existing tracks that are not connected to features in 2nd frame
                numExistTracks = size(trackedFeatureIndx,1);
                indx1U = setdiff(1:numExistTracks,indx1C);
                numRows = length(indx1U);
                
                %determine where to store these tracks in auxiliary matrix
                %extend auxiliary matrices if necessary
                rowStart = rowEnd - numRows + 1;
                if rowStart <= 1
                    trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                        trackedFeatureIndxAux];
                    rowEnd = rowEnd + numTracksWorstCase;
                    rowStart = rowStart + numTracksWorstCase;
                end
                
                %move rows of tracks that are not connected to points in
                %2nd frame to auxilary matrix
                trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx(indx1U,:);
                
                %assign space for new connectivity matrix
                tmp = zeros(numFeaturesFrame2,iFrame+1);

                %fill in the feature numbers in 2nd frame
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';

                %shuffle existing tracks to get the correct connectivity with 2nd frame
                tmp(indx2C,1:iFrame) = trackedFeatureIndx(indx1C,:);
                
                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;

                %update rowEnd to indicate until which row the auxiliary
                %matrices are ampty
                rowEnd = rowStart - 1;

            else %if there are no potential links
                
                %determine where to store the tracks up to 1st frame in
                %auxiliary matrix
                %extend auxiliary matrices if necessary
                numRows = size(trackedFeatureIndx,1);
                rowStart = rowEnd - numRows + 1;
                if rowStart <= 1
                    trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                        trackedFeatureIndxAux];
                    rowEnd = rowEnd + numTracksWorstCase;
                    rowStart = rowStart + numTracksWorstCase;
                end
                
                %move tracks upto 1st frame to auxiliary matrix
                trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx;

                %assign space for new connectivity matrix
                trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);

                %fill in the feature numbers in 2nd frame
                trackedFeatureIndx(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
                
                %update rowEnd to indicate until which row the auxiliary
                %matrices are ampty
                rowEnd = rowStart - 1;

            end %(if any(costMat(:)~=nonlinkMarker))
            
        else %if there are no features in 2nd frame
            
            %determine where to store the tracks up to 1st frame in
            %auxiliary matrix
            %extend auxiliary matrices if necessary
            numRows = size(trackedFeatureIndx,1);
            rowStart = rowEnd - numRows + 1;
            if rowStart <= 1
                trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                    trackedFeatureIndxAux];
                rowEnd = rowEnd + numTracksWorstCase;
                rowStart = rowStart + numTracksWorstCase;
            end
            
            %move tracks upto 1st frame to auxiliary matrix
            trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx;
            
            %update the connectivity matrix "trackedFeatureIndx"
            trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
            %update rowEnd to indicate until which row the auxiliary
            %matrices are ampty
            rowEnd = rowStart - 1;
            
        end %(if numFeaturesFrame2 ~= 0 ... else ...)
        
    else %if there are no feature in 1st frame
        
        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
            
            %assign space for new connectivity matrix
            trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
            
            %fill in the feature numbers in 2nd frame
            trackedFeatureIndx(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
            
        else %if there are no features in 2nd frame

            %assign space for new connectivity matrix
            trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
            
        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    end %(if numFeaturesFrame1 ~= 0 ... else ...)

    %display progress
    if verbose
        progressText(iFrame/(numFrames-1),'Linking frame-to-frame');
    end
    
end %(for iFrame=1:numFrames-1)

%add information from last frame to auxiliary matrices
numRows = size(trackedFeatureIndx,1);
rowStart = rowEnd - numRows + 1;
if rowStart <= 1
    trackedFeatureIndxAux = [zeros(numRows,numFrames); ...
        trackedFeatureIndxAux];
    rowEnd = rowEnd + numRows;
    rowStart = rowStart + numRows;
end
trackedFeatureIndxAux(rowStart:rowEnd,:) = trackedFeatureIndx;

%remove all empty rows
trackedFeatureIndx = trackedFeatureIndxAux(rowStart:end,:);
clear trackedFeatureIndxAux

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

%clear some memory
clear costMat tmp

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system
%trackedFeatureInfo is in sparse format
trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,movieInfo,...
    numFrames,probDim);


%% %%%%% ~~ the end ~~ %%%%%

