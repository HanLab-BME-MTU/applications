function [movieInfoFinal,tracksFinal,xyDrift] = finalizeDetectionSR(...
    movieInfo,maxStd,minDuration,correctDrift)
%FINALIZEDETECTIONSR removes molecule repetition, corrects drift and removes badly localized molecules for super-resolution imaging
%
%SYNPOSIS [movieInfoFinal,tracksFinal,xyDrift] = finalizeDetectionSR(...
%    movieInfo,maxStd,minDuration,correctDrift)
%
%INPUT  movieInfo   : Output of detectSubResFeatures2D_StandAlone.
%       maxStd      : Maximum positional standard deviation allowed, in
%                     pixels. Molecules with larger positional standard
%                     deviation will be discarded. Use Inf to retain
%                     everything.
%                     Optional. Default: 0.5.
%       minDuration : Minimum duration of a molecule, in frames. Molecules
%                     that last for less frames will be discarded.
%                     Optional. Default: 1.
%       correctDrift: 1 to attempt to find fiduciary markers and correct
%                     for drift, 0 otherwise.
%                     Optional. Default: 1.
%
%OUTPUT movieInfoFinal: Like movieInfo, but after processing.
%       tracksFinal   : The surviving molecules' tracks, excluding
%                       fiduciary markers and before drift correction.
%                       These can be basically plotted on top of the raw
%                       images for visual assessment of tracking quality.
%       xyDrift       : Drift in x and y over time. First column for x,
%                       second column for y. Number of row = number of
%                       frames.
%
%REMARKS If a molecule appears more than once, its final position is taken
%as the weighted average of all its positions. The weights are derived from
%the position standard deviations. The new position and amplitude standard
%deviations are taken as the weighted standard errors. Molecules that
%appear only once keep their coordinates but get new standard deviations
%derived from molecules that appear multiple times
%
%Khuloud Jaqaman, August 2011

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--finalizeDetectionSR: Incorrect number of input arguments!');
    return
end

%ssign defaults to input arguments if not supplied
if nargin < 2 || isempty(maxStd)
    maxStd = 0.5;
end

if nargin < 3 || isempty(minDuration)
    minDuration = 1;
end

if nargin < 4 || isempty(correctDrift)
    correctDrift = 1;
end

numFrames = length(movieInfo);

%% Remove bad localizations

if maxStd < Inf
    
    fprintf('Removing localizations with standard deviation > %4.2f pixels ...\n',maxStd)
    
    %copy movieInfo into the output variable movieInfoFinal
    movieInfoFinal = movieInfo;
    
    %go over all coordinates in all frames and remove badly localized molecules
    for iFrame = 1 : numFrames
        
        if ~isempty(movieInfoFinal(iFrame).xCoord)
            
            xCoordStd = movieInfoFinal(iFrame).xCoord(:,2);
            yCoordStd = movieInfoFinal(iFrame).yCoord(:,2);
            
            indxKeep = find(xCoordStd <= maxStd & yCoordStd <= maxStd);
            
            movieInfoFinal(iFrame).xCoord = movieInfoFinal(iFrame).xCoord(indxKeep,:);
            movieInfoFinal(iFrame).yCoord = movieInfoFinal(iFrame).yCoord(indxKeep,:);
            movieInfoFinal(iFrame).amp = movieInfoFinal(iFrame).amp(indxKeep,:);
            
        end
        
    end
end


%% Get molecular tracks

disp('Tracking ...')

%general parameters
gapCloseParam.timeWindow = 1; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.
gapCloseParam.diagnostics = []; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%cost functions

%function name
costMatrices(1).funcName = 'costMatLinearMotionLink2';

%parameters
parameters.linearMotion = 0; %use linear motion Kalman filter.
parameters.minSearchRadius = 1; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = 1; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.
costMatrices(1).parameters = parameters;
clear parameters

%function name
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps2';

%parameters
parameters.linearMotion = 0; %use linear motion Kalman filter.
parameters.minSearchRadius = 1; %minimum allowed search radius.
parameters.maxSearchRadius = 1; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
parameters.brownScaling = [0 0]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
parameters.timeReachConfB = 1; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
parameters.ampRatioLimit = [0.75 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
parameters.linScaling = [0.75 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.
parameters.maxAngleVV = 25; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.
costMatrices(2).parameters = parameters;
clear parameters

%Kalman filter functions
kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

% %function name
% costMatrices(1).funcName = 'costMatStationaryLink';
% 
% %parameters
% parameters.searchRadius = 1; 
% costMatrices(1).parameters = parameters;
% clear parameters
% 
% %function name
% costMatrices(2).funcName = 'costMatStationaryCloseGaps';
% 
% %parameters
% parameters.searchRadius = 1;
% parameters.gapPenalty = 1.5;
% costMatrices(2).parameters = parameters;
% clear parameters
% 
% %Kalman filter functions
% kalmanFunctions = [];

%saveResults
saveResults = 0; %don't save results

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%construct tracks of detected molecules
tracksFinal = trackCloseGapsKalmanSparse(movieInfoFinal,costMatrices,...
    gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

parameterSet.searchRadius = 1;
parameterSet.timeWindow = 1;
parameterSet.gapPenalty = 1.5;

% tracksFinal = trackCloseGapsSR(movieInfoFinal,parameterSet,...
%     probDim,saveResults,verbose);

%convert tracks into matrix format
tracksFinal = convStruct2MatNoMS(tracksFinal);

%get track start, end and lifetime information
trackSEL = getTrackSEL(tracksFinal);

%get number of tracks
numTracks = size(trackSEL,1);

%% Remove tracks that are too short

if minDuration > 1
    
    fprintf('Removing tracks lasting < %i frames ...\n',minDuration)
    
    %remove tracks lasting less than minDuration
    indxKeep = find(trackSEL(:,3) >= minDuration);
    tracksFinal = tracksFinal(indxKeep,:);
    trackSEL = trackSEL(indxKeep,:);
    
    %get number of tracks
    numTracks = length(indxKeep);
    
end

%% Correct for drift

if correctDrift
    
    disp('Correcting for drift ...')
    
    %find objects that last for the whole movie, or say 90% of the movie
    %those can only be fiduciary markers
    indxFiduciary = find(trackSEL(:,3) >= 0.9*numFrames);
    
    if ~isempty(indxFiduciary)
        
        fprintf('Found %i fiduciary markers \n',length(indxFiduciary))
        
        %extract the fiduciary marker tracks
        tracksFiduciary = tracksFinal(indxFiduciary,:);
        
        %calculate drift in x and y
        xCoordFiduciary = tracksFiduciary(:,1:8:end);
        yCoordFiduciary = tracksFiduciary(:,2:8:end);
        xCoordDrift = mean(xCoordFiduciary,1);
        xCoordDrift = xCoordDrift - xCoordDrift(1);
        yCoordDrift = mean(yCoordFiduciary,1);
        yCoordDrift = yCoordDrift - yCoordDrift(1);
        
        %remove tracks of fiduciary markers from list of tracks
        indxKeep = setdiff(1:numTracks,indxFiduciary);
        tracksFinal = tracksFinal(indxKeep,:);
        trackSEL = trackSEL(indxKeep,:);
        numTracks = length(indxKeep);
        
        %eliminate drift from tracks
        tracksFinal(:,1:8:end) = tracksFinal(:,1:8:end) - ...
            repmat(xCoordDrift,numTracks,1);
        tracksFinal(:,2:8:end) = tracksFinal(:,2:8:end) - ...
            repmat(yCoordDrift,numTracks,1);
                
    else
        
        disp('Could not find any fiduciary markers');
        xCoordDrift = zeros(1,numFrames);
        yCoordDrift = zeros(1,numFrames);
        
        
    end
    
else
    
    xCoordDrift = zeros(1,numFrames);
    yCoordDrift = zeros(1,numFrames);
    
end

%for output
xyDrift = [xCoordDrift' yCoordDrift'];

%% Eliminate repetition and calculate final position of each molecule

disp('Eliminating repetition ...')

%initilize matrix storing molecule information
%information stored is: frame of first appearance, number of appearances,
%x-coordinate, y-coordinate, amplitude, x-std, y-std, amplitude-std
moleculeInfo = zeros(numTracks,8);

%get indices of molecules that appear for more than one frame
indxMultiple = find(trackSEL(:,3) > 1);

%go over these molecules
for iMolecule = indxMultiple'
    
    %get appearance and disappearance frames
    moleculeStart = trackSEL(iMolecule,1);
    moleculeEnd   = trackSEL(iMolecule,2);
    
    %extract molecule's coordinates and amplitude and their standard deviations
    coordAmpStd = tracksFinal(iMolecule,(moleculeStart-1)*8+1:moleculeEnd*8)';
    coordAmpStd = reshape(coordAmpStd,8,[])';
    coordAmpStd = coordAmpStd(:,[1 2 4 5 6 8]);
    
    %find and keep only non-NaN rows
    indxNoNaN = find(~isnan(coordAmpStd(:,1)));
    coordAmpStd = coordAmpStd(indxNoNaN,:);
    
    %determine number of times molecule appears
    %it can be different from track lifetime because of gap closing
    numObs = length(indxNoNaN);
    
    %calculate molecule's average standard deviation over x and y
    xyStdAve = mean(coordAmpStd(:,4:5),2);
    
    %calculate each observation's weight, normalized such that the
    %sum of weights = 1
    obsWeight = 1 ./ (xyStdAve.^2);
    obsWeight = obsWeight / sum(obsWeight);
    
    %calculate weighted average coordinates and intensity
    weightedAve = sum( repmat(obsWeight,1,3).*coordAmpStd(:,1:3),1 );
    
    %calculate weighted standard deviations
    sumSqDiff = sum( repmat(obsWeight,1,3) .* ...
        (coordAmpStd(:,1:3)-repmat(weightedAve,numObs,1)).^2,1 );
    denominator = (numObs-1) / numObs;
    weightedStd = sqrt( sumSqDiff / denominator );
    
    %collect weighted average and standard deviation
    coordAmpStd = [weightedAve weightedStd];
    
    %put molecule's information in the molecule information matrix
    %put it in the frame where it first appears
    moleculeInfo(iMolecule,:) = [moleculeStart numObs coordAmpStd];
    
end

%calculate the average standard deviation from all molecules that appear
%more than once in order to assign it to molecules that appear only once
%calculate a weighted average in order to account for the different
%number of appearances for the different molecules
stdMultiple = moleculeInfo(indxMultiple,6:8);
obsWeight = moleculeInfo(indxMultiple,2);
obsWeight = obsWeight / sum(obsWeight);
stdMultiple = sum( repmat(obsWeight,1,3).*stdMultiple,1 );

%get indices of molecules that appear for only one frame
indxOne = find(trackSEL(:,3) == 1);

%give them all the average standard deviations determined above
moleculeInfo(indxOne,6:8) = repmat(stdMultiple,length(indxOne),1);

%go over all molecules that appear for only one frame
for iMolecule = indxOne'
    
    %get appearance frame
    moleculeStart = trackSEL(iMolecule,1);
    
    %extract molecule's coordinates and amplitude and their standard deviations
    coordAmpStd = tracksFinal(iMolecule,(moleculeStart-1)*8+1:moleculeStart*8);
    
    %put molecule's information in the molecule information matrix
    moleculeInfo(iMolecule,1:5) = [moleculeStart 1 coordAmpStd([1 2 4])];
    
end

%for all molecules, divide standard deviation by sqrt(n) in order to
%calculate standard error of the mean
moleculeInfo(:,6:8) = moleculeInfo(:,6:8) ./ repmat(sqrt(moleculeInfo(:,2)),1,3);

%go over all frames ...
for iFrame = 1 : numFrames
    
   %get molecule information for this frame
   moleculeInfoFrame = moleculeInfo(moleculeInfo(:,1)==iFrame,:);
   
   %keep only coordinates belonging to these molecules
   movieInfoFinal(iFrame).xCoord = moleculeInfoFrame(:,[3 6]);
   movieInfoFinal(iFrame).yCoord = moleculeInfoFrame(:,[4 7]);
   movieInfoFinal(iFrame).amp    = moleculeInfoFrame(:,[5 8]);
   movieInfoFinal(iFrame).numAppearance = moleculeInfoFrame(:,2);
   
end

disp('DONE')

%% ~~~ end ~~~

