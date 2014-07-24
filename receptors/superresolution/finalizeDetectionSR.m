function [movieInfoFinal,tracksFinal,xyDrift] = finalizeDetectionSR(...
    movieInfo,param,saveResults)
%FINALIZEDETECTIONSR removes molecule repetition, corrects drift and removes badly localized molecules for super-resolution imaging
%
%SYNPOSIS [movieInfoFinal,tracksFinal,xyDrift] = finalizeDetectionSR(...
%    movieInfo,param,saveResults)
%
%INPUT  movieInfo   : Output of detectSubResFeatures2D_StandAlone.
%       param       : Structure with fields:
%           .maxStd    : Maximum positional standard deviation allowed, in
%                        pixels. Molecules with larger positional standard
%                        deviation will be discarded. Use Inf to retain
%                        everything.
%                        Optional. Default: 0.5.
%           .minDuration : Minimum duration of a molecule, in frames.
%                        Molecules that last for less frames will be
%                        discarded.
%                        Optional. Default: 1.
%           .correctDrift: 1 to attempt to find fiduciary markers and
%                        correct for drift, 0 otherwise.
%                        Optional. Default: 1.
%           .tracking  : Structure listing tracking parameteres:
%               .searchRadius: Search radius for linking molecules between
%                              frames and for gap closing.
%                              Optional. Default: 1.
%               .timeWindow  : Time window for gap closing. A time window
%                              of n allows molecule disappearance for at
%                              most n-1 consecutive frames.
%                              Optional. Default: 3.
%               .gapPenalty  : Penalty for increasing gap length in gap
%                              closing. For a disappearance of n frames,
%                              penalty is defined as gapPenalty^n.
%                              Optional. Default: 1.5.
%                     Whole structure optional. 
%       saveResults : 0 if no saving is requested.
%                     If saving is requested, structure with fields:
%           .dir       : Directory where results should be saved.
%                        Optional. Default: current directory.
%           .filename  : Name of file where results should be saved.
%                        Optional. Default: detectedFeaturesFinal.
%                     Whole structure optional.
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

%check parameters structure
if nargin < 2 || isempty(param)
    param.maxStd = 0.5;
    param.minDuration = 1;
    param.correctDrift = 1;
    param.tracking.searchRadius = 1;
    param.tracking.timeWindow = 3;
    param.tracking.gapPenalty = 1.5;
else
    if ~isfield(param,'maxStd') || isempty(param.maxStd)
        param.maxStd = 0.5;
    end
    if ~isfield(param,'minDuration') || isempty(param.minDuration)
        param.minDuration = 1;
    end
    if ~isfield(param,'correctDrift') || isempty(param.correctDrift)
        param.correctDrift = 1;
    end
    if ~isfield(param,'tracking') || isempty(param.tracking)
        param.tracking.searchRadius = 1;
        param.tracking.timeWindow = 3;
        param.tracking.gapPenalty = 1.5;
    else
        trackingParam = param.tracking;
        if ~isfield(trackingParam,'searchRadius')
            param.tracking.searchRadius = 1;
        end
        if ~isfield(trackingParam,'timeWindow')
            param.tracking.timeWindow = 3;
        end
        if ~isfield(trackingParam,'gapPenalty')
            param.tracking.gapPenalty = 1.5;
        end
    end
end
maxStd = param.maxStd;
minDuration = param.minDuration;
correctDrift = param.correctDrift;
trackingParam = param.tracking;

%check whether to save results and where
if nargin < 3 || isempty(saveResults)
    saveResDir = pwd;
    saveResFile = 'detectedFeaturesFinal.mat';
    saveResults.dir = pwd;
else
    if isstruct(saveResults)
        if ~isfield(saveResults,'dir') || isempty(saveResults.dir)
            saveResDir = pwd;
        else
            saveResDir = saveResults.dir;
        end
        if ~isfield(saveResults,'filename') || isempty(saveResults.filename)
            saveResFile = 'detectedFeatures.mat';
        else
            saveResFile = saveResults.filename;
        end
    else
        saveResults = 0;
    end
end

%get number of frames in movie
numFrames = length(movieInfo);

%copy movieInfo into the output variable movieInfoFinal
movieInfoFinal = movieInfo;

%% Remove bad localizations

if maxStd < Inf
    
    fprintf('Removing localizations with standard deviation > %4.2f pixels ...\n',maxStd)
    
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

tracksFinal = trackCloseGapsSR(movieInfoFinal,trackingParam,2,0,1);

%convert tracks into matrix format
%use sparse matrix for the sake of memory
tracksFinal = convStruct2SparseMatNoMS(tracksFinal);

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
        xCoordFiduciary = full(tracksFiduciary(:,1:8:end));
        yCoordFiduciary = full(tracksFiduciary(:,2:8:end));
        xCoordFiduciary(xCoordFiduciary==0) = NaN;
        yCoordFiduciary(yCoordFiduciary==0) = NaN;
        xCoordDrift = nanmean(xCoordFiduciary,1);
        yCoordDrift = nanmean(yCoordFiduciary,1);
        
        %find first non-NaN value in xCoordDrift and yCoordDrift
        indxFirst = find(~isnan(xCoordDrift),1,'first');
        if indxFirst > 1
            fprintf('Can only estimate drift from %i onwards \n',indxFirst)
        end
        xCoordDrift = xCoordDrift - xCoordDrift(indxFirst);
        yCoordDrift = yCoordDrift - yCoordDrift(indxFirst);
        
        %fill in values for frames without drift estimation
        xCoordDrift(1:indxFirst-1) = 0;
        yCoordDrift(1:indxFirst-1) = 0;
        indxNaN = find(isnan(xCoordDrift));
        for i=indxNaN
            xCoordDrift(i) = xCoordDrift(i-1);
            yCoordDrift(i) = yCoordDrift(i-1);
        end
        
        %remove tracks of fiduciary markers from list of tracks
        indxKeep = setdiff(1:numTracks,indxFiduciary);
        tracksFinal = tracksFinal(indxKeep,:);
        trackSEL = trackSEL(indxKeep,:);
        numTracks = length(indxKeep);
        
        %eliminate drift from tracks
        for iTrack = 1 : numTracks
            xCoordTrack = tracksFinal(iTrack,1:8:end);
            xCoordTrack(trackSEL(iTrack,1):trackSEL(iTrack,2)) = ...
                xCoordTrack(trackSEL(iTrack,1):trackSEL(iTrack,2)) - ...
                xCoordDrift(trackSEL(iTrack,1):trackSEL(iTrack,2));
            tracksFinal(iTrack,1:8:end) = xCoordTrack;
            yCoordTrack = tracksFinal(iTrack,2:8:end);
            yCoordTrack(trackSEL(iTrack,1):trackSEL(iTrack,2)) = ...
                yCoordTrack(trackSEL(iTrack,1):trackSEL(iTrack,2)) - ...
                yCoordDrift(trackSEL(iTrack,1):trackSEL(iTrack,2));
            tracksFinal(iTrack,2:8:end) = yCoordTrack;
        end
                
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
    sumWeights = sum(obsWeight);
    obsWeight = obsWeight / sumWeights;
    
    %calculate weighted average of coordinates
    weightedAve = sum( repmat(obsWeight,1,2).*coordAmpStd(:,1:2),1 );
    
    %sum up the amplitudes
    weightedAve(:,3) = sum(coordAmpStd(:,3),1);
    
    %calculate weighted standard deviation of coordinates
    sumSqDiff = sum( repmat(obsWeight,1,2) .* ...
        (coordAmpStd(:,1:2)-repmat(weightedAve(:,1:2),numObs,1)).^2,1 );
    denominator = (numObs-1) / numObs;
    weightedStd = sqrt( sumSqDiff / denominator );
    
    %calculate amplitude standard deviation
    %note the multiplication by sqrt(numObs) - this is because at this stage
    %we want the sample standard deviation and not the standard error of the
    %mean, which is calculated from the standard deviation at a later step
    weightedStd(:,3) = sqrt(sum(coordAmpStd(:,6).^2)) * sqrt(numObs);
    
    %collect new coordinates and amplitude and their standard deviations
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

%% Save results
if isstruct(saveResults)
    save([saveResDir filesep saveResFile],'movieInfoFinal','param',...
        'tracksFinal','xyDrift');
end

disp('DONE')

%% ~~~ end ~~~

