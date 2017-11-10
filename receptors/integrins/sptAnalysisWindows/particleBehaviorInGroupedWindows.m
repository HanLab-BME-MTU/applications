function [sptPropInWindow,windowTrackAssign,trackWindowAssign,windowSize,analysisParam] = ...
    particleBehaviorInGroupedWindows(tracksFinal,winPositions,winFrames,...
    protSamples,diffAnalysisRes,minLength,bandRange,windowRange,frameRange)
%PARTICLEBEHAVIORINGROUPEDWINDOWS averages single particle behavior in windows grouped based on protrusion activity
%
%SYNOPSIS[sptPropInWindow,windowTrackAssign,trackWindowAssign,windowSize,analysisParam] = ...
%    particleBehaviorInGroupedWindows(tracksFinal,winPositions,winFrames,...
%    protSamples,diffAnalysisRes,minLength,bandRange,windowRange,frameRange)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       winPositions   : A 2D array of the window edges.
%                        Number of rows = number of window frames.
%                        Number of columns = number of windows parallel to
%                        Each entry is the output of Hunter's new
%                        windowing software.
%                        Basically, to make this variable, one puts
%                        together the windows of each frame coming out of
%                        the windowing software.
%       winFrames      : The frames at which there are windows.
%       protSamples    : The protrusion samples as output by the windowing
%                        software.
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%                        Optional. If not input but needed, it will be
%                        calculated within the code.
%       minLength      : Minimum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: 5.
%       bandRange      : Nx2 array indicating range of bands (i.e. slabs
%                        going away from the cell edge) to include in
%                        analysis. Each row indicates bands that will be
%                        grouped together.
%                        Optional. Default: [1 2; 3 4; 5 6; 7 8; 9 10].
%       windowRange    : Vector with 2 entries indicating range of windows
%                        (i.e. window number along cell edge) to include in
%                        analysis.
%                        Optional. Default: [1 (last window)].
%       frameRange     : Vector with 2 entries indicating range of window
%                        frames to include in analysis.
%                        Optional. Default: [1 (last frame)].
%
%OUTPUT sptPropInWindow: Structure storing various single particle
%                        properties, with the following fields:
%               From the tracks directly ...
%           .spDensity      : Single particle density.
%           .f2fDisp        : Average frame-to-frame displacement.
%           .angleMean      : Average angle between direction of motion and
%                             protrusion vector.
%           .angleStd       : STD of angle between direction of motion and
%                             protrusion vector.
%           .dirDisp        : Average frame-to-frame displacement along
%                             direction of motion.
%               From track diffusion analysis ...
%           .fracUnclass    : Fraction of completely unclassified tracks
%                             (i.e. tracks < frames).
%           .fracLin        : Among tracks classifiable by asymmetry
%                             analysis (i.e. tracks >= 5 frames),
%                             fraction of tracks classified as linear.
%           .fracIso        : //, fraction of tracks classified as
%                             isotropic (i.e. 1 - fracLin).
%           .fracIsoUnclass : Among isotropic tracks, fraction of tracks
%                             with unclassified diffusion (i.e. tracks >= 5
%                             frames but < 20 frames).
%           .fracConf       : Among isotropic tracks that are also
%                             classifiable by diffusion analysis (i.e.
%                             tracks >= 20 frames), fraction of confined.
%           .fracBrown      : //, fraction of Brownian.
%           .fracDir        : //, fraction of directed.
%           .diffCoef       : Average diffusion coefficient.
%           .confRad        : Average confinement radius.
%                        Each field contains 2 sub-fields, "values" and
%                        "numPoints" storing the values of the properties
%                        and the number of points contributing to each
%                        value, respectively.
%                        Each of these fields is a 3-D matrix, of
%                        dimensions (number of bands) x (number of windows)
%                        x (number of window frames-1).
%       windowTrackAssign, trackWindowAssign: Output of assignTracks2Windows.
%       windowSize     : 3-D matrix of dimensions (number of bands) x (number
%                        of slices) x (number of window frames - 1)
%                        storing the size of each window. NaN indicates a
%                        window of zero size.
%
%REMARKS This code is designed for experiments where the particle
%        trajectories are sampled much more frequently than the cell edge.
%        It assumes that particle lifetimes are much shorter than the time
%        between cell edge frames.
%
%        For a different scenario where particle lifetimes are longer than
%        the time between cell edge frames, the tracks should not be
%        grouped like this. Rather, each track should get divided into
%        several segments corresponding to the times between cell edge
%        frames and each track segment should be analyzed separately.
%        Something like that.
%
%Khuloud Jaqaman, January 2012

%% Input

if nargin < 4
    disp('--particleBehaviorInGroupedWindows: Missing input arguments!');
    return
end

if nargin < 5 || isempty(diffAnalysisRes)
    diffAnalysisRes = trackDiffusionAnalysis1(tracksFinal,1,2,0,0.05);
end

if nargin < 6 || isempty(minLength)
    minLength = 5;
end

%get number of frames that have windows and number of windows parallel to
%the edge
[numWinFrames,numWinPara] = size(winPositions);

%find number of windows perpendicular to the edge
nBands = cellfun(@(x)(numel(x)),winPositions);
numWinPerp = max(nBands(:));

%determine the number of SPT frames between window frames
numSPTFrames = winFrames(2) - winFrames(1);

if nargin < 7 || isempty(bandRange)
    bandRange1 = (1:2:9)';
    bandRange2 = (2:2:min(10,numWinPerp))';
    bandRange1 = bandRange1(1:length(bandRange2));
    bandRange = [bandRange1 bandRange2];
end
numBandRange = size(bandRange,1);

if nargin < 8 || isempty(windowRange)
    windowRange = [1 numWinPara];
end

if nargin < 9 || isempty(frameRange)
    frameRange = [1 numWinFrames];
end

%% Trajectory pre-processing

%keep only trajectories longer than minLength
criteria.lifeTime.min = minLength;
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx,:);
diffAnalysisRes = diffAnalysisRes(indx);

%get number of trajectories
numTracks = length(indx);

%divide the trajectories among the windows
[windowTrackAssign,trackWindowAssign] = assignTracks2Windows(tracksFinal,winPositions,winFrames,1);

%% Window pre-processing

%initialize array storing window sizes
winSize = NaN(numWinPerp,numWinPara,numWinFrames-1);

%go over all windows and get their sizes
for iFrame = 1 : numWinFrames-1
    for iPara = 1 : numWinPara
        for iPerp = 1 : nBands(iFrame,iPara)
            
            %if this window has a finite size
            if ~isempty(winPositions{iFrame,iPara}{iPerp})
                
                %get the window boundaries
                windowsPoly = [winPositions{iFrame,iPara}{iPerp}{:}];
                winX = windowsPoly(1,:);
                winY = windowsPoly(2,:);
                
                %calculate the window size
                winSize(iPerp,iPara,iFrame) = polyarea(winX,winY);
                
            end
            
        end
    end
end
%make sure that there are no zeros
winSize(winSize==0) = NaN;

%copy winSize into the output variable windowSize and convert NaNs to zeros
%in windowSize
windowSize = winSize;
windowSize(isnan(windowSize)) = 0;

%get protrusion vectors
protVec = protSamples.avgVector;

%for windows/frames outside of the range of interest, convert their values
%to NaN
protVec(1:windowRange(1)-1,:,:) = NaN;
protVec(windowRange(2)+1:end,:,:) = NaN;
protVec(:,1:frameRange(1)-1,:) = NaN;
protVec(:,frameRange(2):end,:) = NaN;

%calculate normalized protrusion vectors
protVecMag = sqrt(sum(protVec.^2,3));
protVecUnit = protVec ./ repmat(protVecMag,[1 1 2]);

%get the activity classification of each window in each frame
windowMotionType = classifyEdgeMotion(protSamples,0,windowRange,frameRange);

%% Particle behavior pre-processing

%get the start, end and lifetime of each track segment
trackLft = getTrackSEL(tracksFinal,1);
% trackSE = trackLft(:,1:2);
trackLft = trackLft(:,3);

%get the number of track segments
numSegments = length(trackLft);

%From asymmetry and diffusion analysis ...

%get trajectory classifications
trajClassDiff = vertcat(diffAnalysisRes.classification);
trajClassAsym = trajClassDiff(:,1); %asymmetry classification
trajClassDiff = trajClassDiff(:,2); %diffusion classification

%WARNING
%MAKE tracjClassAsym ALL NaN, SO THAT THERE IS NO ASYMMETRY CLASSIFICATION
trajClassAsym(:) = NaN;

%process classifications to make categories
if all(isnan(trajClassAsym))
    trajClass = trajClassDiff;
else
    trajClass = NaN(numSegments,1);
    trajClass(trajClassAsym==0&trajClassDiff==1) = 1; %isotropic+confined
    trajClass(trajClassAsym==0&trajClassDiff==2) = 2; %isotropic+free
    trajClass(trajClassAsym==0&trajClassDiff==3) = 3; %isotropic+directed
    trajClass(trajClassAsym==0&isnan(trajClassDiff)) = 4; %isotropic+undetermined
    trajClass(trajClassAsym==1) = 5; %linear+anything
end

%get diffusion coefficients
diffCoefGen = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');

%get confinement radii
confRadAll = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');

%From tracks directly ...

%get direction of motion, angle with protrusion vector and various
%frame-to-frame displacements
[~,angleWithProtTmp,f2fDispTmp,paraDirDispTmp,perpDirDispTmp,...
    paraProtDispTmp,perpProtDispTmp,asymParamTmp] = ...
    trackMotionCharProtrusion(tracksFinal,protVecUnit,trackWindowAssign,minLength);

%calculate ratio of parallel to perpendicular displacements
ratioDispDirTmp = abs( paraDirDispTmp ./ perpDirDispTmp );
ratioDispProtTmp = abs( paraProtDispTmp ./ perpProtDispTmp );

%% Calculate property values per window group

%initialize output variable
sptPropInWindow = repmat(struct('directAll',[],'directPos',[],...
    'directNeg',[],'diffAnalysis',[]),numBandRange,1);

%go over the different band ranges
for iBandRange = 1 : numBandRange
    
    %get current bands to analyze
    bandsCurrent = bandRange(iBandRange,1) : bandRange(iBandRange,2);
    
    %initialize variables storing results
    [...
        f2fDispMag2DAll,angleProtAll,asymParamAll,...
        f2fDispMagParaDirAll,f2fDispMagPerpDirAll,f2fDispSignParaDirAll,f2fDispSignPerpDirAll,...
        f2fDispMagParaProtAll,f2fDispMagPerpProtAll,f2fDispSignParaProtAll,f2fDispSignPerpProtAll,...
        ratioDispMagDirAll,ratioDispSignDirAll,ratioDispMagProtAll,ratioDispSignProtAll,...
        ...
        f2fDispMag2DPos,angleProtPos,asymParamPos,...
        f2fDispMagParaDirPos,f2fDispMagPerpDirPos,f2fDispSignParaDirPos,f2fDispSignPerpDirPos,...
        f2fDispMagParaProtPos,f2fDispMagPerpProtPos,f2fDispSignParaProtPos,f2fDispSignPerpProtPos,...
        ratioDispMagDirPos,ratioDispSignDirPos,ratioDispMagProtPos,ratioDispSignProtPos,...
        ...
        f2fDispMag2DNeg,angleProtNeg,asymParamNeg,...
        f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,...
        f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,...
        ratioDispMagDirNeg,ratioDispSignDirNeg,ratioDispMagProtNeg,ratioDispSignProtNeg,...
        ...
        spDensity,fracNetDispNeg,...
        fracUnclass,fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
        diffCoefAll,diffCoefConf,diffCoefBrown,confRad] = ...
        ...
        deal(struct('mean',NaN(12,1),'std',NaN(12,1),'numPoints',zeros(12,1)));
    
    classType = [...
        1 2 2; ...
        1 2 3; ...
        1 3 2; ...
        1 3 3; ...
        2 1 1; ...
        2 1 3; ...
        2 3 1; ...
        2 3 3; ...
        3 1 1; ...
        3 1 2; ...
        3 2 1; ...
        3 2 2];
    
    for iClass = 1 : 12
        
        %find all windows with the current classification:
        [indxWindow,indxFrame] = find( ...
            windowMotionType(:,:,1)==classType(iClass,1) & ...
            windowMotionType(:,:,2)==classType(iClass,2) & ...
            windowMotionType(:,:,3)==classType(iClass,3) & ...
            windowMotionType(:,:,4)>=2 );
        
        %get number of windows in group
        numWinInGroup = size(indxWindow,1);
        
        %get the tracks that belong to this window group
        bandVec = repmat(bandsCurrent,numWinInGroup,1);
        bandVec = bandVec(:);
        winVec = repmat(indxWindow,length(bandsCurrent),1);
        frameVec = repmat(indxFrame,length(bandsCurrent),1);
        linearInd = sub2ind([numWinPerp numWinPara numWinFrames-1],bandVec,winVec,frameVec);
        tracksCurrent = vertcat(windowTrackAssign{linearInd});
        numTracksCurrent = length(tracksCurrent);
        
        %calculate the sum of all window sizes in order to calculate
        %particle density later
        windowGroupSize = nansum(winSize(linearInd));
        
        if numTracksCurrent == 0 %if there are no tracks in this window ...
            
            %calculate particle density
            %note that here particle density is not automatically 0
            %because it can be that the window size is 0 (stored as
            %NaN), in which case the density will be NaN
            spDensity.mean(iClass) = 0 / windowGroupSize;
            spDensity.numPoints(iClass) = numTracksCurrent;
            
        else %if there are tracks in this window ...
            
            %calculate particle density
            spDensity.mean(iClass) = sum(trackLft(tracksCurrent)) / ...
                windowGroupSize / numSPTFrames;
            spDensity.numPoints(iClass) = numTracksCurrent;
            
            %From the tracks directly ...
            
            %characteristics of all current tracks together
            [angleProtAll,f2fDispMag2DAll,...
                f2fDispSignParaDirAll,f2fDispSignPerpDirAll,...
                f2fDispMagParaDirAll,f2fDispMagPerpDirAll,...
                f2fDispSignParaProtAll,f2fDispSignPerpProtAll,...
                f2fDispMagParaProtAll,f2fDispMagPerpProtAll,...
                ratioDispSignDirAll,ratioDispMagDirAll,...
                ratioDispSignProtAll,ratioDispMagProtAll,...
                asymParamAll] = ...
                putTrackCharTogether(tracksCurrent,...
                numTracksCurrent,angleWithProtTmp,f2fDispTmp,...
                paraDirDispTmp,perpDirDispTmp,...
                paraProtDispTmp,perpProtDispTmp,...
                ratioDispDirTmp,ratioDispProtTmp,...
                asymParamTmp,...
                angleProtAll,f2fDispMag2DAll,...
                f2fDispSignParaDirAll,f2fDispSignPerpDirAll,...
                f2fDispMagParaDirAll,f2fDispMagPerpDirAll,...
                f2fDispSignParaProtAll,f2fDispSignPerpProtAll,...
                f2fDispMagParaProtAll,f2fDispMagPerpProtAll,...
                ratioDispSignDirAll,ratioDispMagDirAll,...
                ratioDispSignProtAll,ratioDispMagProtAll,asymParamAll,iClass);
            
            %divide current tracks into those with net displacement
            %parallel or opposite to protrusion vector
            paraProtDispCurrent = paraProtDispTmp(tracksCurrent,1);
            tracksPos = tracksCurrent(paraProtDispCurrent >= 0);
            tracksNeg = tracksCurrent(paraProtDispCurrent < 0);
            numTracksPos = length(tracksPos);
            numTracksNeg = length(tracksNeg);
            
            %store fraction of tracks with net displacement opposite to
            %protrusion vector (i.e. negative net displacement)
            fracNetDispNeg.mean(iClass) = numTracksNeg / numTracksCurrent;
            fracNetDispNeg.numPoints(iClass) = numTracksCurrent;
            
            %characteristics of tracks with net positive displacement
            [angleProtPos,f2fDispMag2DPos,...
                f2fDispSignParaDirPos,f2fDispSignPerpDirPos,...
                f2fDispMagParaDirPos,f2fDispMagPerpDirPos,...
                f2fDispSignParaProtPos,f2fDispSignPerpProtPos,...
                f2fDispMagParaProtPos,f2fDispMagPerpProtPos,...
                ratioDispSignDirPos,ratioDispMagDirPos,...
                ratioDispSignProtPos,ratioDispMagProtPos,...
                asymParamPos] = ...
                putTrackCharTogether(tracksPos,...
                numTracksPos,angleWithProtTmp,f2fDispTmp,...
                paraDirDispTmp,perpDirDispTmp,...
                paraProtDispTmp,perpProtDispTmp,...
                ratioDispDirTmp,ratioDispProtTmp,...
                asymParamTmp,...
                angleProtPos,f2fDispMag2DPos,...
                f2fDispSignParaDirPos,f2fDispSignPerpDirPos,...
                f2fDispMagParaDirPos,f2fDispMagPerpDirPos,...
                f2fDispSignParaProtPos,f2fDispSignPerpProtPos,...
                f2fDispMagParaProtPos,f2fDispMagPerpProtPos,...
                ratioDispSignDirPos,ratioDispMagDirPos,...
                ratioDispSignProtPos,ratioDispMagProtPos,asymParamPos,iClass);
            
            %characteristics of tracks with net negative displacement
            [angleProtNeg,f2fDispMag2DNeg,...
                f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,...
                f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,...
                f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,...
                f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,...
                ratioDispSignDirNeg,ratioDispMagDirNeg,...
                ratioDispSignProtNeg,ratioDispMagProtNeg,...
                asymParamNeg] = ...
                putTrackCharTogether(tracksNeg,...
                numTracksNeg,angleWithProtTmp,f2fDispTmp,...
                paraDirDispTmp,perpDirDispTmp,...
                paraProtDispTmp,perpProtDispTmp,...
                ratioDispDirTmp,ratioDispProtTmp,...
                asymParamTmp,...
                angleProtNeg,f2fDispMag2DNeg,...
                f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,...
                f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,...
                f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,...
                f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,...
                ratioDispSignDirNeg,ratioDispMagDirNeg,...
                ratioDispSignProtNeg,ratioDispMagProtNeg,asymParamNeg,iClass);
            
            %From the asymmetry and diffusion analysis ...
            
            %calculate the fraction of tracks in each motion category
            
            %first completely unclassified tracks
            trajClassCurrent = trajClass(tracksCurrent);
            fracUnclass.mean(iClass) = ...
                length(find(isnan(trajClassCurrent)))/numTracksCurrent;
            fracUnclass.numPoints(iClass) = numTracksCurrent;
            
            %then tracks with asymmetry analysis
            trajClassCurrent = trajClassCurrent(~isnan(trajClassCurrent));
            numTracksAsym = length(trajClassCurrent);
            if numTracksAsym > 0
                
                %tracks classified as linear or isotropic
                fracLin.mean(iClass) = length(find(trajClassCurrent==5))/numTracksAsym;
                fracLin.numPoints(iClass) = numTracksAsym;
                fracIso.mean(iClass) = length(find(trajClassCurrent~=5))/numTracksAsym;
                fracIso.numPoints(iClass) = numTracksAsym;
                
                %within the tracks classified as isotropic
                trajClassCurrent = trajClassCurrent(trajClassCurrent~=5);
                numTracksIso = length(trajClassCurrent);
                if numTracksIso > 0
                    
                    %tracks without diffusion analysis
                    fracIsoUnclass.mean(iClass) = length(find(trajClassCurrent==4))/numTracksIso;
                    fracIsoUnclass.numPoints(iClass) = numTracksIso;
                    
                    %finally tracks with diffusion analysis
                    trajClassCurrent = trajClassCurrent(trajClassCurrent~=4);
                    numTracksDiff = length(trajClassCurrent);
                    if numTracksDiff > 0
                        
                        fracClass = hist(trajClassCurrent,(1:3))/numTracksDiff;
                        fracConf.mean(iClass) = fracClass(1);
                        fracConf.numPoints(iClass) = numTracksDiff;
                        fracBrown.mean(iClass) = fracClass(2);
                        fracBrown.numPoints(iClass) = numTracksDiff;
                        fracDir.mean(iClass) = fracClass(3);
                        fracDir.numPoints(iClass) = numTracksDiff;
                        
                    end %(if numTracksDiff > 0)
                    
                end %(if numTracksIso > 0)
                
            end %(if numTracksAsym > 0)
            
            %calculate the average diffusion coefficient
            %for everything together, for only confined, and for only
            %Brownian
            tmpVec = diffCoefGen(tracksCurrent);
            tmpVecAll = tmpVec(~isnan(tmpVec));
            tmpVecConf = tmpVec(trajClass(tracksCurrent) == 1);
            tmpVecBrown = tmpVec(trajClass(tracksCurrent) == 2);
            if ~isempty(tmpVecAll)
                diffCoefAll.mean(iClass) = mean(tmpVecAll);
                diffCoefAll.std(iClass) = std(tmpVecAll);
                diffCoefAll.numPoints(iClass) = length(tmpVecAll);
            end
            if ~isempty(tmpVecConf)
                diffCoefConf.mean(iClass) = mean(tmpVecConf);
                diffCoefConf.std(iClass) = std(tmpVecConf);
                diffCoefConf.numPoints(iClass) = length(tmpVecConf);
            end
            if ~isempty(tmpVecBrown)
                diffCoefBrown.mean(iClass) = mean(tmpVecBrown);
                diffCoefBrown.std(iClass) = std(tmpVecBrown);
                diffCoefBrown.numPoints(iClass) = length(tmpVecBrown);
            end
            
            %calculate the average confinement radius for confined
            %particles
            tmpVec = confRadAll(tracksCurrent);
            tmpVec = tmpVec(trajClass(tracksCurrent) == 1);
            if ~isempty(tmpVec)
                confRad.mean(iClass) = mean(tmpVec);
                confRad.std(iClass) = std(tmpVec);
                confRad.numPoints(iClass) = length(tmpVec);
            end
            
        end %(if numTracksCurrent == 0 ... else ...)
        
    end %(for iClass = 1 : 12)
    
    %store all properties in output structure
    
    directAll = struct('spDensity',spDensity,'fracNetDispNeg',fracNetDispNeg,...
        'f2fDispMag2D',f2fDispMag2DAll,'asymParam',asymParamAll,'angleProt',angleProtAll,...
        'f2fDispMagParaDir',f2fDispMagParaDirAll,'f2fDispMagPerpDir',f2fDispMagPerpDirAll,...
        'f2fDispSignParaDir',f2fDispSignParaDirAll,'f2fDispSignPerpDir',f2fDispSignPerpDirAll,...
        'f2fDispMagParaProt',f2fDispMagParaProtAll,'f2fDispMagPerpProt',f2fDispMagPerpProtAll,...
        'f2fDispSignParaProt',f2fDispSignParaProtAll,'f2fDispSignPerpProt',f2fDispSignPerpProtAll,...
        'ratioDispMagDir',ratioDispMagDirAll,'ratioDispSignDir',ratioDispSignDirAll,...
        'ratioDispMagProt',ratioDispMagProtAll,'ratioDispSignProt',ratioDispSignProtAll);
    
    directPos = struct('f2fDispMag2D',f2fDispMag2DPos,'asymParam',asymParamPos,'angleProt',angleProtPos,...
        'f2fDispMagParaDir',f2fDispMagParaDirPos,'f2fDispMagPerpDir',f2fDispMagPerpDirPos,...
        'f2fDispSignParaDir',f2fDispSignParaDirPos,'f2fDispSignPerpDir',f2fDispSignPerpDirPos,...
        'f2fDispMagParaProt',f2fDispMagParaProtPos,'f2fDispMagPerpProt',f2fDispMagPerpProtPos,...
        'f2fDispSignParaProt',f2fDispSignParaProtPos,'f2fDispSignPerpProt',f2fDispSignPerpProtPos,...
        'ratioDispMagDir',ratioDispMagDirPos,'ratioDispSignDir',ratioDispSignDirPos,...
        'ratioDispMagProt',ratioDispMagProtPos,'ratioDispSignProt',ratioDispSignProtPos);
    
    directNeg = struct('f2fDispMag2D',f2fDispMag2DNeg,'asymParam',asymParamNeg,'angleProt',angleProtNeg,...
        'f2fDispMagParaDir',f2fDispMagParaDirNeg,'f2fDispMagPerpDir',f2fDispMagPerpDirNeg,...
        'f2fDispSignParaDir',f2fDispSignParaDirNeg,'f2fDispSignPerpDir',f2fDispSignPerpDirNeg,...
        'f2fDispMagParaProt',f2fDispMagParaProtNeg,'f2fDispMagPerpProt',f2fDispMagPerpProtNeg,...
        'f2fDispSignParaProt',f2fDispSignParaProtNeg,'f2fDispSignPerpProt',f2fDispSignPerpProtNeg,...
        'ratioDispMagDir',ratioDispMagDirNeg,'ratioDispSignDir',ratioDispSignDirNeg,...
        'ratioDispMagProt',ratioDispMagProtNeg,'ratioDispSignProt',ratioDispSignProtNeg);
    
    diffAnalysis = struct('fracUnclass',fracUnclass,'fracLin',fracLin,'fracIso',fracIso,...
        'fracIsoUnclass',fracIsoUnclass,'fracConf',fracConf,'fracBrown',fracBrown,'fracDir',fracDir,...
        'diffCoefAll',diffCoefAll,'diffCoefConf',diffCoefConf,'diffCoefBrown',diffCoefBrown,...
        'confRad',confRad);
    
    sptPropInWindow(iBandRange).directAll = directAll;
    sptPropInWindow(iBandRange).directPos = directPos;
    sptPropInWindow(iBandRange).directNeg = directNeg;
    sptPropInWindow(iBandRange).diffAnalysis = diffAnalysis;
    
end %(for iBandRange = 1 : numBandRange)

%store analysis parameters in output structure for documentation
analysisParam.bandRange = bandRange;
analysisParam.windowRange = windowRange;
analysisParam.frameRange = frameRange;
analysisParam.minTrackLen = minLength;

%% Subfunction

function [angleProt,f2fDispMag2D,f2fDispSignParaDir,f2fDispSignPerpDir,...
    f2fDispMagParaDir,f2fDispMagPerpDir,f2fDispSignParaProt,f2fDispSignPerpProt,...
    f2fDispMagParaProt,f2fDispMagPerpProt,ratioDispSignDir,ratioDispMagDir,...
    ratioDispSignProt,ratioDispMagProt,asymParam] = ...
    putTrackCharTogether(tracksCurrent,numTracksCurrent,angleWithProtTmp,...
    f2fDispTmp,paraDirDispTmp,perpDirDispTmp,paraProtDispTmp,perpProtDispTmp,...
    ratioDispDirTmp,ratioDispProtTmp,asymParamTmp,...
    angleProt,f2fDispMag2D,f2fDispSignParaDir,f2fDispSignPerpDir,...
    f2fDispMagParaDir,f2fDispMagPerpDir,f2fDispSignParaProt,f2fDispSignPerpProt,...
    f2fDispMagParaProt,f2fDispMagPerpProt,ratioDispSignDir,ratioDispMagDir,...
    ratioDispSignProt,ratioDispMagProt,asymParam,iClass)


%angle between direction of motion and protrusion vector
angleProt.mean(iClass) = nanmean(angleWithProtTmp(tracksCurrent));
angleProt.std(iClass) = nanstd(angleWithProtTmp(tracksCurrent));
angleProt.numPoints(iClass) = numTracksCurrent;

%frame-to-frame displacement
f2fDispMag2D.mean(iClass) = nanmean(f2fDispTmp(tracksCurrent));
f2fDispMag2D.std(iClass) = nanstd(f2fDispTmp(tracksCurrent));
f2fDispMag2D.numPoints(iClass) = numTracksCurrent;

%frame-to-frame displacements along and perpendicular to direction of motion
f2fDispSignParaDir.mean(iClass) = nanmean(paraDirDispTmp(tracksCurrent,1));
f2fDispSignParaDir.std(iClass) = nanstd(paraDirDispTmp(tracksCurrent,1));
f2fDispSignParaDir.numPoints(iClass) = numTracksCurrent;

f2fDispSignPerpDir.mean(iClass) = nanmean(perpDirDispTmp(tracksCurrent,1));
f2fDispSignPerpDir.std(iClass) = nanstd(perpDirDispTmp(tracksCurrent,1));
f2fDispSignPerpDir.numPoints(iClass) = numTracksCurrent;

f2fDispMagParaDir.mean(iClass) = nanmean(paraDirDispTmp(tracksCurrent,2));
f2fDispMagParaDir.std(iClass) = nanstd(paraDirDispTmp(tracksCurrent,2));
f2fDispMagParaDir.numPoints(iClass) = numTracksCurrent;

f2fDispMagPerpDir.mean(iClass) = nanmean(perpDirDispTmp(tracksCurrent,2));
f2fDispMagPerpDir.std(iClass) = nanstd(perpDirDispTmp(tracksCurrent,2));
f2fDispMagPerpDir.numPoints(iClass) = numTracksCurrent;

%frame-to-frame displacements along and perpendicular to protrusion vector
f2fDispSignParaProt.mean(iClass) = nanmean(paraProtDispTmp(tracksCurrent,1));
f2fDispSignParaProt.std(iClass) = nanstd(paraProtDispTmp(tracksCurrent,1));
f2fDispSignParaProt.numPoints(iClass) = numTracksCurrent;

f2fDispSignPerpProt.mean(iClass) = nanmean(perpProtDispTmp(tracksCurrent,1));
f2fDispSignPerpProt.std(iClass) = nanstd(perpProtDispTmp(tracksCurrent,1));
f2fDispSignPerpProt.numPoints(iClass) = numTracksCurrent;

f2fDispMagParaProt.mean(iClass) = nanmean(paraProtDispTmp(tracksCurrent,2));
f2fDispMagParaProt.std(iClass) = nanstd(paraProtDispTmp(tracksCurrent,2));
f2fDispMagParaProt.numPoints(iClass) = numTracksCurrent;

f2fDispMagPerpProt.mean(iClass) = nanmean(perpProtDispTmp(tracksCurrent,2));
f2fDispMagPerpProt.std(iClass) = nanstd(perpProtDispTmp(tracksCurrent,2));
f2fDispMagPerpProt.numPoints(iClass) = numTracksCurrent;

%ratio of perpendicular to parallel displacement components
ratioDispSignDir.mean(iClass) = nanmean(ratioDispDirTmp(tracksCurrent,1));
ratioDispSignDir.std(iClass) = nanstd(ratioDispDirTmp(tracksCurrent,1));
ratioDispSignDir.numPoints(iClass) = numTracksCurrent;

ratioDispMagDir.mean(iClass)  = nanmean(ratioDispDirTmp(tracksCurrent,2));
ratioDispMagDir.std(iClass)  = nanstd(ratioDispDirTmp(tracksCurrent,2));
ratioDispMagDir.numPoints(iClass)  = numTracksCurrent;

ratioDispSignProt.mean(iClass) = nanmean(ratioDispProtTmp(tracksCurrent,1));
ratioDispSignProt.std(iClass) = nanstd(ratioDispProtTmp(tracksCurrent,1));
ratioDispSignProt.numPoints(iClass) = numTracksCurrent;

ratioDispMagProt.mean(iClass) = nanmean(ratioDispProtTmp(tracksCurrent,2));
ratioDispMagProt.std(iClass) = nanstd(ratioDispProtTmp(tracksCurrent,2));
ratioDispMagProt.numPoints(iClass)  = numTracksCurrent;

%asymmetry parameter (ratio of maximum to minimum eigenvalues of position
%covariance matrix)
asymParam.mean(iClass) = nanmean(asymParamTmp(tracksCurrent));
asymParam.std(iClass) = nanstd(asymParamTmp(tracksCurrent));
asymParam.numPoints(iClass) = numTracksCurrent;


%% ~~~ the end ~~~
