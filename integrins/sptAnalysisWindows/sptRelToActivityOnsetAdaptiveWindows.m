function [sptPropInWindow,windowTrackAssignExt,windowSize,analysisParam] = ...
    sptRelToActivityOnsetAdaptiveWindows(tracksFinal,diffAnalysisRes,...
    diffModeAnRes,directTrackChar,winPositions,winFrames,protSamples,...
    minLength,indxSlices,frameRange,windowTrackAssignExt,firstMaskFile)
%sptRelToActivityOnsetAdaptiveWindows calculates single particle behavior in adaptive windows grouped based on edge activity
%
%SYNOPSIS [sptPropInWindow,windowTrackAssignExt,windowSize,analysisParam] = ...
%    sptRelToActivityOnsetAdaptiveWindows(tracksFinal,diffAnalysisRes,...
%    diffModeAnRes,directTrackChar,winPositions,winFrames,protSamples,...
%    minLength,indxSlices,frameRange,windowTrackAssignExt,firstMaskFile)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%       diffModeAnRes  : Output of trackDiffModeAnalysis.
%       directTrackChar: Output of trackMotionCharProtrusion.
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
%       minLength      : Minimum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: 5.
%       indxSlices     : Vector with indices of slices (i.e. window number 
%                        along cell edge) to include in analysis.
%                        Optional. Default: all windows.
%       frameRange     : Vector with 2 entries indicating range of window
%                        frames to include in analysis.
%                        Optional. Default: all frames.
%       windowTrackAssignExt: Output of assignTracks2Windows.
%                        Optional. If not input, will be calculated.
%       firstMaskFile  : Name (including full path) of first mask file.
%                        Optional. If not input, will be asked for.
%
%OUTPUT sptPropInWindow: Structure array with N entries (N from bandRange)
%                        storing various single particle properties, with
%                        four fields:
%           .directAll      : Properties calculated directly from tracks,
%                             for all tracks.
%           .directPos      : Properties calculated directly from tracks,
%                             for tracks with net displacement parallel to
%                             protrusion vector.
%           .directNeg      : Properties calculated directly from tracks,
%                             for tracks with net displacement
%                             anti-parallel to protrusion vector.
%           .diffAnalysis   : Properties extracted by diffusion analysis.
%                     ***directAll, directPos and directNeg contain the
%                        following sub-fields:
%               .spDensity     : Single particle density.
%               .fracNetDispNeg: Fraction of tracks with net displacement
%                                anti-parallel to protrusion vector.
%               .f2fDispMag2D  : Frame-to-frame displacement.
%               .asymParam     : Asymmtry measure.
%               .angleProt     : Angle between direction of motion and
%                                protrusion vector.
%               .f2fDispMagParaDir: Frame-to-frame displacement component
%                                parallel to direction of motion (no sign
%                                information per track).
%               .f2fDispMagPerpDir: Frame-to-frame displacement component
%                                perpendicular to direction of motion (no
%                                sign information per track).
%               .f2fDispSignParaDir: Frame-to-frame displacement component
%                                parallel to direction of motion (with sign
%                                information per track).
%               .f2fDispSignPerpDir: Frame-to-frame displacement component
%                                perpendicular to direction of motion (with
%                                sign information per track).
%               .f2fDispMagParaProt: Frame-to-frame displacement component
%                                parallel to protrusion vector (no sign
%                                information per track).
%               .f2fDispMagPerpProt: Frame-to-frame displacement component
%                                perpendicular to protrusion vector (no
%                                sign information per track).
%               .f2fDispSignParaProt: Frame-to-frame displacement component
%                                parallel to protrusion vector (with sign
%                                information per track).
%               .f2fDispSignPerpProt: Frame-to-frame displacement component
%                                perpendicular to protrusion vector (with
%                                sign information per track).
%               .ratioDispMagDir: Ratio of parallel to perpendicular
%                                component.
%               .ratioDispSignDir: Ratio of parallel to perpendicular
%                                component.
%               .ratioDispMagProt: Ratio of parallel to perpendicular
%                                component.
%               .ratioDispSignProt: Ratio of parallel to perpendicular
%                                component.
%                     ***diffAnalysis contains the following sub-fields:
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
%           .diffCoefAll    : Diffusion coefficient of both confined and
%                             Brownian tracks.
%           .diffCoefConf   : Diffusion coefficient of confined tracks.
%           .diffCoefBrown  : Diffusion coefficient of Brownian tracks.
%           .confRad        : Confinement radius of confined tracks.
%                        Each field contains 3 sub-fields, "mean", "std" and
%                        "numPoints" storing the mean and std of each property
%                        and the number of points contributing to each
%                        value.
%                        Each of these fields is a 2D array.
%                        Each row refers to an activity type, as listed in
%                        groupWindowsActivity.
%                        The columns refer to increment from activity
%                        onset, starting with -3.
%       windowTrackAssignExt: Output of assignTracks2Windows.
%                        If input, same as input.
%       windowSize     : 3-D matrix of dimensions (number of bands) x (number
%                        of slices) x (number of window frames - 1)
%                        storing the size of each window. NaN indicates a
%                        window of zero size.
%       analysisParam  : Parameters used in analysis.
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

if nargin < 7
    disp('--sptRelToActivityOnsetAdaptiveWindows: Missing input arguments!');
    return
end

if nargin < 8 || isempty(minLength)
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

if nargin < 9 || isempty(indxSlices)
    indxSlices = (1:numWinPara)';
end

if nargin < 10 || isempty(frameRange)
    frameRange = [1 numWinFrames];
end

if nargin < 11 || isempty(windowTrackAssignExt)
    windowTrackAssignExt = [];
end

if nargin < 12 || iesmpty(firstMaskFile)
    [fName,dirName] = uigetfile('*.tif','specify first cell mask file to use in analysis.');
    firstMaskFile = fullfile(dirName,fName);
end

%% Trajectory pre-processing

%keep only trajectories longer than minLength
criteria.lifeTime.min = minLength;
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx,:);
diffAnalysisRes = diffAnalysisRes(indx);
diffModeAnRes = diffModeAnRes(indx);
directTrackChar = directTrackChar(indx);

%divide trajectories among windows
if isempty(windowTrackAssignExt)
    [~,~,~,windowTrackAssignExt] = assignTracks2Windows(...
        tracksFinal,winPositions,winFrames,1);
end

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

%group window slices based on activity type
sliceActivityGroup = groupWindowsActivity(protSamples,0,indxSlices,frameRange);

%generate adaptive window combinations
protrusionCombinedWindows = combineWindowsProtrusion(sliceActivityGroup,...
    winPositions,firstMaskFile);

%put together the tracks that belong to each window combination
[eventCombinedTracks,eventWindowsList] = combineTracksManyWindows(protrusionCombinedWindows,...
    windowTrackAssignExt);

%% Particle behavior pre-processing

%get the lifetime of each track segment
trackLft = getTrackSEL(tracksFinal,1);
trackLft = trackLft(:,3);

%get the number of track segments
numSegments = length(trackLft);

% FROM ASYMMETRY AND DIFFUSION ANALYSIS ...

%get trajectory classifications
trajClassDiff = vertcat(diffAnalysisRes.classification);
trajClassAsym = trajClassDiff(:,1); %asymmetry classification
trajClassDiff = trajClassDiff(:,2); %diffusion classification

%WARNING
%MAKE trajClassAsym ALL NaN, SO THAT THERE IS NO ASYMMETRY CLASSIFICATION
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

% FROM DIFFUSION MODE ANALYSIS ...

%get trajectory diffusion modes
trajDiffMode = vertcat(diffModeAnRes.diffMode);

%get trajectory diffusion coefficient
trajDiffCoef2 = vertcat(diffModeAnRes.diffCoef);

%get number of diffusion modes
numDiffMode = max(trajDiffMode);

% FROM TRACKS DIRECTLY ...

%get direction of motion, angle with protrusion vector and various
%frame-to-frame displacement measures
angleWithProtTmp = vertcat(directTrackChar.angleWithProt);
f2fDispTmp       = vertcat(directTrackChar.f2fDisp);
paraDirDispTmp   = vertcat(directTrackChar.paraDirDisp);
perpDirDispTmp   = vertcat(directTrackChar.perpDirDisp);
paraProtDispTmp  = vertcat(directTrackChar.paraProtDisp);
perpProtDispTmp  = vertcat(directTrackChar.perpProtDisp);
asymParamTmp     = vertcat(directTrackChar.asymParam);

%calculate ratio of parallel to perpendicular displacements
ratioDispDirTmp = abs( paraDirDispTmp ./ perpDirDispTmp );
ratioDispProtTmp = abs( paraProtDispTmp ./ perpProtDispTmp );

%% Calculate property values per window group

%number of activity types
numTypes = length(eventCombinedTracks);

%get fields in the structure eventCombinedTracks
eventField = fieldnames(eventCombinedTracks);
numFields = length(eventField);

%initialize output variable
sptPropInWindow = repmat(struct('directAll',[],'directPos',[],...
    'directNeg',[],'diffMSSAnalysis',[],'diffModeAnalysis',[]),numTypes,1);

%go over activity types
for iType = 1 : numTypes
    
    %get current activity
    activityCurrent = eventCombinedTracks(iType);
    
    %check whether there are any events in this activity
    tmp = activityCurrent.(eventField{1});
    eventsExist = ~isempty(tmp{:});
    
    %if there are events
    if eventsExist
        
        %initialize variables storing results
        for iField = 1 : numFields
            nanTmp = NaN(size(activityCurrent.(eventField{iField})));
            spDensity.(eventField{iField}).mean = nanTmp;
            spDensity.(eventField{iField}).std = nanTmp;
            spDensity.(eventField{iField}).numPoints = nanTmp;
        end        
        
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
            fracNetDispNeg,...
            fracUnclass,fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
            diffCoefAll,diffCoefConf,diffCoefBrown,confRad,diffCoefModeAll] = ...
            ...
            deal(spDensity);
        
        for iField = 1 : numFields
            nanTmp = NaN([size(activityCurrent.(eventField{iField})) numDiffMode+1]);
            fracMode.(eventField{iField}).mean = nanTmp;
            fracMode.(eventField{iField}).std = nanTmp;
            fracMode.(eventField{iField}).numPoints = nanTmp;
        end
        modeDensity = fracMode;
        
        for iField = 1 : numFields
            nanTmp = NaN([size(activityCurrent.(eventField{iField})) numDiffMode]);
            diffCoefModeInd.(eventField{iField}).mean = nanTmp;
            diffCoefModeInd.(eventField{iField}).std = nanTmp;
            diffCoefModeInd.(eventField{iField}).numPoints = nanTmp;
        end        
        
        %go over all fields in activityCurrent (indicating various
        %series relative to activity onset)
        for iField = 1 : numFields
            
            %get the tracks corresponding to this activity type and field
            tracksCurrentAll = activityCurrent.(eventField{iField});
            
            %get the windows they belong to
            windowListCurrentAll = eventWindowsList(iType).(eventField{iField});
            
            %extract maximum increment (whether positive or negative)
            [maxInc,numCol] = size(tracksCurrentAll);
            
            %go over the increments and calculate particle properties in
            %windows at those increments
            for iCol = 1 : numCol
                for iInc = 1 : maxInc
                    
                    %get tracks
                    tracksCurrent = tracksCurrentAll{iInc,iCol};
                    numTracksCurrent = length(tracksCurrent);
                    
                    %get windows
                    windowListCurrent = windowListCurrentAll{iInc,iCol};
                    
                    %if there are tracks
                    if numTracksCurrent > 0
                        
                        %get particle behavior
                        [spDensity,fracNetDispNeg,...
                            ...
                            angleProtAll,asymParamAll,f2fDispMag2DAll,...
                            f2fDispSignParaDirAll,f2fDispSignPerpDirAll,f2fDispMagParaDirAll,f2fDispMagPerpDirAll,...
                            f2fDispSignParaProtAll,f2fDispSignPerpProtAll,f2fDispMagParaProtAll,f2fDispMagPerpProtAll,...
                            ratioDispSignDirAll,ratioDispMagDirAll,ratioDispSignProtAll,ratioDispMagProtAll,...
                            ...
                            angleProtPos,asymParamPos,f2fDispMag2DPos,...
                            f2fDispSignParaDirPos,f2fDispSignPerpDirPos,f2fDispMagParaDirPos,f2fDispMagPerpDirPos,...
                            f2fDispSignParaProtPos,f2fDispSignPerpProtPos,f2fDispMagParaProtPos,f2fDispMagPerpProtPos,...
                            ratioDispSignDirPos,ratioDispMagDirPos,ratioDispSignProtPos,ratioDispMagProtPos,...
                            ...
                            angleProtNeg,asymParamNeg,f2fDispMag2DNeg,...
                            f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,...
                            f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,...
                            ratioDispSignDirNeg,ratioDispMagDirNeg,ratioDispSignProtNeg,ratioDispMagProtNeg,...
                            ...
                            fracUnclass,fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
                            diffCoefAll,diffCoefConf,diffCoefBrown,confRad,...
                            modeDensity,fracMode,diffCoefModeInd,diffCoefModeAll] = ...
                            ...
                            getSliceActivityGroupChar...
                            ...
                            (eventField{iField},iInc,iCol,numSPTFrames,...
                            winSize,windowListCurrent,tracksCurrent,trackLft,...
                            trajClass,diffCoefGen,confRadAll,angleWithProtTmp,asymParamTmp,f2fDispTmp,...
                            paraDirDispTmp,perpDirDispTmp,paraProtDispTmp,perpProtDispTmp,...
                            ratioDispDirTmp,ratioDispProtTmp,trajDiffMode,trajDiffCoef2,...
                            ...
                            spDensity,fracNetDispNeg,...
                            ...
                            angleProtAll,asymParamAll,f2fDispMag2DAll,...
                            f2fDispSignParaDirAll,f2fDispSignPerpDirAll,f2fDispMagParaDirAll,f2fDispMagPerpDirAll,...
                            f2fDispSignParaProtAll,f2fDispSignPerpProtAll,f2fDispMagParaProtAll,f2fDispMagPerpProtAll,...
                            ratioDispSignDirAll,ratioDispMagDirAll,ratioDispSignProtAll,ratioDispMagProtAll,...
                            ...
                            angleProtPos,asymParamPos,f2fDispMag2DPos,...
                            f2fDispSignParaDirPos,f2fDispSignPerpDirPos,f2fDispMagParaDirPos,f2fDispMagPerpDirPos,...
                            f2fDispSignParaProtPos,f2fDispSignPerpProtPos,f2fDispMagParaProtPos,f2fDispMagPerpProtPos,...
                            ratioDispSignDirPos,ratioDispMagDirPos,ratioDispSignProtPos,ratioDispMagProtPos,...
                            ...
                            angleProtNeg,asymParamNeg,f2fDispMag2DNeg,...
                            f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,...
                            f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,...
                            ratioDispSignDirNeg,ratioDispMagDirNeg,ratioDispSignProtNeg,ratioDispMagProtNeg,...
                            ...
                            fracUnclass,fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
                            diffCoefAll,diffCoefConf,diffCoefBrown,confRad,...
                            modeDensity,fracMode,diffCoefModeInd,diffCoefModeAll);
                        
                    end %(if numTracksCurrent > 0)
                    
                end %(for iInc = 1 : maxInc)
            end %(for iCol = 1 : numCol)
            
        end %(for iField = 1 : numFields)
        
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
        
        modeAnalysis = struct('modeDensity',modeDensity,'fracMode',fracMode,...
            'diffCoefModeInd',diffCoefModeInd,'diffCoefModeAll',diffCoefModeAll);
        
        sptPropInWindow(iType).directAll = directAll;
        sptPropInWindow(iType).directPos = directPos;
        sptPropInWindow(iType).directNeg = directNeg;
        sptPropInWindow(iType).diffMSSAnalysis = diffAnalysis;
        sptPropInWindow(iType).diffModeAnalysis = modeAnalysis;
        
    end %(if eventsExist)
    
end %(for iType = 1 : numTypes)

%store analysis parameters in output structure for documentation
analysisParam.indxSlices = indxSlices;
analysisParam.frameRange = frameRange;
analysisParam.minTrackLen = minLength;


%% Subfunction 1 "getsliceActivityGroupChar"

function [spDensity,fracNetDispNeg,...
    ...
    angleProtAll,asymParamAll,f2fDispMag2DAll,...
    f2fDispSignParaDirAll,f2fDispSignPerpDirAll,f2fDispMagParaDirAll,f2fDispMagPerpDirAll,...
    f2fDispSignParaProtAll,f2fDispSignPerpProtAll,f2fDispMagParaProtAll,f2fDispMagPerpProtAll,...
    ratioDispSignDirAll,ratioDispMagDirAll,ratioDispSignProtAll,ratioDispMagProtAll,...
    ...
    angleProtPos,asymParamPos,f2fDispMag2DPos,...
    f2fDispSignParaDirPos,f2fDispSignPerpDirPos,f2fDispMagParaDirPos,f2fDispMagPerpDirPos,...
    f2fDispSignParaProtPos,f2fDispSignPerpProtPos,f2fDispMagParaProtPos,f2fDispMagPerpProtPos,...
    ratioDispSignDirPos,ratioDispMagDirPos,ratioDispSignProtPos,ratioDispMagProtPos,...
    ...
    angleProtNeg,asymParamNeg,f2fDispMag2DNeg,...
    f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,...
    f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,...
    ratioDispSignDirNeg,ratioDispMagDirNeg,ratioDispSignProtNeg,ratioDispMagProtNeg,...
    ...
    fracUnclass,fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
    diffCoefAll,diffCoefConf,diffCoefBrown,confRad,...
    modeDensity,fracMode,diffCoefModeInd,diffCoefModeAll] = ...
    ...
    getSliceActivityGroupChar...
    ...
    (eventField,iInc,iCol,numSPTFrames,...
    winSize,windowListCurrent,tracksCurrent,trackLft,...
    trajClass,diffCoefGen,confRadAll,angleWithProtTmp,asymParamTmp,f2fDispTmp,...
    paraDirDispTmp,perpDirDispTmp,paraProtDispTmp,perpProtDispTmp,...
    ratioDispDirTmp,ratioDispProtTmp,trajDiffMode,trajDiffCoef2,...
    ...
    spDensity,fracNetDispNeg,...
    ...
    angleProtAll,asymParamAll,f2fDispMag2DAll,...
    f2fDispSignParaDirAll,f2fDispSignPerpDirAll,f2fDispMagParaDirAll,f2fDispMagPerpDirAll,...
    f2fDispSignParaProtAll,f2fDispSignPerpProtAll,f2fDispMagParaProtAll,f2fDispMagPerpProtAll,...
    ratioDispSignDirAll,ratioDispMagDirAll,ratioDispSignProtAll,ratioDispMagProtAll,...
    ...
    angleProtPos,asymParamPos,f2fDispMag2DPos,...
    f2fDispSignParaDirPos,f2fDispSignPerpDirPos,f2fDispMagParaDirPos,f2fDispMagPerpDirPos,...
    f2fDispSignParaProtPos,f2fDispSignPerpProtPos,f2fDispMagParaProtPos,f2fDispMagPerpProtPos,...
    ratioDispSignDirPos,ratioDispMagDirPos,ratioDispSignProtPos,ratioDispMagProtPos,...
    ...
    angleProtNeg,asymParamNeg,f2fDispMag2DNeg,...
    f2fDispSignParaDirNeg,f2fDispSignPerpDirNeg,f2fDispMagParaDirNeg,f2fDispMagPerpDirNeg,...
    f2fDispSignParaProtNeg,f2fDispSignPerpProtNeg,f2fDispMagParaProtNeg,f2fDispMagPerpProtNeg,...
    ratioDispSignDirNeg,ratioDispMagDirNeg,ratioDispSignProtNeg,ratioDispMagProtNeg,...
    ...
    fracUnclass,fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
    diffCoefAll,diffCoefConf,diffCoefBrown,confRad,...
    modeDensity,fracMode,diffCoefModeInd,diffCoefModeAll)

%get number of tracks
numTracksCurrent = length(tracksCurrent);

%calculate the sum of all window sizes in this window group
[numWinPerp,numWinPara,numWinFramesM1] = size(winSize);
linearInd = sub2ind([numWinPerp numWinPara numWinFramesM1],...
    windowListCurrent(:,1),windowListCurrent(:,2),windowListCurrent(:,3));
windowGroupSize = nansum(winSize(linearInd));

if numTracksCurrent == 0 %if there are no tracks in this window ...
    
    %calculate particle density
    %note that here particle density is not automatically 0
    %because it can be that the windows are collapsed (stored as
    %NaN), in which case the density will be NaN
    spDensity.(eventField).mean(iInc,iCol) = 0 / windowGroupSize;
    spDensity.(eventField).numPoints(iInc,iCol) = numTracksCurrent;
    
else %if there are tracks in this window ...
    
    % FROM THE TRACKS DIRECTLY ...
    
    %calculate overall particle density
    %(the density standard deviation will be calculated later in combination
    %with the mode density standard deviation)
    spDensity.(eventField).mean(iInc,iCol) = sum(trackLft(tracksCurrent)) / ...
        windowGroupSize / numSPTFrames;
    spDensity.(eventField).numPoints(iInc,iCol) = numTracksCurrent;
    
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
        ratioDispSignProtAll,ratioDispMagProtAll,asymParamAll,...
        iInc,iCol,eventField);
    
    %divide current tracks into those with net displacement
    %parallel or opposite to protrusion vector
    paraProtDispCurrent = paraProtDispTmp(tracksCurrent,1);
    tracksPos = tracksCurrent(paraProtDispCurrent >= 0);
    tracksNeg = tracksCurrent(paraProtDispCurrent < 0);
    numTracksPos = length(tracksPos);
    numTracksNeg = length(tracksNeg);
    
    %store fraction of tracks with net displacement opposite to
    %protrusion vector (i.e. negative net displacement)
    fracNetDispNeg.(eventField).mean(iInc,iCol) = numTracksNeg / numTracksCurrent;
    fracNetDispNeg.(eventField).numPoints(iInc,iCol) = numTracksCurrent;
    
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
        ratioDispSignProtPos,ratioDispMagProtPos,asymParamPos,...
        iInc,iCol,eventField);
    
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
        ratioDispSignProtNeg,ratioDispMagProtNeg,asymParamNeg,...
        iInc,iCol,eventField);
    
    % FROM THE ASYMMETRY AND DIFFUSION ANALYSIS ...
    
    %calculate the fraction of tracks in each motion category
    
    %first completely unclassified tracks
    trajClassCurrent = trajClass(tracksCurrent);
    fracUnclass.(eventField).mean(iInc,iCol) = ...
        length(find(isnan(trajClassCurrent)))/numTracksCurrent;
    fracUnclass.(eventField).numPoints(iInc,iCol) = numTracksCurrent;
    
    %then tracks with asymmetry analysis
    trajClassCurrent = trajClassCurrent(~isnan(trajClassCurrent));
    numTracksAsym = length(trajClassCurrent);
    if numTracksAsym > 0
        
        %tracks classified as linear or isotropic
        fracLin.(eventField).mean(iInc,iCol) = ...
            length(find(trajClassCurrent==5))/numTracksAsym;
        fracLin.(eventField).numPoints(iInc,iCol) = numTracksAsym;
        fracIso.(eventField).mean(iInc,iCol) = ...
            length(find(trajClassCurrent~=5))/numTracksAsym;
        fracIso.(eventField).numPoints(iInc,iCol) = numTracksAsym;
        
        %within the tracks classified as isotropic
        trajClassCurrent = trajClassCurrent(trajClassCurrent~=5);
        numTracksIso = length(trajClassCurrent);
        if numTracksIso > 0
            
            %tracks without diffusion analysis
            fracIsoUnclass.(eventField).mean(iInc,iCol) = ...
                length(find(trajClassCurrent==4))/numTracksIso;
            fracIsoUnclass.(eventField).numPoints(iInc,iCol) = numTracksIso;
            
            %finally tracks with diffusion analysis
            trajClassCurrent = trajClassCurrent(trajClassCurrent~=4);
            numTracksDiff = length(trajClassCurrent);
            if numTracksDiff > 0
                
                fracClass = hist(trajClassCurrent,(1:3))/numTracksDiff;
                fracConf.(eventField).mean(iInc,iCol) = fracClass(1);
                fracConf.(eventField).numPoints(iInc,iCol) = numTracksDiff;
                fracBrown.(eventField).mean(iInc,iCol) = fracClass(2);
                fracBrown.(eventField).numPoints(iInc,iCol) = numTracksDiff;
                fracDir.(eventField).mean(iInc,iCol) = fracClass(3);
                fracDir.(eventField).numPoints(iInc,iCol) = numTracksDiff;
                
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
        diffCoefAll.(eventField).mean(iInc,iCol) = mean(tmpVecAll);
        diffCoefAll.(eventField).std(iInc,iCol) = std(tmpVecAll);
        diffCoefAll.(eventField).numPoints(iInc,iCol) = length(tmpVecAll);
    end
    if ~isempty(tmpVecConf)
        diffCoefConf.(eventField).mean(iInc,iCol) = mean(tmpVecConf);
        diffCoefConf.(eventField).std(iInc,iCol) = std(tmpVecConf);
        diffCoefConf.(eventField).numPoints(iInc,iCol) = length(tmpVecConf);
    end
    if ~isempty(tmpVecBrown)
        diffCoefBrown.(eventField).mean(iInc,iCol) = mean(tmpVecBrown);
        diffCoefBrown.(eventField).std(iInc,iCol) = std(tmpVecBrown);
        diffCoefBrown.(eventField).numPoints(iInc,iCol) = length(tmpVecBrown);
    end
    
    %calculate the average confinement radius for confined
    %particles
    tmpVec = confRadAll(tracksCurrent);
    tmpVec = tmpVec(trajClass(tracksCurrent) == 1);
    if ~isempty(tmpVec)
        confRad.(eventField).mean(iInc,iCol) = mean(tmpVec);
        confRad.(eventField).std(iInc,iCol) = std(tmpVec);
        confRad.(eventField).numPoints(iInc,iCol) = length(tmpVec);
    end
    
    % FROM THE DIFFUSION MODE ANALYSIS ...
    
    %get number of modes
    numDiffMode = size(fracMode.(eventField).mean,3) - 1;
    
    %get diffusion coefficient and diffusion mode of all tracks
    trajModeCurrent = trajDiffMode(tracksCurrent);
    trajModeCurrent(isnan(trajModeCurrent)) = numDiffMode + 1;
    trajD2Current = trajDiffCoef2(tracksCurrent);
    
    %calculate the average diffusion coefficient for all tracks
    %regardless of their mode
    tmpVec = trajD2Current(~isnan(trajD2Current));
    if ~isempty(tmpVec)
        diffCoefModeAll.(eventField).mean(iInc,iCol) = mean(tmpVec);
        diffCoefModeAll.(eventField).std(iInc,iCol) = std(tmpVec);
        diffCoefModeAll.(eventField).numPoints(iInc,iCol) = length(tmpVec);
    end
    
    %calculate the average diffusion coefficient per mode
    for iMode = 1 : numDiffMode
        tmpVec = trajD2Current(trajModeCurrent==iMode);
        if ~isempty(tmpVec)
            diffCoefModeInd.(eventField).mean(iInc,iCol,iMode) = mean(tmpVec);
            diffCoefModeInd.(eventField).std(iInc,iCol,iMode) = std(tmpVec);
            diffCoefModeInd.(eventField).numPoints(iInc,iCol,iMode) = length(tmpVec);
        end
    end
    
    %calculate the density of particles in each diffusion mode
    %also calculate the fraction of particles in each mode
    for iMode = 1 : numDiffMode+1
        tracksMode = tracksCurrent(trajModeCurrent==iMode);
        if isempty(tracksMode)
            modeDensity.(eventField).mean(iInc,iCol,iMode) = 0 / windowGroupSize;
            modeDensity.(eventField).numPoints(iInc,iCol,iMode) = numTracksCurrent;
        else
            modeDensity.(eventField).mean(iInc,iCol,iMode) = ...
                sum(trackLft(tracksMode)) / windowGroupSize / numSPTFrames;
            modeDensity.(eventField).numPoints(iInc,iCol,iMode) = numTracksCurrent;
        end
        fracMode.(eventField).mean(iInc,iCol,iMode) = ...
            modeDensity.(eventField).mean(iInc,iCol,iMode) / spDensity.(eventField).mean(iInc,iCol);
        fracMode.(eventField).numPoints(iInc,iCol,iMode) = numTracksCurrent;
    end
    
    %estimate the standard deviations of overall density, mode density and
    %mode fraction from a bootstrap sample
    densityBoot = NaN(200,numDiffMode+1);
    for iBoot = 1 : 200
        indxBoot = randsample(numTracksCurrent,numTracksCurrent,true);
        tracksBoot = tracksCurrent(indxBoot);
        modeBoot = trajModeCurrent(indxBoot);
        for iMode = 1 : numDiffMode+1
            tracksMode = tracksBoot(modeBoot==iMode);
            if isempty(tracksMode)
                densityBoot(iBoot,iMode) = 0 / windowGroupSize;
            else
                densityBoot(iBoot,iMode) = ...
                    sum(trackLft(tracksMode)) / windowGroupSize / numSPTFrames;
            end
        end
    end
    spDensityBoot = sum(densityBoot,2);
    fracModeBoot = densityBoot ./ repmat(spDensityBoot,1,numDiffMode+1);
    modeDensity.(eventField).std(iInc,iCol,:) = std(densityBoot) * sqrt(numTracksCurrent);
    fracMode.(eventField).std(iInc,iCol,:) = std(fracModeBoot) * sqrt(numTracksCurrent);    
    spDensity.(eventField).std(iInc,iCol) = std(spDensityBoot) * sqrt(numTracksCurrent);
    
        
    %     %calculate the fraction of tracks in each diffusion mode
    %     %also calculate the average diffusion coefficient per mode
    %
    %     %first completely unclassified tracks
    %     trajModeCurrent = trajDiffMode(tracksCurrent);
    %     fracMode.(eventField).mean(iInc,iCol,end) = ...
    %         length(find(isnan(trajModeCurrent)))/numTracksCurrent;
    %     fracMode.(eventField).numPoints(iInc,iCol,end) = numTracksCurrent;
    %
    %     %then classified tracks
    %     tracksCurrentModeClass = tracksCurrent(~isnan(trajModeCurrent));
    %     numTracksModeClass = length(tracksCurrentModeClass);
    %
    %     if numTracksModeClass > 0
    %
    %         %extract diffusion mode and diffusion coefficient for relevant
    %         %tracks
    %         tmpVecMode = trajDiffMode(tracksCurrentModeClass);
    %         tmpVecCoef = trajDiffCoef2(tracksCurrentModeClass);
    %
    %         %generate a bootstrap sample to calculate standard deviation of
    %         %fractions
    %         modeVecBoot = NaN(numTracksModeClass,200);
    %         for iBoot = 1 : 200
    %             indxBoot = randsample(numTracksModeClass,numTracksModeClass,true);
    %             modeVecBoot(:,iBoot) = tmpVecMode(indxBoot);
    %         end
    %
    %         %calculate fraction of tracks in each mode
    %         n = hist(tmpVecMode,1:numDiffMode);
    %         fracMode.(eventField).mean(iInc,iCol,1:end-1) = n/numTracksModeClass;
    %         fracMode.(eventField).numPoints(iInc,iCol,1:end-1) = numTracksModeClass;
    %
    %         %get standard deviations from bootstrap sample
    %         n = hist(modeVecBoot,1:numDiffMode)';
    %         fracBoot = n/numTracksModeClass;
    %         fracMode.(eventField).std(iInc,iCol,1:end-1) = ...
    %             std(fracBoot) * sqrt(numTracksModeClass);
    %
    %         %calculate the diffusion coefficient of each mode
    %         for iMode = 1 : numDiffMode
    %             tmpTmpVec = tmpVecCoef(tmpVecMode==iMode);
    %             diffCoefModeInd.(eventField).mean(iInc,iCol,iMode) = mean(tmpTmpVec);
    %             diffCoefModeInd.(eventField).std(iInc,iCol,iMode) = std(tmpTmpVec);
    %             diffCoefModeInd.(eventField).numPoints(iInc,iCol,iMode) = length(tmpTmpVec);
    %         end
    %
    %     end %(if numTracksModeClass > 0)
    
end %(if numTracksCurrent == 0)


%% Subfunction 2 "putTrackCharTogether"

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
    ratioDispSignProt,ratioDispMagProt,asymParam,iInc,iCol,eventField)


%angle between direction of motion and protrusion vector
angleProt.(eventField).mean(iInc,iCol) = nanmean(angleWithProtTmp(tracksCurrent));
angleProt.(eventField).std(iInc,iCol) = nanstd(angleWithProtTmp(tracksCurrent));
angleProt.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

%frame-to-frame displacement
f2fDispMag2D.(eventField).mean(iInc,iCol) = nanmean(f2fDispTmp(tracksCurrent));
f2fDispMag2D.(eventField).std(iInc,iCol) = nanstd(f2fDispTmp(tracksCurrent));
f2fDispMag2D.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

%frame-to-frame displacements along and perpendicular to direction of motion
f2fDispSignParaDir.(eventField).mean(iInc,iCol) = nanmean(paraDirDispTmp(tracksCurrent,1));
f2fDispSignParaDir.(eventField).std(iInc,iCol) = nanstd(paraDirDispTmp(tracksCurrent,1));
f2fDispSignParaDir.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

f2fDispSignPerpDir.(eventField).mean(iInc,iCol) = nanmean(perpDirDispTmp(tracksCurrent,1));
f2fDispSignPerpDir.(eventField).std(iInc,iCol) = nanstd(perpDirDispTmp(tracksCurrent,1));
f2fDispSignPerpDir.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

f2fDispMagParaDir.(eventField).mean(iInc,iCol) = nanmean(paraDirDispTmp(tracksCurrent,2));
f2fDispMagParaDir.(eventField).std(iInc,iCol) = nanstd(paraDirDispTmp(tracksCurrent,2));
f2fDispMagParaDir.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

f2fDispMagPerpDir.(eventField).mean(iInc,iCol) = nanmean(perpDirDispTmp(tracksCurrent,2));
f2fDispMagPerpDir.(eventField).std(iInc,iCol) = nanstd(perpDirDispTmp(tracksCurrent,2));
f2fDispMagPerpDir.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

%frame-to-frame displacements along and perpendicular to protrusion vector
f2fDispSignParaProt.(eventField).mean(iInc,iCol) = nanmean(paraProtDispTmp(tracksCurrent,1));
f2fDispSignParaProt.(eventField).std(iInc,iCol) = nanstd(paraProtDispTmp(tracksCurrent,1));
f2fDispSignParaProt.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

f2fDispSignPerpProt.(eventField).mean(iInc,iCol) = nanmean(perpProtDispTmp(tracksCurrent,1));
f2fDispSignPerpProt.(eventField).std(iInc,iCol) = nanstd(perpProtDispTmp(tracksCurrent,1));
f2fDispSignPerpProt.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

f2fDispMagParaProt.(eventField).mean(iInc,iCol) = nanmean(paraProtDispTmp(tracksCurrent,2));
f2fDispMagParaProt.(eventField).std(iInc,iCol) = nanstd(paraProtDispTmp(tracksCurrent,2));
f2fDispMagParaProt.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

f2fDispMagPerpProt.(eventField).mean(iInc,iCol) = nanmean(perpProtDispTmp(tracksCurrent,2));
f2fDispMagPerpProt.(eventField).std(iInc,iCol) = nanstd(perpProtDispTmp(tracksCurrent,2));
f2fDispMagPerpProt.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

%ratio of perpendicular to parallel displacement components
ratioDispSignDir.(eventField).mean(iInc,iCol) = nanmean(ratioDispDirTmp(tracksCurrent,1));
ratioDispSignDir.(eventField).std(iInc,iCol) = nanstd(ratioDispDirTmp(tracksCurrent,1));
ratioDispSignDir.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

ratioDispMagDir.(eventField).mean(iInc,iCol)  = nanmean(ratioDispDirTmp(tracksCurrent,2));
ratioDispMagDir.(eventField).std(iInc,iCol)  = nanstd(ratioDispDirTmp(tracksCurrent,2));
ratioDispMagDir.(eventField).numPoints(iInc,iCol)  = numTracksCurrent;

ratioDispSignProt.(eventField).mean(iInc,iCol) = nanmean(ratioDispProtTmp(tracksCurrent,1));
ratioDispSignProt.(eventField).std(iInc,iCol) = nanstd(ratioDispProtTmp(tracksCurrent,1));
ratioDispSignProt.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

ratioDispMagProt.(eventField).mean(iInc,iCol) = nanmean(ratioDispProtTmp(tracksCurrent,2));
ratioDispMagProt.(eventField).std(iInc,iCol) = nanstd(ratioDispProtTmp(tracksCurrent,2));
ratioDispMagProt.(eventField).numPoints(iInc,iCol)  = numTracksCurrent;

%asymmetry parameter (ratio of maximum to minimum eigenvalues of position
%covariance matrix)
asymParam.(eventField).mean(iInc,iCol) = nanmean(asymParamTmp(tracksCurrent));
asymParam.(eventField).std(iInc,iCol) = nanstd(asymParamTmp(tracksCurrent));
asymParam.(eventField).numPoints(iInc,iCol) = numTracksCurrent;

%% ~~~ the end ~~~
