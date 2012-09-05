function windowNumbersAssignExt = assignNumbers2Windows(tracksFinal,...
    diffAnalysisRes,diffModeAnRes,directTrackChar,winPositions,winFrames,...
    windowTrackAssignExt,lengthMinMax,prop2calc,windowNumbersAssignExtIn)
%assignNumbers2Windows calculates various particle properties within spatial and temporal windows derived from the cell edge
%
%SYNOPSIS [windowNumbersAssignExt] = assignNumbers2Windows(tracksFinal,...
%    diffAnalysisRes,diffModeAnRes,directTrackChar,winPositions,winFrames,...
%    windowTrackAssignExt,lengthMinMax)
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
%       winFrames      : The "SPT frames" at which there are windows.
%       windowTrackAssignExt: Output of assignTracks2Windows.
%       lengthMinMax   : Row vector with 2 entries indicating minimum and
%                        maximum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: [5 99].
%       prop2calc      : Vector of flags indicating which properties to
%                        calculate.
%                        1 = MSS diffusion analysis.
%                        2 = diffusion mode analysis.
%                        3 = direct track properties.
%                        4 = merging and splitting properties.
%                        To calculate all, enter [1 2 3 4].
%                        Optional. Default: [1 2 3 4].
%       windowNumbersAssignExtIn: Output of this code, from a previous run,
%                        in the case when an update is desired instead of a
%                        full calculate from scratch.
%
%OUTPUT windowNumbersAssignExt: 
%                          *** FIX ***
%                          5D array of dimensions (number of bands) x
%                          (number of slices) x (number of window frames-1)
%                          x (number of window frames-1) x 2 storing for
%                          each window in each frame the number of merges
%                          and splits that fall in it, not only in its
%                          proper frame range but throughout all frame
%                          ranges, as indicated by the 4th dimension. The
%                          5th dimension indicates merges (index 1) and
%                          splits (index 2).
%
%Khuloud Jaqaman, August 2012

%% Input

if nargin < 7
    disp('--assignNumbers2Windows: Incorrect number of input arguments!');
    return
end

if nargin < 8 || isempty(lengthMinMax)
    lengthMinMax = [5 99];
end

if nargin < 9 || isempty(prop2calc);
    prop2calc = [1 2 3 4];
end

if nargin < 10 || isempty(windowNumbersAssignExtIn)
    windowNumbersAssignExtIn = [];
end

%generate winFrameMin and winFramesMax for parfor loop
winFramesMin = winFrames(1:end-1);
winFramesMax = winFrames(2:end);

%% Pre-processing

%get number of frames that have windows and number of windows parallel to
%the edge
[numWinFramesM1,numWinPara] = size(winPositions);
numWinFramesM1 = numWinFramesM1 - 1;

%find number of windows perpendicular to the edge
nBands = cellfun(@(x)(numel(x)),winPositions);
numWinPerp = max(nBands(:));

%collect various particle motion information
[trackLft,trajClass,~,~,trajDiffMode,~,numDiffMode,~,~,~,~,paraProtDispTmp] = ...
    collectParticleBehavior(tracksFinal,diffAnalysisRes,diffModeAnRes,directTrackChar);

%for merging and splitting calculation, keep only tracks satisfying minimum
%length condition and containing more then 1 segment
%note that the maximum length condition is not imposed here because of the
%complication of compound tracks
criteria.lifeTime.min = lengthMinMax(1);
criteria.numSegments.min = 2;
indxGood = chooseTracks(tracksFinal,criteria);
tracksFinalMS = tracksFinal(indxGood);

%find merges and splits in compound tracks
[mergesInfo,splitsInfo,mergesInfoSpace,splitsInfoSpace] = ...
    findMergesSplits(tracksFinalMS,2,1,0,0);

%make a list of merges and splits, storing only their times and locations
mergesInfoTime = mergesInfo(:,4:end);
mergesInfoTime = mergesInfoTime(:);
mergesInfoSpaceX = mergesInfoSpace(:,1:2:end);
mergesInfoSpaceX = mergesInfoSpaceX(:);
mergesInfoSpaceY = mergesInfoSpace(:,2:2:end);
mergesInfoSpaceY = mergesInfoSpaceY(:);
indxKeep = find(mergesInfoTime~=0);
mergesInfoTime = mergesInfoTime(indxKeep);
mergesInfoSpaceX = mergesInfoSpaceX(indxKeep);
mergesInfoSpaceY = mergesInfoSpaceY(indxKeep);
splitsInfoTime = splitsInfo(:,4:end);
splitsInfoTime = splitsInfoTime(:);
splitsInfoSpaceX = splitsInfoSpace(:,1:2:end);
splitsInfoSpaceX = splitsInfoSpaceX(:);
splitsInfoSpaceY = splitsInfoSpace(:,2:2:end);
splitsInfoSpaceY = splitsInfoSpaceY(:);
indxKeep = find(splitsInfoTime~=0);
splitsInfoTime = splitsInfoTime(indxKeep);
splitsInfoSpaceX = splitsInfoSpaceX(indxKeep);
splitsInfoSpaceY = splitsInfoSpaceY(indxKeep);

%% Number of particles in various categories, and merges and splits

%group merges and splits based on what window frame range they fall into
msGroup = repmat(struct('indxM',[],'indxS',[]),numWinFramesM1,1);
parfor iWinFrame = 1 : numWinFramesM1
    
    %get current spt frame number and next spt frame number
    minFrame = winFramesMin(iWinFrame);
    maxFrame = winFramesMax(iWinFrame);
    
    %find merges and splits that happen in this frame range
    indxFrameRangeM = find(mergesInfoTime>=minFrame & mergesInfoTime<maxFrame);
    indxFrameRangeS = find(splitsInfoTime>=minFrame & splitsInfoTime<maxFrame);
    
    %store this information for later use
    msGroup(iWinFrame).indxM = indxFrameRangeM;
    msGroup(iWinFrame).indxS = indxFrameRangeS;
    
end

%initialize arrays
[numTotalPerWindow,numNetDispPosPerWindow,numNetDispNegPerWindow,...
    numUnclassPerWindow,numLinPerWindow,numIsoPerWindow,...
    numIsoUnclassPerWindow,numConfPerWindow,numBrownPerWindow,numDirPerWindow,...
    numMergePerWindow,numSplitPerWindow] = ...
    deal(NaN(numWinPerp,numWinPara,numWinFramesM1,numWinFramesM1));
numModePerWindow = NaN(numWinPerp,numWinPara,numWinFramesM1,...
    numDiffMode+1,numWinFramesM1);

parfor iWinFrameExt = 1 : numWinFramesM1 %window frames to fetch tracks
    
    %initialize temporary matrices
    [numTotalTmp,numNetDispPosTmp,numNetDispNegTmp,numUnclassTmp,numLinTmp,...
        numIsoTmp,numIsoUnclassTmp,numConfTmp,numBrownTmp,numDirTmp,...
        numMergeTmp,numSplitTmp] = deal(zeros(numWinPerp,numWinPara,numWinFramesM1));
    numModeTmp = zeros(numWinPerp,numWinPara,numWinFramesM1,numDiffMode+1);
    
    %find merges and splits in this time window
    indxFrameRangeM = msGroup(iWinFrameExt).indxM;
    indxFrameRangeS = msGroup(iWinFrameExt).indxS;
    
    %get the positions of these merges and splits
    xCoordMergeFR = mergesInfoSpaceX(indxFrameRangeM);
    yCoordMergeFR = mergesInfoSpaceY(indxFrameRangeM);
    xCoordSplitFR = splitsInfoSpaceX(indxFrameRangeS);
    yCoordSplitFR = splitsInfoSpaceY(indxFrameRangeS);

    for iWinFrame = 1 : numWinFramesM1 %window frames to fetch windows
        for iPara = 1 : numWinPara %slices
            for iPerp = 1 : nBands(iWinFrame,iPara) %bands
                
                %if the window in this frame has a finite size
                if ~isempty(winPositions{iWinFrame,iPara}{iPerp})
                    
                    %get current tracks
                    tracksCurrent = windowTrackAssignExt{iPerp,iPara,iWinFrame,iWinFrameExt};
                    
                    %keep only tracks whose length is within the required range
                    trackLftCurrent = trackLft(tracksCurrent);
                    tracksCurrent = tracksCurrent(trackLftCurrent>=lengthMinMax(1) & ...
                        trackLftCurrent<=lengthMinMax(2));
                    numTracksCurrent = length(tracksCurrent);
                    
                    if any(prop2calc==3)
                        
                        %total number of particles
                        numTotalTmp(iPerp,iPara,iWinFrame) = numTracksCurrent;
                        
                        %number of particles with net displacement parallel or
                        %anti-parallel to protrusion vector
                        paraProtDispCurrent = paraProtDispTmp(tracksCurrent,1);
                        tracksPos = tracksCurrent(paraProtDispCurrent >= 0);
                        tracksNeg = tracksCurrent(paraProtDispCurrent < 0);
                        numNetDispPosTmp(iPerp,iPara,iWinFrame) = length(tracksPos);
                        numNetDispNegTmp(iPerp,iPara,iWinFrame) = length(tracksNeg);
                        
                    end
                    
                    %asym+MSS analysis classification
                    if any(prop2calc==1)
                        trajClassCurrent = trajClass(tracksCurrent);
                        n = hist(trajClassCurrent,1:5);
                        numUnclassTmp(iPerp,iPara,iWinFrame) = numTracksCurrent - sum(n);
                        numLinTmp(iPerp,iPara,iWinFrame) = n(5);
                        numIsoTmp(iPerp,iPara,iWinFrame) = sum(n(1:4));
                        numIsoUnclassTmp(iPerp,iPara,iWinFrame) = n(4);
                        numConfTmp(iPerp,iPara,iWinFrame) = n(1);
                        numBrownTmp(iPerp,iPara,iWinFrame) = n(2);
                        numDirTmp(iPerp,iPara,iWinFrame) = n(3);
                    end
                    
                    %mode analysis classification
                    if any(prop2calc==2)
                        trajModeCurrent = trajDiffMode(tracksCurrent);
                        trajModeCurrent(isnan(trajModeCurrent)) = numDiffMode+1;
                        numModeTmp(iPerp,iPara,iWinFrame,:) = hist(trajModeCurrent,1:numDiffMode+1);
                    end
                    
                    %get the window boundaries
                    if any(prop2calc==4)
                        windowsPoly = [winPositions{iWinFrame,iPara}{iPerp}{:}];
                        winX = windowsPoly(1,:);
                        winY = windowsPoly(2,:);
                        
                        %find the merges whose position lies in this window
                        indxWin = inpolygon(xCoordMergeFR,yCoordMergeFR,winX,winY);
                        numMergeTmp(iPerp,iPara,iWinFrame) = length(find(indxWin));
                        
                        %find the splits whose position lies in this window
                        indxWin = inpolygon(xCoordSplitFR,yCoordSplitFR,winX,winY);
                        numSplitTmp(iPerp,iPara,iWinFrame) = length(find(indxWin));
                    end

                end %(if ~isempty(winPositions{iWinFrame,iPara}{iPerp}))
                
            end
        end
    end
    
    %store extended assignment in overall matrix
    numTotalPerWindow(:,:,:,iWinFrameExt) = numTotalTmp;
    numNetDispPosPerWindow(:,:,:,iWinFrameExt) = numNetDispPosTmp;
    numNetDispNegPerWindow(:,:,:,iWinFrameExt) = numNetDispNegTmp;
    numUnclassPerWindow(:,:,:,iWinFrameExt) = numUnclassTmp;
    numLinPerWindow(:,:,:,iWinFrameExt) = numLinTmp;
    numIsoPerWindow(:,:,:,iWinFrameExt) = numIsoTmp;
    numIsoUnclassPerWindow(:,:,:,iWinFrameExt) = numIsoUnclassTmp;
    numConfPerWindow(:,:,:,iWinFrameExt) = numConfTmp;
    numBrownPerWindow(:,:,:,iWinFrameExt) = numBrownTmp;
    numDirPerWindow(:,:,:,iWinFrameExt) = numDirTmp;
    numModePerWindow(:,:,:,:,iWinFrameExt) = numModeTmp;
    numMergePerWindow(:,:,:,iWinFrameExt) = numMergeTmp;
    numSplitPerWindow(:,:,:,iWinFrameExt) = numSplitTmp;
    
end %(parfor iWinFrameExt = 1 : numWinFramesM1)

%put in structure for output
if isempty(windowNumbersAssignExtIn) %make brand new structure
    
    windowNumbersAssignExt = struct('Particles',numTotalPerWindow,...
        'NetDispPos',numNetDispPosPerWindow,'NetDispNeg',numNetDispNegPerWindow,...
        'Unclass',numUnclassPerWindow,'Lin',numLinPerWindow,'Iso',numIsoPerWindow,...
        'IsoUnclass',numIsoUnclassPerWindow,'Conf',numConfPerWindow,...
        'Brown',numBrownPerWindow,'Dir',numDirPerWindow,'ModeClass',numModePerWindow,...
        'Merge',numMergePerWindow,'Split',numSplitPerWindow);
    
else %update input structure
    
    windowNumbersAssignExt = windowNumbersAssignExtIn;
    if any(prop2calc==1)
        windowNumbersAssignExt.Unclass = numUnclassPerWindow;
        windowNumbersAssignExt.Lin = numLinPerWindow;
        windowNumbersAssignExt.Iso = numIsoPerWindow;
        windowNumbersAssignExt.IsoUnclass = numIsoUnclassPerWindow;
        windowNumbersAssignExt.Conf = numConfPerWindow;
        windowNumbersAssignExt.Brown = numBrownPerWindow;
        windowNumbersAssignExt.Dir = numDirPerWindow;
    end
    if any(prop2calc==2)
        windowNumbersAssignExt.ModeClass = numModePerWindow;
    end
    if any(prop2calc==3)
        windowNumbersAssignExt.Particles = numTotalPerWindow;
        windowNumbersAssignExt.NetDispPos = numNetDispPosPerWindow;
        windowNumbersAssignExt.NetDispNeg = numNetDispNegPerWindow;
    end
    if any(prop2calc==4)
        windowNumbersAssignExt.Merge = numMergePerWindow;
        windowNumbersAssignExt.Split = numSplitPerWindow;
    end
    
end

%% ~~~ the end ~~~

