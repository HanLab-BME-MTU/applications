function [sptPropInWindow,windowTrackAssign,trackWindowAssign,windowSize] = ...
    particleBehaviorInWindows(tracksFinal,winPositions,winFrames,...
    protSamples,diffAnalysisRes,minLength)
%PARTICLEBEHAVIORINWINDOWS averages single particle behavior in windows based on cell edge segmentation
%
%SYNOPSIS [sptPropInWindow,windowTrackAssign,trackWindowAssign,windowSize] = ...
%    particleBehaviorInWindows(tracksFinal,winPositions,winFrames,...
%    protSamples,diffAnalysisRes,minLength)
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
%       protPerWindow  : The protrusion samples as output by the windowing
%                        software.
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%                        Optional. If not input but needed, it will be
%                        calculated within the code.
%       minLength      : Minimum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: 5.
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
%Khuloud Jaqaman, May 2010

%% Input

if nargin < 4
    disp('--particleBehaviorInWindows: Missing input arguments!');
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

%get normalized protrusion vectors
protVec = protSamples.avgVector;
protVecMag = sqrt(sum(protVec.^2,3));
protVecUnit = protVec ./ repmat(protVecMag,[1 1 2]);

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
ratioDispDirTmp = paraDirDispTmp ./ perpDirDispTmp;
ratioDispProtTmp = paraProtDispTmp ./ perpProtDispTmp;

%% Calculate property values per window

%initialize output variables
[spDensity,f2fDispMag2D,angleMean,angleStd,...
    f2fDispMagParaDir,f2fDispMagPerpDir,f2fDispSignParaDir,f2fDispSignPerpDir,...
    f2fDispMagParaProt,f2fDispMagPerpProt,f2fDispSignParaProt,f2fDispSignPerpProt,...
    ratioDispMagDir,ratioDispSignDir,ratioDispMagProt,ratioDispSignProt,...
    asymParam,f2fDispMag2DLin,angleMeanLin,angleStdLin,f2fDispMagParaDirLin,...
    f2fDispMagPerpDirLin,f2fDispSignParaDirLin,f2fDispSignPerpDirLin,...
    f2fDispMagParaProtLin,f2fDispMagPerpProtLin,f2fDispSignParaProtLin,...
    f2fDispSignPerpProtLin,ratioDispMagDirLin,ratioDispSignDirLin,...
    ratioDispMagProtLin,ratioDispSignProtLin,asymParamLin,fracUnclass,...
    fracLin,fracIso,fracIsoUnclass,fracConf,fracBrown,fracDir,...
    diffCoef,confRad] = ...
    deal(struct('values',NaN(numWinPerp,numWinPara,numWinFrames-1),...
    'numPoints',zeros(numWinPerp,numWinPara,numWinFrames-1)));

%go over all windows and calculate particle properties in each
for iFrame = 1 : numWinFrames-1
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %get the tracks belonging to this window
            tracksCurrent = windowTrackAssign{iPerp,iPara,iFrame};
            numTracksCurrent = length(tracksCurrent);
            
            %if there are tracks in this window ...
            if numTracksCurrent ~= 0
                
                %From the tracks directly ...
                
                %calculate the average particle density
                spDensity.values(iPerp,iPara,iFrame) = sum(trackLft(tracksCurrent)) / ...
                    winSize(iPerp,iPara,iFrame) / numSPTFrames;
                spDensity.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %calculate the mean and std of angle between direction of
                %motion and protrusion vector
                angleMean.values(iPerp,iPara,iFrame) = nanmean(angleWithProtTmp(tracksCurrent));
                angleMean.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                if numTracksCurrent > 1
                    angleStd.values(iPerp,iPara,iFrame) = nanstd(angleWithProtTmp(tracksCurrent));
                end
                angleStd.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %calculate the average frame-to-frame displacement
                f2fDispMag2D.values(iPerp,iPara,iFrame) = nanmean(f2fDispTmp(tracksCurrent));
                f2fDispMag2D.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %calculate the average frame-to-frame displacements along
                %and perpendicular to the direction of motion
                f2fDispSignParaDir.values(iPerp,iPara,iFrame) = nanmean(paraDirDispTmp(tracksCurrent,1));
                f2fDispSignParaDir.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                f2fDispSignPerpDir.values(iPerp,iPara,iFrame) = nanmean(perpDirDispTmp(tracksCurrent,1));
                f2fDispSignPerpDir.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                f2fDispMagParaDir.values(iPerp,iPara,iFrame) = nanmean(paraDirDispTmp(tracksCurrent,2));
                f2fDispMagParaDir.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                f2fDispMagPerpDir.values(iPerp,iPara,iFrame) = nanmean(perpDirDispTmp(tracksCurrent,2));
                f2fDispMagPerpDir.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %calculate the average frame-to-frame displacements along
                %and perpendicular to the protrusion vector
                f2fDispSignParaProt.values(iPerp,iPara,iFrame) = nanmean(paraProtDispTmp(tracksCurrent,1));
                f2fDispSignParaProt.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                f2fDispSignPerpProt.values(iPerp,iPara,iFrame) = nanmean(perpProtDispTmp(tracksCurrent,1));
                f2fDispSignPerpProt.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                f2fDispMagParaProt.values(iPerp,iPara,iFrame) = nanmean(paraProtDispTmp(tracksCurrent,2));
                f2fDispMagParaProt.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                f2fDispMagPerpProt.values(iPerp,iPara,iFrame) = nanmean(perpProtDispTmp(tracksCurrent,2));
                f2fDispMagPerpProt.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %calculate the ratio of perpendicular to parallel
                %displacement components
                ratioDispSignDir.values(iPerp,iPara,iFrame) = nanmean(ratioDispDirTmp(tracksCurrent,1));
                ratioDispSignDir.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                ratioDispMagDir.values(iPerp,iPara,iFrame)  = nanmean(ratioDispDirTmp(tracksCurrent,2));
                ratioDispMagDir.numPoints(iPerp,iPara,iFrame)  = numTracksCurrent;
                ratioDispSignProt.values(iPerp,iPara,iFrame) = nanmean(ratioDispProtTmp(tracksCurrent,1));
                ratioDispSignProt.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                ratioDispMagProt.values(iPerp,iPara,iFrame)  = nanmean(ratioDispProtTmp(tracksCurrent,2));
                ratioDispMagProt.numPoints(iPerp,iPara,iFrame)  = numTracksCurrent;
                
                %calculate the asymmetry parameter (ratio of maximum to
                %minimum eigenvalues of position covariance matrix)
                asymParam.values(iPerp,iPara,iFrame) = nanmean(asymParamTmp(tracksCurrent));
                asymParam.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %repeat the above but only for tracks classified as
                %asymmetric
                trajClassCurrent = trajClass(tracksCurrent);
                tracksLin = tracksCurrent(trajClassCurrent==5);
                numTracksLin = length(tracksLin);
                if numTracksLin ~= 0
                    angleMeanLin.values(iPerp,iPara,iFrame) = nanmean(angleWithProtTmp(tracksLin));
                    angleMeanLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    if numTracksLin > 1
                        angleStdLin.values(iPerp,iPara,iFrame) = nanstd(angleWithProtTmp(tracksLin));
                    end
                    angleStdLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispMag2DLin.values(iPerp,iPara,iFrame) = nanmean(f2fDispTmp(tracksLin));
                    f2fDispMag2DLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispSignParaDirLin.values(iPerp,iPara,iFrame) = nanmean(paraDirDispTmp(tracksLin,1));
                    f2fDispSignParaDirLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispSignPerpDirLin.values(iPerp,iPara,iFrame) = nanmean(perpDirDispTmp(tracksLin,1));
                    f2fDispSignPerpDirLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispMagParaDirLin.values(iPerp,iPara,iFrame) = nanmean(paraDirDispTmp(tracksLin,2));
                    f2fDispMagParaDirLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispMagPerpDirLin.values(iPerp,iPara,iFrame) = nanmean(perpDirDispTmp(tracksLin,2));
                    f2fDispMagPerpDirLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispSignParaProtLin.values(iPerp,iPara,iFrame) = nanmean(paraProtDispTmp(tracksLin,1));
                    f2fDispSignParaProtLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispSignPerpProtLin.values(iPerp,iPara,iFrame) = nanmean(perpProtDispTmp(tracksLin,1));
                    f2fDispSignPerpProtLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispMagParaProtLin.values(iPerp,iPara,iFrame) = nanmean(paraProtDispTmp(tracksLin,2));
                    f2fDispMagParaProtLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    f2fDispMagPerpProtLin.values(iPerp,iPara,iFrame) = nanmean(perpProtDispTmp(tracksLin,2));
                    f2fDispMagPerpProtLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    ratioDispSignDirLin.values(iPerp,iPara,iFrame) = nanmean(ratioDispDirTmp(tracksLin,1));
                    ratioDispSignDirLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    ratioDispMagDirLin.values(iPerp,iPara,iFrame)  = nanmean(ratioDispDirTmp(tracksLin,2));
                    ratioDispMagDirLin.numPoints(iPerp,iPara,iFrame)  = numTracksLin;
                    ratioDispSignProtLin.values(iPerp,iPara,iFrame) = nanmean(ratioDispProtTmp(tracksLin,1));
                    ratioDispSignProtLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                    ratioDispMagProtLin.values(iPerp,iPara,iFrame)  = nanmean(ratioDispProtTmp(tracksLin,2));
                    ratioDispMagProtLin.numPoints(iPerp,iPara,iFrame)  = numTracksLin;
                    asymParamLin.values(iPerp,iPara,iFrame) = nanmean(asymParamTmp(tracksLin));
                    asymParamLin.numPoints(iPerp,iPara,iFrame) = numTracksLin;
                end
                
                %From the asymmetry and diffusion analysis ...
                
                %calculate the fraction of tracks in each motion category
                
                %first completely unclassified tracks
                trajClassCurrent = trajClass(tracksCurrent);
                fracUnclass.values(iPerp,iPara,iFrame) = ...
                    length(find(isnan(trajClassCurrent)))/numTracksCurrent;
                fracUnclass.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %then tracks with asymmetry analysis
                trajClassCurrent = trajClassCurrent(~isnan(trajClassCurrent));
                numTracksAsym = length(trajClassCurrent);
                if numTracksAsym > 0
                    
                    %tracks classified as linear or isotropic
                    fracLin.values(iPerp,iPara,iFrame) = ...
                        length(find(trajClassCurrent==5))/numTracksAsym;
                    fracLin.numPoints(iPerp,iPara,iFrame) = numTracksAsym;
                    fracIso.values(iPerp,iPara,iFrame) = ...
                        length(find(trajClassCurrent~=5))/numTracksAsym;
                    fracIso.numPoints(iPerp,iPara,iFrame) = numTracksAsym;
                    
                    %within the tracks classified as isotropic
                    trajClassCurrent = trajClassCurrent(trajClassCurrent~=5);
                    numTracksIso = length(trajClassCurrent);
                    if numTracksIso > 0
                        
                        %tracks without diffusion analysis
                        fracIsoUnclass.values(iPerp,iPara,iFrame) = ...
                            length(find(trajClassCurrent==4))/numTracksIso;
                        fracIsoUnclass.numPoints(iPerp,iPara,iFrame) = numTracksIso;
                        
                        %finally tracks with diffusion analysis
                        trajClassCurrent = trajClassCurrent(trajClassCurrent~=4);
                        numTracksDiff = length(trajClassCurrent);
                        if numTracksDiff > 0
                            
                            fracClass = hist(trajClassCurrent,(1:3))/numTracksDiff;
                            fracConf.values(iPerp,iPara,iFrame) = fracClass(1);
                            fracConf.numPoints(iPerp,iPara,iFrame) = numTracksDiff;
                            fracBrown.values(iPerp,iPara,iFrame) = fracClass(2);
                            fracBrown.numPoints(iPerp,iPara,iFrame) = numTracksDiff;
                            fracDir.values(iPerp,iPara,iFrame) = fracClass(3);
                            fracDir.numPoints(iPerp,iPara,iFrame) = numTracksDiff;
                            
                        end %(if numTracksDiff > 0)
                        
                    end %(if numTracksIso > 0)
                    
                end %(if numTracksAsym > 0)
                
                %calculate the average diffusion coefficient
                tmpVec = diffCoefGen(tracksCurrent);
                tmpVec = tmpVec(~isnan(tmpVec));
                diffCoef.values(iPerp,iPara,iFrame) = mean(tmpVec);
                diffCoef.numPoints(iPerp,iPara,iFrame) = length(tmpVec);
                
                %calculate the average confinement radius
                tmpVec = confRadAll(tracksCurrent);
                tmpVec = tmpVec(~isnan(tmpVec));
                confRad.values(iPerp,iPara,iFrame) = mean(tmpVec);
                confRad.numPoints(iPerp,iPara,iFrame) = length(tmpVec);

            end %(if numTracksCurrent ~= 0)
            
        end %(for iPerp = 1 : numWinPerp)
    end %(for iPara = 1 : numWinPara)
end %(for iFrame = 1 : numWinFrames-1)

%store all properties in output structure
sptPropInWindow = struct('spDensity',spDensity,'f2fDispMag2D',f2fDispMag2D,...
    'angleMean',angleMean,'angleStd',angleStd,...
    'f2fDispMagParaDir',f2fDispMagParaDir,'f2fDispMagPerpDir',f2fDispMagPerpDir,...
    'f2fDispSignParaDir',f2fDispSignParaDir,'f2fDispSignPerpDir',f2fDispSignPerpDir,...
    'f2fDispMagParaProt',f2fDispMagParaProt,'f2fDispMagPerpProt',f2fDispMagPerpProt,...
    'f2fDispSignParaProt',f2fDispSignParaProt,'f2fDispSignPerpProt',f2fDispSignPerpProt,...
    'ratioDispMagDir',ratioDispMagDir,'ratioDispSignDir',ratioDispSignDir,...
    'ratioDispMagProt',ratioDispMagProt,'ratioDispSignProt',ratioDispSignProt,...    
    'asymParam',asymParam,'f2fDispMag2DLin',f2fDispMag2DLin,...
    'angleMeanLin',angleMeanLin,'angleStdLin',angleStdLin,...
    'f2fDispMagParaDirLin',f2fDispMagParaDirLin,'f2fDispMagPerpDirLin',f2fDispMagPerpDirLin,...
    'f2fDispSignParaDirLin',f2fDispSignParaDirLin,'f2fDispSignPerpDirLin',f2fDispSignPerpDirLin,...
    'f2fDispMagParaProtLin',f2fDispMagParaProtLin,'f2fDispMagPerpProtLin',f2fDispMagPerpProtLin,...
    'f2fDispSignParaProtLin',f2fDispSignParaProtLin,'f2fDispSignPerpProtLin',f2fDispSignPerpProtLin,...
    'ratioDispMagDirLin',ratioDispMagDirLin,'ratioDispSignDirLin',ratioDispSignDirLin,...
    'ratioDispMagProtLin',ratioDispMagProtLin,'ratioDispSignProtLin',ratioDispSignProtLin,...    
    'asymParamLin',asymParamLin,'fracUnclass',fracUnclass,'fracLin',fracLin,...
    'fracIso',fracIso,'fracIsoUnclass',fracIsoUnclass,'fracConf',fracConf,...
    'fracBrown',fracBrown,'fracDir',fracDir,'diffCoef',diffCoef,'confRad',confRad);

%% ~~~ the end ~~~

    %     %extract the x- and y-coordinates from the track matrix
    %     xCoord = tracksFinal(:,1:8:end);
    %     yCoord = tracksFinal(:,2:8:end);
    %
    %     %calculate the average frame-to-frame displacement per track
    %     frame2frameDisp = nanmean( sqrt( diff(xCoord,[],2).^2 + diff(yCoord,[],2).^2 ) ,2);
    %
    %     %calculate direction of motion using the position covariance matrix
    %     %also calculate frame-to-frame displacement along direction of motion
    %     motionDirection   = NaN(numTracks,2);
    %     parallelDisp      = NaN(numTracks,1);
    %     perpendicularDisp = NaN(numTracks,1);
    %     for iTrack = 1 : numTracks
    %         posCov = nancov([xCoord(iTrack,:); yCoord(iTrack,:)]');
    %         [eigVec,eigVal] = eig(posCov);
    %         eigVal = diag(eigVal);
    %         eigValMax = max(eigVal);
    %         eigValMin = min(eigVal);
    %         motionDirection(iTrack,:) = eigVec(:,eigVal==eigValMax)';
    %         parallelDisp(iTrack)      = sqrt( eigValMax / (trackLft(iTrack)-1) );
    %         perpendicularDisp(iTrack) = sqrt( eigValMin / (trackLft(iTrack)-1) );
    %     end
    
