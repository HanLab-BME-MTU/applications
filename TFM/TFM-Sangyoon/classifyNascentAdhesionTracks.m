function [idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9]=classifyNascentAdhesionTracks(pathForColocalization,varargin)
% 
% this function classifyNascentAdhesionTracks(pathForColocalization,varargin) reads from Colocalization folder 
% and classify NA tracks based on fluorescence signal amplitude evolution
% and distance to the edge. 
% group 1: NAs that form at the edge but say there as the edge protrude and
% turn over:
% group 2: NAs that form and mature as the edge protrude: 
% group 3: NAs that move along the edge and quickly turn over

%% input reading
ip =inputParser;
ip.addRequired('pathForColocalization',@ischar)
ip.addParamValue('tracksNA',[],@isstruct); % selcted track ids
ip.parse(pathForColocalization,varargin{:});
pathForColocalization=ip.Results.pathForColocalization;
tracksNA=ip.Results.tracksNA;

%% Load processed data
disp('Loading raw files ...')
tic
if isempty(tracksNA)
    tracksNA = load([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNA = tracksNA.tracksNA;
end
imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
imgMap = imgMap.paxImgStack;
tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
tMap = tMap.tMap;
outputPath = [pathForColocalization filesep 'trackAnalysis'];
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end
numFrames = size(imgMap,3);
startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
% movieData to find out pixel size
coloPath = fileparts(pathForColocalization);
MDPath = fileparts(coloPath);
MDfilePath = [MDPath filesep 'movieData.mat'];
MD = load(MDfilePath,'MD');
MD = MD.MD;
pixSize = MD.pixelSize_; % nm/pixel
tInterval = MD.timeInterval_; % time interval in sec
scaleBar = 1; %micron
toc
%% Look at distToEdge and ampTotal
% How to find a threshold for distToEdge?
% % Plot histogram
tIntervalMin = tInterval/60; % in min
periodMin = 1;
periodFrames = floor(periodMin/tIntervalMin); % early period in frames
% disp(['Extrapolate by ' num2str(periodMin) ' min...'])
% tic
% for k=1:numel(tracksNA)
%     tracksNA(k) = readIntensityFromTracks(tracksNA(k),imgMap,1,'extraLength',10);
%     tracksNA(k) = readIntensityFromTracks(tracksNA(k),tMap,2,'extraLength',10);
% end
% toc
%% regression for slopes
% disp('Regression...')
% tic
% for k=1:numel(tracksNA)
%     curEndFrame = min(tracksNA(k).startingFrameExtra+periodFrames-1,tracksNA(k).endingFrame);
%     curEarlyPeriod = curEndFrame - tracksNA(k).startingFrameExtra+1;
%     [~,curM] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:curEndFrame));
%     tracksNA(k).earlyAmpSlope = curM; % in a.u./min
% 
%     curStartFrame = max(tracksNA(k).startingFrame,tracksNA(k).endingFrameExtra-periodFrames+1);
%     curLatePeriod = tracksNA(k).endingFrameExtra - curStartFrame+1;
%     [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),tracksNA(k).ampTotal(curStartFrame:tracksNA(k).endingFrameExtra));
%     tracksNA(k).lateAmpSlope = curMlate; % in a.u./min
% 
%     curEndFrame = min(tracksNA(k).startingFrame+periodFrames-1,tracksNA(k).endingFrame);
%     curEarlyPeriod = curEndFrame - tracksNA(k).startingFrame+1;
%     [~,curMdist] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).distToEdge(tracksNA(k).startingFrame:curEndFrame));
%     tracksNA(k).distToEdgeSlope = curMdist; % in a.u./min
%     tracksNA(k).distToEdgeChange = (tracksNA(k).distToEdge(end)-tracksNA(k).distToEdge(tracksNA(k).startingFrame)); % in pixel
%     sF=tracksNA(k).startingFrame;
%     eF=tracksNA(k).endingFrame;
%     % Determining protrusion/retraction based on closestBdPoint and [xCoord
%     % yCoord]. If the edge at one time point does not cross the adhesion
%     % tracks over entire life time, there is no problem. But that's not
%     % always the case: adhesion tracks crosses the cell boundary at first
%     % or last time point. And we don't know if the adhesion tracks are
%     % moving in the direction of protrusion or retraction. Thus, what I
%     % will do is to calculate the inner product of vectors from the
%     % adhesion tracks from the closestBdPoint at the first frame. If the
%     % product is positive, it means both adhesion points are in the same
%     % side. And if distance from the boundary to the last point is larger
%     % than the distance to the first point, it means the adhesion track is
%     % retracting. In the opposite case, the track is protruding (but it
%     % will happen less likely because usually the track would cross the
%     % first frame boundary if it is in the protrusion direction). 
%     
%     % We need to find out boundary points projected on the line of adhesion track 
%     edge = [tracksNA(k).xCoord(sF) tracksNA(k).yCoord(sF) tracksNA(k).xCoord(eF) tracksNA(k).yCoord(eF)];
%     line = edgeToLine(edge);
%     firstBdPoint = [tracksNA(k).closestBdPoint(sF,1) tracksNA(k).closestBdPoint(sF,2)];
%     firstBdPointProjected = projPointOnLine(firstBdPoint, line);
%     lastBdPoint = [tracksNA(k).closestBdPoint(eF,1) tracksNA(k).closestBdPoint(eF,2)];
%     lastBdPointProjected = projPointOnLine(lastBdPoint, line);
%     
%     firstAdhToFirstBdPoint = [tracksNA(k).xCoord(sF)-firstBdPointProjected(1), tracksNA(k).yCoord(sF)-firstBdPointProjected(2)];
%     lastAdhToFirstBdPoint = [tracksNA(k).xCoord(eF)-firstBdPointProjected(1), tracksNA(k).yCoord(eF)-firstBdPointProjected(2)];
%     firstAdhToLastBdPoint = [tracksNA(k).xCoord(sF)-lastBdPointProjected(1), tracksNA(k).yCoord(sF)-lastBdPointProjected(2)];
%     lastAdhToLastBdPoint = [tracksNA(k).xCoord(eF)-lastBdPointProjected(1), tracksNA(k).yCoord(eF)-lastBdPointProjected(2)];
%     if firstAdhToFirstBdPoint*lastAdhToFirstBdPoint'>0 % both adhesion points are in the same side
%         tracksNA(k).advanceDist = (firstAdhToFirstBdPoint(1)^2 + firstAdhToFirstBdPoint(2)^2)^0.5 - ...
%                                                         (lastAdhToFirstBdPoint(1)^2 + lastAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%         tracksNA(k).edgeAdvanceDist = (lastAdhToLastBdPoint(1)^2 + lastAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                         (lastAdhToFirstBdPoint(1)^2 + lastAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%     else
%         if firstAdhToLastBdPoint*lastAdhToLastBdPoint'>0 % both adhesion points are in the same side w.r.t. last boundary point
%             tracksNA(k).advanceDist = (firstAdhToLastBdPoint(1)^2 + firstAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                             (lastAdhToLastBdPoint(1)^2 + lastAdhToLastBdPoint(2)^2)^0.5; % in pixel
%             tracksNA(k).edgeAdvanceDist = (firstAdhToLastBdPoint(1)^2 + firstAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                             (firstAdhToFirstBdPoint(1)^2 + firstAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%         else
%             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. This does not make sense!'])
%         end
%     end
%         
% %     tracksNA(k).advanceDist = ((tracksNA(k).xCoord(tracksNA(k).endingFrame)-tracksNA(k).xCoord(tracksNA(k).startingFrame))^2 +...
% %                                             (tracksNA(k).xCoord(tracksNA(k).endingFrame)-tracksNA(k).xCoord(tracksNA(k).startingFrame))^2)^0.5; % in pixel
% end
%% Let's start fresh. Let's use adhesion tracks and edge protrusion
% NAs in group 1 should experience pretty static adhesion movement while
% significant edge protrusion

% NAs in group 2 slides in retrograde direction and edge protrude. How can
% I know if the adhesion goes front or back? It can be decided by edge
% normal direction. Actually I can just use tracksNA.state

% NAs in group 3 roll with edge protrusion, so constant distace between
% adhesion and edge is expected

% % There should be multiple steps for outlier detection and exclusion
% % 1. Adhesions that has small adhesion density should be thrown away
% tracksNA = calculateAdhesionSpatialDensity(tracksNA,35); % This need to be more smart if time allows..
% numNeisAll = arrayfun(@(x) x.numNeighbors,tracksNA);
% % figure, histogram(numNeisAll,100)
% thresNumNeis = multithresh(numNeisAll,2);
% idxNAsHiDensity = numNeisAll>thresNumNeis(1);
% 
% % 2. tracks whose amplitudes are decreasing in the early phase (first one
% % minute)
% initialAmpSlope = arrayfun(@(x) x.earlyAmpSlope,tracksNA);
% 
% % set critical slope as local minimum in the negative value closest to zero
% [Namp,edgesHistAmp]= histcounts(initialAmpSlope);
% idxUnderZeroAmpSlope = edgesHistAmp<0;
% idxLocMinAmpSlope = false(size(idxUnderZeroAmpSlope));
% LocMinAmpSlope = locmin1d(Namp);
% idxLocMinAmpSlope(LocMinAmpSlope)=true;
% idxLocMinSlopeAmpClosestToZero = find(idxLocMinAmpSlope & idxUnderZeroAmpSlope, 1, 'last');
% critAmpSlope = edgesHistAmp(idxLocMinSlopeAmpClosestToZero)/3;
% % critAmpSlope = 0;
% idxNAsIncreasingAmp = initialAmpSlope>critAmpSlope;
% 
% % 3.tracks whose closest edge is not moving over the period of the track
% % make 0 to NaN for closesBdPoint
% for ii=1:numel(tracksNA)
%     idxZeros = tracksNA(ii).closestBdPoint(:)==0;
%     tracksNA(ii).closestBdPoint(idxZeros)=NaN;
% end
% % edgeMSD = arrayfun(@(x) nanmean((x.closestBdPoint(:,1)'-nanmean(x.closestBdPoint(:,1))).^2+...
% %     (x.closestBdPoint(:,2)'-nanmean(x.closestBdPoint(:,2))).^2),tracksNA(idxNAsHiDensity&idxNAsIncreasingAmp));
% edgeMSD = arrayfun(@(x) nanmean((x.closestBdPoint(:,1)'-nanmean(x.closestBdPoint(:,1))).^2+...
%     (x.closestBdPoint(:,2)'-nanmean(x.closestBdPoint(:,2))).^2),tracksNA);
% % figure, histogram(edgeMSD,1000)
% thresEdgeMSD = 2;%multithresh(edgeMSD,2);
% idxNAsEdgeMSD = edgeMSD>thresEdgeMSD;
% 
% % 4. Filter out NAs that formed inside the edge
% distToEdgeStart = arrayfun(@(x) x.distToEdge(x.startingFrame),tracksNA);
% thresDistToEdge = multithresh(distToEdgeStart,1);
% idxNAsFormedInside = arrayfun(@(x) x.distToEdge(x.startingFrame)<thresDistToEdge,tracksNA);
% % tracksNA = tracksNA(idxNAsFormedInside);
% 
% % 0. Adhesions whose lifetime is less than 20 sec
% thresMinLT = 20/tInterval;
% ltNAs = arrayfun(@(x) x.lifeTime,tracksNA); 
% idxNAsLongLT = ltNAs>thresMinLT;
% 
% 
% % figure, imshow(imgMap(:,:,end),[])
% % hold on
% % arrayfun(@(x) plot(x.xCoord,x.yCoord,'g'),tracksNA);
% % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'go'),tracksNA);
% % arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),tracksNA(idxNAsHiDensity));
% % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),tracksNA(idxNAsHiDensity));
% % arrayfun(@(x) plot(x.xCoord,x.yCoord,'c'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity));
% % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'co'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity));
% % arrayfun(@(x) plot(x.xCoord,x.yCoord,'m'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity &idxNAsEdgeMSD));
% % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'mo'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity & idxNAsEdgeMSD));
% 
% % Now filter out all the outliers
% idValidNAs = idxNAsIncreasingAmp & idxNAsHiDensity & idxNAsEdgeMSD & idxNAsFormedInside & idxNAsLongLT;
% tracksNA = tracksNA(idValidNAs);
% idxValidNAs = find(idValidNAs);

% Then, separate group 1 and 3 by dist to Edge at the ending frame
% group 3
% idxAdvancingNAs =  arrayfun(@(x) x.distToEdge(x.endingFrame)<thresDistToEdge,tracksNA);%arrayfun(@(x) x.edgeAdvanceDist<thresEdgeAdv,tracksNA);
% idxMatureNAs = idxMatureNAs & idxNAsFormedInside; % refinement of group 2 by eliminating ones at the edge
% group 1: static
% figure, histogram(asymTracks,100)
% figure, histogram(advanceDistNAs,200)
% thresMSD = multithresh(MSD,3);

% train data labeling
outputFile=[pathForColocalization filesep 'data' filesep 'selectedGroups.mat'];
reuseSelectedGroups = 'n';
% See if you can use existing tracks
if exist(outputFile,'file')
    reuseSelectedGroups=input('selectedGroups.mat is found. Do you want to use it (y), add some more on top of it (a), or start over(n or enter)?: ','s');
    if strcmp(reuseSelectedGroups, 'y') || strcmp(reuseSelectedGroups, 'a')
        idGroups = load(outputFile);
        idGroup1Selected = idGroups.idGroup1Selected;
        idGroup2Selected = idGroups.idGroup2Selected;
        idGroup3Selected = idGroups.idGroup3Selected;
        idGroup4Selected = idGroups.idGroup4Selected;
        idGroup5Selected = idGroups.idGroup5Selected;
        idGroup6Selected = idGroups.idGroup6Selected;
        idGroup7Selected = idGroups.idGroup7Selected;
        idGroup8Selected = idGroups.idGroup8Selected;
        idGroup9Selected = idGroups.idGroup9Selected;
%         newTracksNA=tracksNA(~idxMatureNAs);
        if strcmp(reuseSelectedGroups, 'a')
            % save to idTracks and iGroup
            idTracks = [idGroup1Selected idGroup2Selected idGroup3Selected idGroup4Selected idGroup5Selected ...
                idGroup6Selected idGroup7Selected idGroup8Selected idGroup9Selected];
            iGroup = [ones(size(idGroup1Selected))*1 ones(size(idGroup2Selected))*2 ones(size(idGroup3Selected))*3 ones(size(idGroup4Selected))*4 ones(size(idGroup5Selected))*5 ...
                ones(size(idGroup6Selected))*6 ones(size(idGroup7Selected))*7 ones(size(idGroup8Selected))*8 ones(size(idGroup9Selected))*9];
        end
    end
end
if ~exist(outputFile,'file') || strcmp(reuseSelectedGroups, 'n') || isempty(reuseSelectedGroups) || strcmp(reuseSelectedGroups, 'a')
    if strcmp(reuseSelectedGroups, 'a')
        display('Click additional tracks that belong to each group ...')
        [idTracksAdditional, iGroupAdditional] = showAdhesionTracks(pathForColocalization,'all','tracksNA',tracksNA);
        idTracks = [idTracks idTracksAdditional];
        iGroup = [iGroup iGroupAdditional];
    else
        display('Click tracks that belong to each group ...')
    %     newTracksNA=tracksNA(~idxMatureNAs);
    %     idNAs = find(~idxMatureNAs);
        [idTracks, iGroup] = showAdhesionTracks(pathForColocalization,'all','tracksNA',tracksNA);
    end
    % for duplicate ids that were assigned to multiple different groups,
    % assign them to the latest group
    [uniqIdTracks, ia, ic]=unique(idTracks,'stable');
    uniqIGroup = zeros(size(uniqIdTracks));
    for jj=1:numel(ia)
        laterIdx = find(ic==jj,1,'last');
        uniqIGroup(jj)=iGroup(laterIdx);
    end
    idGroup1Selected = uniqIdTracks(uniqIGroup==1);
    idGroup2Selected = uniqIdTracks(uniqIGroup==2);
    idGroup3Selected = uniqIdTracks(uniqIGroup==3);
    idGroup4Selected = uniqIdTracks(uniqIGroup==4);
    idGroup5Selected = uniqIdTracks(uniqIGroup==5);
    idGroup6Selected = uniqIdTracks(uniqIGroup==6);
    idGroup7Selected = uniqIdTracks(uniqIGroup==7);
    idGroup8Selected = uniqIdTracks(uniqIGroup==8);
    idGroup9Selected = uniqIdTracks(uniqIGroup==9);
    % figure, plot(asymTracks',MSDall','.')
    % hold on
    % figure,plot(advanceDistNAs(idGroup1)',distToEdgeLastNAs(idGroup1)','ro'), hold on
    % xlabel('asymmetry')
    % ylabel('MSD')
%     display('Click tracks that belong to group 3 (moving along with protruding edge)')
%     [idGroup3Selected] = showAdhesionTracks(pathForColocalization,'all','tracksNA',newTracksNA);
    % figure, plot(distToEdgeLastNAs(idGroup3)',distToEdgeChangeNAs(idGroup3)','go')
    % plot(advanceDistNAs(idGroup3)',distToEdgeLastNAs(idGroup3)','go')
    % figure, histogram(advanceDistNAs(idGroup1),1:30); hold on; histogram(advanceDistNAs(idGroup3),1:30)
    % figure, plot(asymTracks(idGroup3)',MSDall(idGroup3)','go')
    % create linear discriminator
    save([pathForColocalization filesep 'data' filesep 'selectedGroups.mat'],'idGroup1Selected',...
        'idGroup2Selected','idGroup3Selected','idGroup4Selected','idGroup5Selected','idGroup6Selected',...
        'idGroup7Selected','idGroup8Selected','idGroup9Selected');
end

%% feature extraction
%#1
maxIntensityNAs = arrayfun(@(x) nanmax(x.ampTotal),tracksNA); %this should be high for group 2
% startingIntensityNAs = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),tracksNA); %this should be high for group 2
endingIntensityNAs = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),tracksNA); %this should be high for group 2
decayingIntensityNAs = maxIntensityNAs-endingIntensityNAs; % this will differentiate group 1 vs group 2.
%#2
edgeAdvanceDistNAs = arrayfun(@(x) x.edgeAdvanceDist(x.endingFrameExtra),tracksNA); %this should be also low for group 3
%#3
advanceDistNAs = arrayfun(@(x) x.advanceDist(x.endingFrameExtra),tracksNA); %this should be also low for group 3
%#4
lifeTimeNAs = arrayfun(@(x) x.lifeTime,tracksNA); %this should be low for group 6
%#5
meanIntensityNAs = arrayfun(@(x) nanmean(x.amp),tracksNA); %this should be high for group 2
%#6
distToEdgeFirstNAs = arrayfun(@(x) x.distToEdge(x.startingFrameExtra),tracksNA); %this should be low for group 3
%#7
startingIntensityNAs = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),tracksNA); %this should be high for group 2
%#8
distToEdgeChangeNAs = arrayfun(@(x) x.distToEdgeChange,tracksNA); %this should be low for group 3 and group 5
%#9
distToEdgeLastNAs = arrayfun(@(x) x.distToEdge(x.endingFrameExtra),tracksNA); %this should be low for group 3 and group 7
%#10
edgeAdvanceDistLastChangeNAs =  arrayfun(@(x) x.advanceDistChange2min(x.endingFrameExtra),tracksNA); %this should be negative for group 5 and for group 7
%#11
maxEdgeAdvanceDistChangeNAs =  arrayfun(@(x) x.maxEdgeAdvanceDistChange,tracksNA); %this should be negative for group 5 and for group 7
% Some additional features - will be commented out eventually
asymTracks=arrayfun(@(x) asymDetermination([x.xCoord(logical(x.presence))', x.yCoord(logical(x.presence))']),tracksNA);
MSDall=arrayfun(@(x) sum((x.xCoord(logical(x.presence))'-mean(x.xCoord(logical(x.presence)))).^2+...
    (x.yCoord(logical(x.presence))'-mean(x.yCoord(logical(x.presence)))).^2),tracksNA());

MSDrate = MSDall./lifeTimeNAs;
advanceSpeedNAs = advanceDistNAs./lifeTimeNAs; %this should be also low for group 3
edgeAdvanceSpeedNAs = edgeAdvanceDistNAs./lifeTimeNAs; %this should be also low for group 3
relMovWRTEdge = distToEdgeChangeNAs./lifeTimeNAs;
%% Building classifier...
meas = [decayingIntensityNAs(idGroup1Selected) edgeAdvanceSpeedNAs(idGroup1Selected) advanceSpeedNAs(idGroup1Selected) ...
     lifeTimeNAs(idGroup1Selected) meanIntensityNAs(idGroup1Selected) distToEdgeFirstNAs(idGroup1Selected) ...
     startingIntensityNAs(idGroup1Selected) distToEdgeChangeNAs(idGroup1Selected) distToEdgeLastNAs(idGroup1Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup1Selected) maxEdgeAdvanceDistChangeNAs(idGroup1Selected);
     decayingIntensityNAs(idGroup2Selected) edgeAdvanceSpeedNAs(idGroup2Selected) advanceSpeedNAs(idGroup2Selected) ...
     lifeTimeNAs(idGroup2Selected) meanIntensityNAs(idGroup2Selected) distToEdgeFirstNAs(idGroup2Selected) ...
     startingIntensityNAs(idGroup2Selected) distToEdgeChangeNAs(idGroup2Selected) distToEdgeLastNAs(idGroup2Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup2Selected) maxEdgeAdvanceDistChangeNAs(idGroup2Selected);
     decayingIntensityNAs(idGroup3Selected) edgeAdvanceSpeedNAs(idGroup3Selected) advanceSpeedNAs(idGroup3Selected) ...
     lifeTimeNAs(idGroup3Selected) meanIntensityNAs(idGroup3Selected) distToEdgeFirstNAs(idGroup3Selected) ...
     startingIntensityNAs(idGroup3Selected) distToEdgeChangeNAs(idGroup3Selected) distToEdgeLastNAs(idGroup3Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup3Selected) maxEdgeAdvanceDistChangeNAs(idGroup3Selected);
     decayingIntensityNAs(idGroup4Selected) edgeAdvanceSpeedNAs(idGroup4Selected) advanceSpeedNAs(idGroup4Selected) ...
     lifeTimeNAs(idGroup4Selected) meanIntensityNAs(idGroup4Selected) distToEdgeFirstNAs(idGroup4Selected) ...
     startingIntensityNAs(idGroup4Selected) distToEdgeChangeNAs(idGroup4Selected) distToEdgeLastNAs(idGroup4Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup4Selected) maxEdgeAdvanceDistChangeNAs(idGroup4Selected);
     decayingIntensityNAs(idGroup5Selected) edgeAdvanceSpeedNAs(idGroup5Selected) advanceSpeedNAs(idGroup5Selected) ...
     lifeTimeNAs(idGroup5Selected) meanIntensityNAs(idGroup5Selected) distToEdgeFirstNAs(idGroup5Selected) ...
     startingIntensityNAs(idGroup5Selected) distToEdgeChangeNAs(idGroup5Selected) distToEdgeLastNAs(idGroup5Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup5Selected) maxEdgeAdvanceDistChangeNAs(idGroup5Selected);
     decayingIntensityNAs(idGroup6Selected) edgeAdvanceSpeedNAs(idGroup6Selected) advanceSpeedNAs(idGroup6Selected) ...
     lifeTimeNAs(idGroup6Selected) meanIntensityNAs(idGroup6Selected) distToEdgeFirstNAs(idGroup6Selected) ...
     startingIntensityNAs(idGroup6Selected) distToEdgeChangeNAs(idGroup6Selected) distToEdgeLastNAs(idGroup6Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup6Selected) maxEdgeAdvanceDistChangeNAs(idGroup6Selected);
     decayingIntensityNAs(idGroup7Selected) edgeAdvanceSpeedNAs(idGroup7Selected) advanceSpeedNAs(idGroup7Selected) ...
     lifeTimeNAs(idGroup7Selected) meanIntensityNAs(idGroup7Selected) distToEdgeFirstNAs(idGroup7Selected) ...
     startingIntensityNAs(idGroup7Selected) distToEdgeChangeNAs(idGroup7Selected) distToEdgeLastNAs(idGroup7Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup7Selected) maxEdgeAdvanceDistChangeNAs(idGroup7Selected);
     decayingIntensityNAs(idGroup8Selected) edgeAdvanceSpeedNAs(idGroup8Selected) advanceSpeedNAs(idGroup8Selected) ...
     lifeTimeNAs(idGroup8Selected) meanIntensityNAs(idGroup8Selected) distToEdgeFirstNAs(idGroup8Selected) ...
     startingIntensityNAs(idGroup8Selected) distToEdgeChangeNAs(idGroup8Selected) distToEdgeLastNAs(idGroup8Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup8Selected) maxEdgeAdvanceDistChangeNAs(idGroup8Selected);
     decayingIntensityNAs(idGroup9Selected) edgeAdvanceSpeedNAs(idGroup9Selected) advanceSpeedNAs(idGroup9Selected) ...
     lifeTimeNAs(idGroup9Selected) meanIntensityNAs(idGroup9Selected) distToEdgeFirstNAs(idGroup9Selected) ...
     startingIntensityNAs(idGroup9Selected) distToEdgeChangeNAs(idGroup9Selected) distToEdgeLastNAs(idGroup9Selected) ...
     edgeAdvanceDistLastChangeNAs(idGroup9Selected) maxEdgeAdvanceDistChangeNAs(idGroup9Selected)];
% meas = [advanceDistNAs(idGroup4Selected) edgeAdvanceDistNAs(idGroup4Selected);
%     advanceDistNAs(nonGroup24Selected) edgeAdvanceDistNAs(nonGroup24Selected)];
nG1=length(idGroup1Selected);
nG2=length(idGroup2Selected);
nG3=length(idGroup3Selected);
nG4=length(idGroup4Selected);
nG5=length(idGroup5Selected);
nG6=length(idGroup6Selected);
nG7=length(idGroup7Selected);
nG8=length(idGroup8Selected);
nG9=length(idGroup9Selected);
nTotalG = nG1+nG2+nG3+nG4+nG5+nG6+nG7+nG8+nG9;
species = cell(nTotalG,1);
for ii=1:nTotalG
    if ii<=nG1
        species{ii} = 'Group1';
    elseif ii<=nG1+nG2
        species{ii} = 'Group2';
    elseif ii<=nG1+nG2+nG3
        species{ii} = 'Group3';
    elseif ii<=nG1+nG2+nG3+nG4
        species{ii} = 'Group4';
    elseif ii<=nG1+nG2+nG3+nG4+nG5
        species{ii} = 'Group5';
    elseif ii<=nG1+nG2+nG3+nG4+nG5+nG6
        species{ii} = 'Group6';
    elseif ii<=nG1+nG2+nG3+nG4+nG5+nG6+nG7
        species{ii} = 'Group7';
    elseif ii<=nG1+nG2+nG3+nG4+nG5+nG6+nG7+nG8
        species{ii} = 'Group8';
    elseif ii<=nTotalG
        species{ii} = 'Group9';
    end
end
T = table(meas(:,1),meas(:,2),meas(:,3),meas(:,4),meas(:,5),meas(:,6),meas(:,7),meas(:,8),meas(:,9),meas(:,10),meas(:,11),species,...
    'VariableNames',{'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', ...
    'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', ...
    'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistLastChangeNAs',...
    'maxEdgeAdvanceDistChangeNAs','Group'});
%% visualize feature space and p-dist (similarity) matrix
features = meas;
nGroups = 9;
% normalize features
for i = 1 : size(features,2)
    features(:,i) = (features(:,i) - min(features(:,i)))./(max(features(:,i)) - min(features(:,i)));
end
figure; imagesc(features');hold on
c = colorbar;
c.Label.String = 'feature value';

print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'featureSpace.eps']);

% find(strcmp(response,'Group1'))
% find(strcmp(response,'Group2'))
% find(strcmp(response,'Group3'))
% find(strcmp(response,'Group4'))
% find(strcmp(response,'Group5'))
% find(strcmp(response,'Group6'))
% find(strcmp(response,'Group7'))
% find(strcmp(response,'Group8'))
% find(strcmp(response,'Group9'))
D = pdist(features);
D1 =  squareform(D);
figure; imagesc(D1);
title('similarityAmongTrainedData')
c = colorbar;
c.Label.String = 'p-dist';
for ii=1:nGroups
    x0 = find(strcmp(species,['Group' num2str(ii)]),1)-1;
    w = sum(strcmp(species,['Group' num2str(ii)]));
    rectangle('Position',[x0 x0 w w],'EdgeColor','w','LineWidth',0.5)
end

print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'similarityAmongTrainedData.eps']);

Dfeats = pdist(features');
Dfeats1 =  squareform(Dfeats);
figure; imagesc(Dfeats1); title('similarityAmongFeatures')
c = colorbar;
c.Label.String = 'p-dist';

print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'similarityAmongFeatures.eps']);

%% classifier training or import
useDefinedClassifier=input('Do you want to use already defined classifier? [y/(n)]: ','s');
if isempty(useDefinedClassifier)
    useDefinedClassifier='n';
end

if strcmp(useDefinedClassifier,'n')
    [trainedClassifier, validationAccuracy, C, order] = trainClassifier(T);
    disp(['Validation accuracy is ' num2str(validationAccuracy) '.'])
    save([pathForColocalization filesep 'data' filesep 'trainedClassifier.mat'],'trainedClassifier')

    % normalize confusion matrix
    for ii=1:size(C,1)
        C(ii,:) = C(ii,:)/sum(C(ii,:));
    end
    figure; imagesc(C); title('Confusion Matrix')
    c = colorbar;
    c.Label.String = 'normalized prediction';
    print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'confusionMatrix.eps']);
    disp('The order is :')
    disp(order)
elseif strcmp(useDefinedClassifier,'y')
    [fileDefinedClassifier, pathDefinedClassifier]=uigetfile('*.mat', 'Select the mat file for trained classifier');
    trainedClassifier=load(fullfile(pathDefinedClassifier,fileDefinedClassifier));
    trainedClassifier=trainedClassifier.trainedClassifier;
    disp(['Classifier at ' fullfile(pathDefinedClassifier,fileDefinedClassifier) ' is used for cross-validation.'])
    [validationAccuracy, C, order] = validateClassifier(trainedClassifier,T);
    disp(['Cross-validation accuracy is ' num2str(validationAccuracy) '.'])
    
    % normalize confusion matrix
    for ii=1:size(C,1)
        C(ii,:) = C(ii,:)/sum(C(ii,:));
    end
    figure; imagesc(C); title('Confusion Matrix')
    c = colorbar;
    c.Label.String ='normalized prediction';
    print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'confusionMatrix_otherClassifier.eps']);
    disp('The order is :')
    disp(order)
end

allData = [decayingIntensityNAs edgeAdvanceSpeedNAs advanceSpeedNAs ...
     lifeTimeNAs meanIntensityNAs distToEdgeFirstNAs ...
     startingIntensityNAs distToEdgeChangeNAs distToEdgeLastNAs ...
     edgeAdvanceDistLastChangeNAs maxEdgeAdvanceDistChangeNAs];
allDataClass = predict(trainedClassifier,allData);
idGroup1 = strcmp(allDataClass,'Group1');
idGroup2 = strcmp(allDataClass,'Group2');
idGroup3 = strcmp(allDataClass,'Group3');
idGroup4 = strcmp(allDataClass,'Group4');
idGroup5 = strcmp(allDataClass,'Group5');
idGroup6 = strcmp(allDataClass,'Group6');
idGroup7 = strcmp(allDataClass,'Group7');
idGroup8 = strcmp(allDataClass,'Group8');
idGroup9 = strcmp(allDataClass,'Group9');

colors = distinguishable_colors(9,'k');
figure, imshow(imgMap(:,:,end),[])
hold on
htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(1,:)),tracksNA(idGroup1),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(1,:)),tracksNA(idGroup1));
htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(2,:)),tracksNA(idGroup2),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(2,:)),tracksNA(idGroup2));
htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(3,:)),tracksNA(idGroup3),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(3,:)),tracksNA(idGroup3));
htrackG4=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(4,:)),tracksNA(idGroup4),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(4,:)),tracksNA(idGroup4));
htrackG5=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(5,:)),tracksNA(idGroup5),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(5,:)),tracksNA(idGroup5));
htrackG6=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(6,:)),tracksNA(idGroup6),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(6,:)),tracksNA(idGroup6));
htrackG7=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(7,:)),tracksNA(idGroup7),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(7,:)),tracksNA(idGroup7));
htrackG8=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(8,:)),tracksNA(idGroup8),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(8,:)),tracksNA(idGroup8));
htrackG9=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(9,:)),tracksNA(idGroup9),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(9,:)),tracksNA(idGroup9));
legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1} htrackG9{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
    'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
legend('boxoff')
if strcmp(useDefinedClassifier,'n')
    print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.eps']);
    save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','tracksNA')
else
    print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified_otherClassifier.eps']);
    save([pathForColocalization filesep 'data' filesep 'idsClassified_otherClassifier.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','tracksNA')
end
end
function [trainedClassifier, validationAccuracy,C,order] = trainClassifier(datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Train a classifier
template = templateSVM('KernelFunction', 'polynomial', 'PolynomialOrder', 2, 'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1);
trainedClassifier = fitcecoc(predictors, response,'FitPosterior',1, 'Learners', template, 'Coding', 'onevsone', 'PredictorNames', {'decayingIntensityNAs' 'edgeAdvanceSpeedNAs' 'advanceSpeedNAs' 'lifeTimeNAs' 'meanIntensityNAs' 'distToEdgeFirstNAs' 'startingIntensityNAs' 'distToEdgeChangeNAs' 'distToEdgeLastNAs' 'edgeAdvanceDistLastChangeNAs' 'maxEdgeAdvanceDistChangeNAs'}, 'ResponseName', 'Group', 'ClassNames', {'Group1' 'Group2' 'Group3' 'Group4' 'Group5' 'Group6' 'Group7' 'Group8' 'Group9'});

% Perform cross-validation
partitionedModel = crossval(trainedClassifier, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

% confusion matrix
predictedLabels = trainedClassifier.predict(predictors);
[C,order] = confusionmat(response,predictedLabels);

%% Uncomment this section to compute validation predictions and scores:
% % Compute validation predictions and scores
% [validationPredictions, validationScores] = kfoldPredict(partitionedModel);
end
function  [validationAccuracy,C,order] = validateClassifier(trainedClassifier,datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Perform cross-validation
% figure; imagesc(predictors);
% predictedLabels = predict(trainedClassifier,predictors);
% [predictedLabels,NegLoss,PBScore] = trainedClassifier.predict(predictors);
predictedLabels = trainedClassifier.predict(predictors);
results = nan(1,numel(predictedLabels));
for i = 1 : numel(predictedLabels)
    results(i) = strcmp(predictedLabels{i},response{i});
end
validationAccuracy=sum(results)/length(results);
% confusion matrix
[C,order] = confusionmat(response,predictedLabels);
end
%% classificationLearner
% SVMmodel = fitcecoc(meas,species);
% isLoss = resubLoss(SVMmodel)
% CVmodel = crossval(SVMmodel)
% oosLoss = kfoldLoss(CVmodel)
% % %% Here let's filter out group 4 and 5
% % % Look at each class - first is expected to have very negative
% % % advance
% % idGroup4 = [idGroup4Selected];
% % nonGroup135=[idGroup1Selected idGroup3Selected idGroup5Selected];
% % meas = [advanceDistNAs(idGroup4) edgeAdvanceDistNAs(idGroup4);
% %     advanceDistNAs(nonGroup135) edgeAdvanceDistNAs(nonGroup135)];
% % % meas = [advanceDistNAs(idGroup4Selected) edgeAdvanceDistNAs(idGroup4Selected);
% % %     advanceDistNAs(nonGroup24Selected) edgeAdvanceDistNAs(nonGroup24Selected)];
% % nG4=length(idGroup4);
% % nGn4 = length(nonGroup135);
% % nTotalG = nG4+nGn4;
% % species = cell(nTotalG,1);
% % for ii=1:nTotalG
% %     if ii<=nG4
% %         species{ii} = 'Group4';
% %     else 
% %         species{ii} = 'NonGroup4';
% %     end
% % %     elseif ii<=nG4+nGn4
% % %         species{ii} = 'NonGroup4';
% % %     else 
% % %         species{ii} = 'Group3';
% % %     end
% % end
% %% some manual group reassignment
% % species{1}= 'NonGroup4';
% % species{37}= 'Group4';
% % species{36}= 'Group4';
% % species{32}= 'Group4';
% 
% % %% plotting classes
% % mea1 = meas(:,1);
% % mea2 = meas(:,2);
% % figure
% % h1 = gscatter(mea1,mea2,species,'krb','ov^',[],'off');
% % h1(1).LineWidth = 2;
% % h1(2).LineWidth = 2;
% % % h1(3).LineWidth = 2;
% % legend('Group 4','Non-Group 4','Location','best')
% % hold on
% % cls = fitcdiscr(meas,species);
% % % % Plot the classification boundaries.
% % K = cls.Coeffs(1,2).Const; % First retrieve the coefficients for the linear
% % L = cls.Coeffs(1,2).Linear;% boundary between the second and third classes
% %                            % (versicolor and virginica).
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2 = ezplot(f,[min(mea1(1:nG4+nGn4)) max(mea1(1:nG4+nGn4)) min(mea2(1:nG4+nGn4)) max(mea2(1:nG4+nGn4))]); 
% % h2.Color = 'r';
% % h2.LineWidth = 2;
% % xlabel('advanceDistNAs')
% % ylabel('edgeAdvanceDistNAs')
% 
% % showAdhesionTracks(pathForColocalization,70,'tracksNA',tracksNA);
% %% filter out group 4 from the population
% % % allData = [advanceDistNAs edgeAdvanceDistNAs];
% % % allDataClass = predict(cls,allData);
% % % idxGroup4 = strcmp(allDataClass,'Group4');
% % % idxGroup4(idGroup4Selected) = true;
% % % idGroup4 = find(idxGroup4);
% % % %% Working on group 5 - long life time, steadyly low distToEdgeChange expected...
% % % % group 6 - just noise...
% % % % Just get rid of them
% % % % disp('to do here')
% % % idGroup56 = [idGroup5Selected idGroup6Selected];
% % % 
% % % %% Group 2 - maturing
% % % idxMatureNAs = idGroup2Selected;
% % % % maxIntensityNAs(idGroup2Selected)
% % % % maxIntensityNAs(idxMatureNAs)
% % % % figure, histogram(maxIntensityNAs,50)
% % % % % Then, let's filter out group 2 using state information
% % % % idxMatureNAs =  arrayfun(@(x) x.maturing==1,tracksNA);
% % % % % showAdhesionTracks(pathForColocalization,idGroup6Selected,'tracksNA',tracksNA);
% % % % 
% % % % % Here, it is possible that wrong FA segmentation affected wrong group 2
% % % % % definition. We need to rescue tracks whose edge protrude and adhesion
% % % % % also moves along to the edge. In other words, distToEdge relatively
% % % % % constant AND edge protruded.
% % % % distToEdgeLastMatureNAs = arrayfun(@(x) x.distToEdge(x.endingFrame),tracksNA(idxMatureNAs));
% % % % % edgeAdvanceDistMatureNAs = arrayfun(@(x) x.edgeAdvanceDist,tracksNA(idxMatureNAs));
% % % % adhAdvanceDistMatureNAs = arrayfun(@(x) x.advanceDist,tracksNA(idxMatureNAs));
% % % % % figure, histogram(distToEdgeLastMatureNAs,50)
% % % % % figure, histogram(edgeAdvanceDistMatureNAs,50)
% % % % % figure, histogram(adhAdvanceDistMatureNAs,50)
% % % % % figure, plot(distToEdgeChangeMatureNAs,edgeAdvanceDistMatureNAs,'k.')
% % % % % figure, plot(adhAdvanceDistMatureNAs,edgeAdvanceDistMatureNAs,'k.')
% % % % idxRetractingAdhs = adhAdvanceDistMatureNAs<0 & distToEdgeLastMatureNAs>10;
% % % % % matureTracks = tracksNA(idxMatureNAs);
% % % % % reallyMatureTracks = matureTracks(idxRetractingAdhs);
% % % % p=0;
% % % % for ii=find(idxMatureNAs)'
% % % %    p= p+1;
% % % %    idxMatureNAs(ii)=idxRetractingAdhs(p);
% % % % end
% % % 
% % % %% other processes
% % % % meas = [advanceDistNAs(idGroup1Selected);
% % % %     advanceDistNAs(idGroup3Selected)];
% % % % meas = [MSDrate(idGroup1Selected) ;
% % % %     MSDrate(idGroup3Selected) ];
% % % % meas = [relMovWRTEdge(idGroup1Selected) MSDrate(idGroup1Selected) ;
% % % %     relMovWRTEdge(idGroup3Selected) MSDrate(idGroup3Selected) ];
% % % % meas = [MSDrate(idGroup1Selected) relMovWRTEdge(idGroup1Selected);
% % % %     MSDrate(idGroup3Selected) relMovWRTEdge(idGroup3Selected)];
% % % % meas = [MSDrate(idGroup1Selected) advanceDistNAs(idGroup1Selected);
% % % %     MSDrate(idGroup3Selected) advanceDistNAs(idGroup3Selected)];
% % % % meas = [advanceDistNAs(idGroup1Selected), asymTracks(idGroup1Selected);
% % % %     advanceDistNAs(idGroup3Selected),asymTracks(idGroup3Selected)];
% % % % meas = [advanceDistNAs(idGroup1), distToEdgeLastNAs(idGroup1);
% % % %     advanceDistNAs(idGroup3),distToEdgeLastNAs(idGroup3)];
% % % % meas = [asymTracks(idGroup1), MSDall(idGroup1);
% % % %     asymTracks(idGroup3),MSDall(idGroup3)];
% % % 
% % % % nG1=length(idGroup1Selected);
% % % % nG3 = length(idGroup3Selected);
% % % % nTotalG = nG1+nG3;
% % % % species = cell(nTotalG,1);
% % % % for ii=1:nTotalG
% % % %     if ii<=nG1
% % % %         species{ii} = 'group1';
% % % %     else
% % % %         species{ii} = 'group3';
% % % %     end
% % % % end
% % % 
% % % % In case for more than two features:
% % % % linclass = fitcdiscr(meas,species);
% % % % mea1 = meas(:,1);
% % % % mea2 = meas(:,2);
% % % % figure
% % % % h1 = gscatter(mea1,mea2,species,'krb','ov^',[],'off');
% % % % h1(1).LineWidth = 2;
% % % % h1(2).LineWidth = 2;
% % % % legend('Group 1','Group 3','Location','best')
% % % % hold on
% % % % % Plot the classification boundaries.
% % % % K = linclass.Coeffs(1,2).Const; % First retrieve the coefficients for the linear
% % % % L = linclass.Coeffs(1,2).Linear;% boundary between the second and third classes
% % % %                            % (versicolor and virginica).
% % % % % Plot the curve K + [x,y]*L  = 0.
% % % % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % % % h2 = ezplot(f,[0 max(mea1) 0 max(mea2)]);
% % % % h2.Color = 'r';
% % % % h2.LineWidth = 2;
% % % % ylabel('MSD rate')
% % % % xlabel('Adhesion movement relative to edge rate')
% % % 
% % % %% Now separate group1 and 3 using MSD rate
% % % % In case of one feature, I 'll use mean of the interfacing tales of the
% % % % two distributions:
% % % sortedMSDrateG1 = sort(MSDrate(idGroup1Selected));
% % % sortedMSDrateG3 = sort(MSDrate(idGroup3Selected));
% % % thresMSDrate = mean([sortedMSDrateG1(end-1:end); sortedMSDrateG3(1:2)]);
% % % % additional filtering of maturing adhesions and adhesions in the
% % % % retracting edge ....
% % % arrayfun(@(x) plot(x.xCoord,x.yCoord,'b'),tracksNA(idGroup1Selected));
% % % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'bo'),tracksNA(idGroup1Selected));
% % % arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),tracksNA(idGroup3Selected));
% % % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),tracksNA(idGroup3Selected));
% % % showAdhesionTracks(pathForColocalization,idGroup2Selected,'tracksNA',tracksNA);
% % % % predict now
% % % allData = [MSDrate ];
% % % % allData = [relMovWRTEdge MSDrate ];
% % % % allData = [MSDrate relMovWRTEdge];
% % % % allData = [advanceDistNAs MSDall];
% % % % allData = [advanceDistNAs, distToEdgeLastNAs];
% % % % allDataClass = predict(linclass,allData);
% % % % idxStaticNAs = strcmp(allDataClass,'group1');
% % % % idxAdvancingNAs= strcmp(allDataClass,'group3');
% % % idxStaticNAs =allData<thresMSDrate;
% % % idxAdvancingNAs= allData>=thresMSDrate;
% % % idStaticNAs = find(idxStaticNAs);
% % % idStaticNAs = setdiff(idStaticNAs, [idxMatureNAs'; idGroup4; idGroup56']);
% % % idAdvancingNAs = find(idxAdvancingNAs);
% % % idAdvancingNAs = setdiff(idAdvancingNAs, [idxMatureNAs'; idGroup4; idGroup56']);
% % % % For advancing NAs, there is sub-group where adhesions are in stalling
% % % % edge. These should be separated. These should have low edgeAdvanceDist,
% % % % distToEdge, and low distToEdgeChange...
% % % % Look at MSDrate first
% % % figure, histogram(MSDrate(idGroup3Selected),0:1:15)
% % % % For those who have low MSDrate, check edgeAdvanceDist compared to those with high
% % % % MSDrate
% % % MSDrateG3 = MSDrate(idGroup3Selected);
% % % edgeAdvanceG3 = edgeAdvanceDistNAs(idGroup3Selected);
% % % distToEdgeChangeG3 =distToEdgeChangeNAs(idGroup3Selected);
% % % opts = statset('maxIter', 200);
% % % for n = 1:3
% % %     objMSDG3{n} = gmdistribution.fit(MSDrateG3, n, 'Options', opts);
% % % end
% % % [~,idx] = min(cellfun(@(i) i.BIC, objMSDG3));
% % % objMSDG3 = objMSDG3{idx};
% % % [mu,idx] = sort(objMSDG3.mu);
% % % svec = sqrt(squeeze(objMSDG3.Sigma(:,:,idx)));
% % % thresMSDrate= mu(idx(1))+2*svec(idx(1));
% % % idxLowMSDrateG3 = MSDrateG3<thresMSDrate;
% % % MSDrateG3(idxLowMSDrateG3)
% % % MSDrateG3(~idxLowMSDrateG3)
% % % edgeAdvanceG3(idxLowMSDrateG3)
% % % edgeAdvanceG3(~idxLowMSDrateG3)
% % % figure, scatter(MSDrate,edgeAdvanceDistNAs)
% % % hold on
% % % scatter(MSDrate(idGroup4),edgeAdvanceDistNAs(idGroup4),'ro')
% % % scatter(MSDrate(idGroup5Selected),edgeAdvanceDistNAs(idGroup5Selected),'co')
% % % scatter(MSDrate(idGroup6Selected),edgeAdvanceDistNAs(idGroup6Selected),'go')
% % % scatter(MSDrate(idGroup3Selected),edgeAdvanceDistNAs(idGroup3Selected),'ko')
% % % scatter(MSDrate(idGroup3Selected(idxLowMSDrateG3)),edgeAdvanceDistNAs(idGroup3Selected(idxLowMSDrateG3)),'go')
% % % scatter(MSDrate(idGroup1Selected),edgeAdvanceDistNAs(idGroup1Selected),'mo')
% % % 
% % % idStallingNAs;
% % confirm
% % figure, histogram(allData(idxStaticNAs),0:0.1:15); hold on; histogram(allData(idxAdvancingNAs),0:0.1:15)
% % figure, histogram(MSDrate(idGroup1Selected),0:0.1:15); hold on, histogram(MSDrate(idGroup3Selected),0:0.1:15)
% % figure
% % gscatter(allData(:,1),allData(:,2),allDataClass,'krb','ov^',[],'off');
% % legend('Group 1','Group 3','Location','best')
% % hold on
% % h1 = gscatter(mea1,mea2,species,'krb','ov^',[],'off');
% % h1(1).LineWidth = 2;
% % h1(2).LineWidth = 2;
% % % Plot the classification boundaries.
% % h2 = ezplot(f,[min(allData(:,1)) max(allData(:,1)) min(allData(:,2)) max(allData(:,2))]);
% % h2.Color = 'r';
% % h2.LineWidth = 2;
% % xlabel('MSD rate')
% % ylabel('Adhesion movement relative to edge rate')
% 
% figure, imshow(imgMap(:,:,end),[])
% hold on
% htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'g'),tracksNA(idStaticNAs),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'go'),tracksNA(idStaticNAs));
% htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),tracksNA(idAdvancingNAs),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),tracksNA(idAdvancingNAs));
% htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'c'),tracksNA(idxMatureNAs),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'co'),tracksNA(idxMatureNAs));
% htrackG4=arrayfun(@(x) plot(x.xCoord,x.yCoord,'y'),tracksNA(idxGroup4),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'yo'),tracksNA(idxGroup4));
% htrackG5=arrayfun(@(x) plot(x.xCoord,x.yCoord,'m'),tracksNA(idGroup5Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'mo'),tracksNA(idGroup5Selected));
% htrackG6=arrayfun(@(x) plot(x.xCoord,x.yCoord,'w'),tracksNA(idGroup6Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'wo'),tracksNA(idGroup6Selected));
% htrackG7=arrayfun(@(x) plot(x.xCoord,x.yCoord,'w'),tracksNA(idGroup7Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'wo'),tracksNA(idGroup7Selected));
% htrackG8=arrayfun(@(x) plot(x.xCoord,x.yCoord,'w'),tracksNA(idGroup8Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'wo'),tracksNA(idGroup8Selected));
% % arrayfun(@(x) plot(x.xCoord,x.yCoord,'y'),reallyMatureTracks);
% % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'yo'),reallyMatureTracks);
% % idGroup1 = idNAs(idxStaticNAs);
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1}],{'G1:turn-over','G2:maturing','G3:moving along edge',...
%     'G4:retracting','G5:stable','G6:noise or very transient','G7:adhesions at stalling edge','G8:uncartegorized'},'TextColor','w','Location','best')
% legend('boxoff')
% 
% idGroup1 = idxValidNAs(idStaticNAs);
% idGroup2 = idxValidNAs(idxMatureNAs);
% % idGroup3 = idNAs(idxAdvancingNAs);
% idGroup3 = idxValidNAs(idAdvancingNAs);
% idGroup4 = idxValidNAs(idxGroup4);
% idGroup5 = idxValidNAs(idGroup5Selected);
% idGroup6 = idxValidNAs(idGroup6Selected);
% idGroup7 = idxValidNAs(idStallingNAs);
% disp('classified.')
% save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','tracksNA')
% % title('{\bf Linear Classification with Fisher Training Data}')
% 
% % return ids of 
% 
% %% Filtering is necessary
% % Filtering FAs out: tracks whose amplitude goes down during the first
% % minute or start 
% % initialAmpSlope = arrayfun(@(x) x.earlyAmpSlope,tracksNA);
% % thresInitSlope = multithresh(initialAmpSlope,4);
% % minSlope = thresInitSlope(find(thresInitSlope<0,1,'last'));
% % idxSlope = initialAmpSlope>minSlope;
% % %backup old tracks
% % 
% % initialAmpTotal=arrayfun(@(x) x.ampTotal(x.startingFrame),tracksNA);
% % thresAmp = multithresh(initialAmpTotal,2);
% % idxInitAmp = initialAmpTotal<thresAmp(1);
% % 
% % meanDist=arrayfun(@(x) x.advanceDist,tracksNA);
% % thresDist = multithresh(meanDist,10);
% % idxDist = meanDist>thresDist(1);
% % 
% % isNA = idxSlope & idxInitAmp & idxDist;
% % for ii=find(isNA)'
% %     tracksNA(ii).idxNA = true;
% % end
% % for ii=find(~isNA)'
% %     tracksNA(ii).idxNA = false;
% % end
% % newTracksNA = tracksNA(isNA);
% % % figure,hist(initialAmpTotal,200)
% % 
% % % % starting with median, find a edge disconnected with two consequtive
% % % % zeros.
% % % medC = median(pstruct.c);
% % % idxAfterMedC=find(edges>medC);
% % % qq=idxAfterMedC(1);
% % % while N(qq)>0 || N(qq+1)>0
% % %     qq=qq+1;
% % %     if qq>=length(edges)-1
% % %         break
% % %     end
% % % end
% % % idx = pstruct.c<edges(qq);
% % 
% % toc
% %% Make it to feature matrix and save
% % disp('Making it to feature matrix and saving...')
% % numFeatures = 5; % I have five features
% % naFeatMat = zeros(numFeatures,numel(newTracksNA)); 
% % for k=1:numel(newTracksNA)
% %     naFeatMat(:,k) = [newTracksNA(k).earlyAmpSlope; 
% %                                 newTracksNA(k).lateAmpSlope;
% %                                 newTracksNA(k).distToEdgeSlope;
% %                                 newTracksNA(k).edgeAdvanceDist;
% %                                 newTracksNA(k).advanceDist];
% % end
% % %throwing away observation containing NaN
% % idxNan=any(isnan(naFeatMat),1);
% % naFeatMat(:,idxNan)=[];
% % origNAFeatMat = naFeatMat;
% %% scaling of the X to -1 and 1 for slopes and -2 to 2 for distance
% % scaleMat = zeros(size(naFeatMat,1),1);
% % for jj=1:size(naFeatMat,1)
% %     if jj==1 
% %         % Find 2 percentile and 98 percentile
% %         minRow = quantile(naFeatMat(jj,:),.02);
% %         maxRow = quantile(naFeatMat(jj,:),.98);
% %         % Find which one is larger in magnitude
% %         biggerMaxSlop=max(abs(minRow), abs(maxRow));
% %         % Scale with biggerMaxSlop
% %         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
% %         relMinRow = minRow/biggerMaxSlop;
% %         relMaxRow = maxRow/biggerMaxSlop;
% %         %Make data below 2 percentile relMinRow
% %         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
% %         %Make data above 98 percentile one
% %         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
% %         scaleMat(jj)=biggerMaxSlop;
% %     elseif jj==2
% %         % use the same scale as column 1
% %         minRow = quantile(naFeatMat(jj,:),.02);
% %         maxRow = quantile(naFeatMat(jj,:),.98);
% %         % Scale with biggerMaxSlop
% %         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
% %         relMinRow = minRow/biggerMaxSlop;
% %         relMaxRow = maxRow/biggerMaxSlop;
% %         %Make data below 2 percentile relMinRow
% %         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
% %         %Make data above 98 percentile one
% %         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
% %         scaleMat(jj)=biggerMaxSlop;
% %     else
% %         % Find 2 percentile and 98 percentile
% %         minRow = quantile(naFeatMat(jj,:),.02);
% %         maxRow = quantile(naFeatMat(jj,:),.98);
% %         % Find which one is larger in magnitude
% %         biggerMaxSlop = 0.5*max(abs(minRow), abs(maxRow));
% %         % Scale with biggerMaxSlop
% %         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
% %         relMinRow = minRow/biggerMaxSlop;
% %         relMaxRow = maxRow/biggerMaxSlop;
% %         %Make data below 2 percentile relMinRow
% %         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
% %         %Make data above 98 percentile one
% %         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
% %         scaleMat(jj)=biggerMaxSlop;
% %     end
% % end
% % % save
% % save([pathForColocalization filesep 'data' filesep 'naFeatMat.mat'],'origNAFeatMat','naFeatMat','newTracksNA','isNA','tracksNA','scaleMat')
% disp('Done!')
% % Let's put a hard threshold first to see this criteria work. But later I
% % would like to do some type of machine learning 
% % group 1: distance increases AND ampTotal fall down with a negative slope
% 

% case study:
% % group 1: distance increases AND ampTotal fall down to noise level
% ii=166; fRange=tracksNA(ii).startingFrame:tracksNA(ii).endingFrame; figure, plot(fRange,tracksNA(ii).distToEdge(fRange))
% fRange=tracksNA(ii).startingFrameExtra:tracksNA(ii).endingFrameExtra; figure, plot(fRange,tracksNA(ii).ampTotal(fRange),'r');
% % group 2: distance increases AND ampTotal keep increasing
% ii=449; fRange=tracksNA(ii).startingFrame:tracksNA(ii).endingFrame; figure, plot(fRange,tracksNA(ii).distToEdge(fRange))
% fRange=tracksNA(ii).startingFrameExtra:tracksNA(ii).endingFrameExtra; figure, plot(fRange,tracksNA(ii).ampTotal(fRange),'r');
% % group 3: distance does not increase AND ampTotal does not increase


% see histogram
% figure,hist(arrayfun(@(x) x.ampSlope,tracksNA),200)

