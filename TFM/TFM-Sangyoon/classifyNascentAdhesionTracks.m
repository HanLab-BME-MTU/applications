function [idGroup1,idGroup2,idGroup3]=classifyNascentAdhesionTracks(pathForColocalization,varargin)
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

% There should be multiple steps for outlier detection and exclusion
% 1. Adhesions that has small adhesion density should be thrown away
tracksNA = calculateAdhesionSpatialDensity(tracksNA,35); % This need to be more smart if time allows..
numNeisAll = arrayfun(@(x) x.numNeighbors,tracksNA);
% figure, histogram(numNeisAll,100)
thresNumNeis = multithresh(numNeisAll,2);
idxNAsHiDensity = numNeisAll>thresNumNeis(1);

% 2. tracks whose amplitudes are decreasing in the early phase (first one
% minute)
initialAmpSlope = arrayfun(@(x) x.earlyAmpSlope,tracksNA);

% set critical slope as local minimum in the negative value closest to zero
[Namp,edgesHistAmp]= histcounts(initialAmpSlope);
idxUnderZeroAmpSlope = edgesHistAmp<0;
idxLocMinAmpSlope = false(size(idxUnderZeroAmpSlope));
LocMinAmpSlope = locmin1d(Namp);
idxLocMinAmpSlope(LocMinAmpSlope)=true;
idxLocMinSlopeAmpClosestToZero = find(idxLocMinAmpSlope & idxUnderZeroAmpSlope, 1, 'last');
critAmpSlope = edgesHistAmp(idxLocMinSlopeAmpClosestToZero)/3;
% critAmpSlope = 0;
idxNAsIncreasingAmp = initialAmpSlope>critAmpSlope;

% 3.tracks whose closest edge is not moving over the period of the track
% make 0 to NaN for closesBdPoint
for ii=1:numel(tracksNA)
    idxZeros = tracksNA(ii).closestBdPoint(:)==0;
    tracksNA(ii).closestBdPoint(idxZeros)=NaN;
end
% edgeMSD = arrayfun(@(x) nanmean((x.closestBdPoint(:,1)'-nanmean(x.closestBdPoint(:,1))).^2+...
%     (x.closestBdPoint(:,2)'-nanmean(x.closestBdPoint(:,2))).^2),tracksNA(idxNAsHiDensity&idxNAsIncreasingAmp));
edgeMSD = arrayfun(@(x) nanmean((x.closestBdPoint(:,1)'-nanmean(x.closestBdPoint(:,1))).^2+...
    (x.closestBdPoint(:,2)'-nanmean(x.closestBdPoint(:,2))).^2),tracksNA);
% figure, histogram(edgeMSD,1000)
thresEdgeMSD = 2;%multithresh(edgeMSD,2);
idxNAsEdgeMSD = edgeMSD>thresEdgeMSD;

% 4. Filter out NAs that formed inside the edge
distToEdgeStart = arrayfun(@(x) x.distToEdge(x.startingFrame),tracksNA);
thresDistToEdge = multithresh(distToEdgeStart,1);
idxNAsFormedInside = arrayfun(@(x) x.distToEdge(x.startingFrame)<thresDistToEdge,tracksNA);
% tracksNA = tracksNA(idxNAsFormedInside);

% figure, imshow(imgMap(:,:,end),[])
% hold on
% arrayfun(@(x) plot(x.xCoord,x.yCoord,'g'),tracksNA);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'go'),tracksNA);
% arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),tracksNA(idxNAsHiDensity));
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),tracksNA(idxNAsHiDensity));
% arrayfun(@(x) plot(x.xCoord,x.yCoord,'c'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity));
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'co'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity));
% arrayfun(@(x) plot(x.xCoord,x.yCoord,'m'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity &idxNAsEdgeMSD));
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'mo'),tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity & idxNAsEdgeMSD));

% Now filter out all the outliers
tracksNA = tracksNA(idxNAsIncreasingAmp & idxNAsHiDensity & idxNAsEdgeMSD & idxNAsFormedInside);

% Then, let's filter out group 2 using state information
idxMatureNAs =  arrayfun(@(x) x.maturing==1,tracksNA);

% Here, it is possible that wrong FA segmentation affected wrong group 2
% definition. We need to rescue tracks whose edge protrude and adhesion
% also moves along to the edge. In other words, distToEdge relatively
% constant AND edge protruded.
distToEdgeLastMatureNAs = arrayfun(@(x) x.distToEdge(x.endingFrame),tracksNA(idxMatureNAs));
% edgeAdvanceDistMatureNAs = arrayfun(@(x) x.edgeAdvanceDist,tracksNA(idxMatureNAs));
adhAdvanceDistMatureNAs = arrayfun(@(x) x.advanceDist,tracksNA(idxMatureNAs));
% figure, histogram(distToEdgeLastMatureNAs,50)
% figure, histogram(edgeAdvanceDistMatureNAs,50)
% figure, histogram(adhAdvanceDistMatureNAs,50)
% figure, plot(distToEdgeChangeMatureNAs,edgeAdvanceDistMatureNAs,'k.')
% figure, plot(adhAdvanceDistMatureNAs,edgeAdvanceDistMatureNAs,'k.')
idxRetractingAdhs = adhAdvanceDistMatureNAs<0 & distToEdgeLastMatureNAs>10;
% matureTracks = tracksNA(idxMatureNAs);
% reallyMatureTracks = matureTracks(idxRetractingAdhs);
p=0;
for ii=find(idxMatureNAs)'
   p= p+1;
   idxMatureNAs(ii)=idxRetractingAdhs(p);
end
% Then, separate group 1 and 3 by dist to Edge at the ending frame
% group 3
% idxAdvancingNAs =  arrayfun(@(x) x.distToEdge(x.endingFrame)<thresDistToEdge,tracksNA);%arrayfun(@(x) x.edgeAdvanceDist<thresEdgeAdv,tracksNA);
% idxMatureNAs = idxMatureNAs & idxNAsFormedInside; % refinement of group 2 by eliminating ones at the edge
% group 1: static
% asymTracks=arrayfun(@(x) asymDetermination([x.xCoord(logical(x.presence))', x.yCoord(logical(x.presence))']),tracksNA(~idxMatureNAs));
% MSDall=arrayfun(@(x) sum((x.xCoord(logical(x.presence))'-mean(x.xCoord(logical(x.presence)))).^2+...
%     (x.yCoord(logical(x.presence))'-mean(x.yCoord(logical(x.presence)))).^2),tracksNA(~idxMatureNAs));
% distToEdgeLastNAs = arrayfun(@(x) x.distToEdge(x.endingFrame),tracksNA(~idxMatureNAs)); %this should be low for group 3
% distToEdgeChangeNAs = arrayfun(@(x) x.distToEdgeChange,tracksNA(~idxMatureNAs)); %this should be also low for group 3
advanceDistNAs = arrayfun(@(x) x.advanceDist,tracksNA(~idxMatureNAs)); %this should be also low for group 3
% figure, histogram(asymTracks,100)
% figure, histogram(advanceDistNAs,200)
% thresMSD = multithresh(MSD,3);

% apply linear discriminant analysis
% data training.
outputFile=[pathForColocalization filesep 'data' filesep 'selectedGroups.mat'];
reuseSelectedGroups = 'n';
% See if you can use existing tracks
if exist(outputFile,'file')
    reuseSelectedGroups=input('selectedGroups.mat is found. Do you want to use it? (y/n): ');
    if strcmp(reuseSelectedGroups, 'y')
        idGroups = load(outputFile,'idGroup1Selected','idGroup3Selected');
        idGroup1Selected = idGroups.idGroup1Selected;
        idGroup3Selected = idGroups.idGroup3Selected;
    end
end
if ~exist(outputFile,'file') || strcmp(reuseSelectedGroups, 'n')
    display('Click tracks that belong to group 1 (forming and disassembling as edge protrude)')
    newTracksNA=tracksNA(~idxMatureNAs);
    idNAs = find(~idxMatureNAs);
    [idGroup1Selected] = showAdhesionTracks(pathForColocalization,'all','tracksNA',newTracksNA);
    % figure, plot(asymTracks',MSDall','.')
    % hold on
    % figure,plot(advanceDistNAs(idGroup1)',distToEdgeLastNAs(idGroup1)','ro'), hold on
    % xlabel('asymmetry')
    % ylabel('MSD')
    display('Click tracks that belong to group 3 (moving along with protruding edge)')
    [idGroup3Selected] = showAdhesionTracks(pathForColocalization,'all','tracksNA',newTracksNA);
    % figure, plot(distToEdgeLastNAs(idGroup3)',distToEdgeChangeNAs(idGroup3)','go')
    % plot(advanceDistNAs(idGroup3)',distToEdgeLastNAs(idGroup3)','go')
    % figure, histogram(advanceDistNAs(idGroup1),1:30); hold on; histogram(advanceDistNAs(idGroup3),1:30)
    % figure, plot(asymTracks(idGroup3)',MSDall(idGroup3)','go')
    % create linear discriminator
    save([pathForColocalization filesep 'data' filesep 'selectedGroups.mat'],'idGroup1Selected','idGroup3Selected')
end

meas = [advanceDistNAs(idGroup1Selected);
    advanceDistNAs(idGroup3Selected)];
% meas = [advanceDistNAs(idGroup1), distToEdgeLastNAs(idGroup1);
%     advanceDistNAs(idGroup3),distToEdgeLastNAs(idGroup3)];
% meas = [asymTracks(idGroup1), MSDall(idGroup1);
%     asymTracks(idGroup3),MSDall(idGroup3)];
nG1=length(idGroup1Selected);
nG3 = length(idGroup3Selected);
nTotalG = nG1+nG3;
species = cell(nTotalG,1);
for ii=1:nTotalG
    if ii<=nG1
        species{ii} = 'group1';
    else
        species{ii} = 'group3';
    end
end

linclass = fitcdiscr(meas,species);
distEdge = meas(:,1);
% relDist = meas(:,2);
% figure
% h1 = gscatter(distEdge,relDist,species,'krb','ov^',[],'off');
% h1(1).LineWidth = 2;
% h1(2).LineWidth = 2;
% legend('Group 1','Group 3','Location','best')
% hold on
% % Plot the classification boundaries.
% K = linclass.Coeffs(1,2).Const; % First retrieve the coefficients for the linear
% L = linclass.Coeffs(1,2).Linear;% boundary between the second and third classes
%                            % (versicolor and virginica).
% % Plot the curve K + [x,y]*L  = 0.
% f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% h2 = ezplot(f,[0 max(distEdge) 0 max(relDist)]);
% h2.Color = 'r';
% h2.LineWidth = 2;
% predict now
allData = [advanceDistNAs];
% allData = [advanceDistNAs, distToEdgeLastNAs];
allDataClass = predict(linclass,allData);
idxStaticNAs = strcmp(allDataClass,'group1');
idxAdvancingNAs= strcmp(allDataClass,'group3');
% confirm
figure, histogram(allData(idxStaticNAs),1:30); hold on; histogram(allData(idxAdvancingNAs),1:30)
figure, imshow(imgMap(:,:,end),[])
hold on
arrayfun(@(x) plot(x.xCoord,x.yCoord,'g'),newTracksNA(idxStaticNAs));
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'go'),newTracksNA(idxStaticNAs));
arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),newTracksNA(idxAdvancingNAs));
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),newTracksNA(idxAdvancingNAs));
arrayfun(@(x) plot(x.xCoord,x.yCoord,'c'),tracksNA(idxMatureNAs));
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'co'),tracksNA(idxMatureNAs));
% arrayfun(@(x) plot(x.xCoord,x.yCoord,'y'),reallyMatureTracks);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'yo'),reallyMatureTracks);
idGroup1 = idNAs(idxStaticNAs);
idGroup2 = find(idxMatureNAs);
idGroup3 = idNAs(idxAdvancingNAs);
disp('classified.')
save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','newTracksNA','tracksNA')
% title('{\bf Linear Classification with Fisher Training Data}')

% return ids of 

%% Filtering is necessary
% Filtering FAs out: tracks whose amplitude goes down during the first
% minute or start 
% initialAmpSlope = arrayfun(@(x) x.earlyAmpSlope,tracksNA);
% thresInitSlope = multithresh(initialAmpSlope,4);
% minSlope = thresInitSlope(find(thresInitSlope<0,1,'last'));
% idxSlope = initialAmpSlope>minSlope;
% %backup old tracks
% 
% initialAmpTotal=arrayfun(@(x) x.ampTotal(x.startingFrame),tracksNA);
% thresAmp = multithresh(initialAmpTotal,2);
% idxInitAmp = initialAmpTotal<thresAmp(1);
% 
% meanDist=arrayfun(@(x) x.advanceDist,tracksNA);
% thresDist = multithresh(meanDist,10);
% idxDist = meanDist>thresDist(1);
% 
% isNA = idxSlope & idxInitAmp & idxDist;
% for ii=find(isNA)'
%     tracksNA(ii).idxNA = true;
% end
% for ii=find(~isNA)'
%     tracksNA(ii).idxNA = false;
% end
% newTracksNA = tracksNA(isNA);
% % figure,hist(initialAmpTotal,200)
% 
% % % starting with median, find a edge disconnected with two consequtive
% % % zeros.
% % medC = median(pstruct.c);
% % idxAfterMedC=find(edges>medC);
% % qq=idxAfterMedC(1);
% % while N(qq)>0 || N(qq+1)>0
% %     qq=qq+1;
% %     if qq>=length(edges)-1
% %         break
% %     end
% % end
% % idx = pstruct.c<edges(qq);
% 
% toc
%% Make it to feature matrix and save
% disp('Making it to feature matrix and saving...')
% numFeatures = 5; % I have five features
% naFeatMat = zeros(numFeatures,numel(newTracksNA)); 
% for k=1:numel(newTracksNA)
%     naFeatMat(:,k) = [newTracksNA(k).earlyAmpSlope; 
%                                 newTracksNA(k).lateAmpSlope;
%                                 newTracksNA(k).distToEdgeSlope;
%                                 newTracksNA(k).edgeAdvanceDist;
%                                 newTracksNA(k).advanceDist];
% end
% %throwing away observation containing NaN
% idxNan=any(isnan(naFeatMat),1);
% naFeatMat(:,idxNan)=[];
% origNAFeatMat = naFeatMat;
%% scaling of the X to -1 and 1 for slopes and -2 to 2 for distance
% scaleMat = zeros(size(naFeatMat,1),1);
% for jj=1:size(naFeatMat,1)
%     if jj==1 
%         % Find 2 percentile and 98 percentile
%         minRow = quantile(naFeatMat(jj,:),.02);
%         maxRow = quantile(naFeatMat(jj,:),.98);
%         % Find which one is larger in magnitude
%         biggerMaxSlop=max(abs(minRow), abs(maxRow));
%         % Scale with biggerMaxSlop
%         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
%         relMinRow = minRow/biggerMaxSlop;
%         relMaxRow = maxRow/biggerMaxSlop;
%         %Make data below 2 percentile relMinRow
%         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
%         %Make data above 98 percentile one
%         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
%         scaleMat(jj)=biggerMaxSlop;
%     elseif jj==2
%         % use the same scale as column 1
%         minRow = quantile(naFeatMat(jj,:),.02);
%         maxRow = quantile(naFeatMat(jj,:),.98);
%         % Scale with biggerMaxSlop
%         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
%         relMinRow = minRow/biggerMaxSlop;
%         relMaxRow = maxRow/biggerMaxSlop;
%         %Make data below 2 percentile relMinRow
%         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
%         %Make data above 98 percentile one
%         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
%         scaleMat(jj)=biggerMaxSlop;
%     else
%         % Find 2 percentile and 98 percentile
%         minRow = quantile(naFeatMat(jj,:),.02);
%         maxRow = quantile(naFeatMat(jj,:),.98);
%         % Find which one is larger in magnitude
%         biggerMaxSlop = 0.5*max(abs(minRow), abs(maxRow));
%         % Scale with biggerMaxSlop
%         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
%         relMinRow = minRow/biggerMaxSlop;
%         relMaxRow = maxRow/biggerMaxSlop;
%         %Make data below 2 percentile relMinRow
%         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
%         %Make data above 98 percentile one
%         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
%         scaleMat(jj)=biggerMaxSlop;
%     end
% end
% % save
% save([pathForColocalization filesep 'data' filesep 'naFeatMat.mat'],'origNAFeatMat','naFeatMat','newTracksNA','isNA','tracksNA','scaleMat')
disp('Done!')
% Let's put a hard threshold first to see this criteria work. But later I
% would like to do some type of machine learning 
% group 1: distance increases AND ampTotal fall down with a negative slope


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

