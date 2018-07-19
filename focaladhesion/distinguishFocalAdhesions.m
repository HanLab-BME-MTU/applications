function [tracksAD,indAll] = distinguishFocalAdhesions(tracksAD, MD, outputPath)
%  [tracksAD,indNAatEdge,indMaturingNA,indMaturingFA,indDecayingFA] =
%  distinguishFocalAdhesions(tracksAD, MD, outputPath)
% distinguishes four obviously different adhesion tracks from input
% tracksAD: 1) non-maturing NAs at the edge; 2) maturing NAs to FAs; 3)
% further growing FAs from existing FA state; 4) decaying FAs in terms of
% area, and 5) all the others. These feature will be encoded in the member
% called 'dynamicType'.
% indAll={indNAatEdge,indMaturingNA,[],indDecayingFA,indMaturingFA, indOthers}; %To be compatible with 9-class-based

% Sangyoon Han Feb 2016
if nargin<2
    outputPath = [];
end
% Set up the output file path
if ~isempty(outputPath)
    if ~exist(outputPath,'dir') 
        mkdir(outputPath);
    end
end
% numTracks = numel(tracksAD);
% tInterval = MD.timeInterval_;
pixSize = MD.pixelSize_;
imgWidth=MD.imSize_(2);
imgHeight=MD.imSize_(1);
%% 1) non-maturing NAs at the edge
% For this I will use the criteria that i) they are located within 1 um from the edge
% [~,adhVel]=arrayfun(@(x) regression(x.iFrame(x.startingFrame:x.endingFrame), x.advanceDist(x.startingFrame:x.endingFrame)),tracksAD);
% idxPositiveVel = adhVel>=0;
idxNoBoundary = arrayfun(@(x) isempty(intersect(nanmean(x.closestBdPoint), [0 1 imgWidth imgHeight])), tracksAD);
idxCloseToEdge = arrayfun(@(x) mean(x.distToEdge)<1000/pixSize, tracksAD);
% Check where they are
% FirstImgStack=MD.channels_(1).loadImage(1); 
% figure, imshow(FirstImgStack,[]), hold on
% arrayfun(@(x) plot(x.xCoord(x.startingFrame),x.yCoord(x.startingFrame),'ro'),tracksAD)
% arrayfun(@(x) plot(x.xCoord(x.startingFrame),x.yCoord(x.startingFrame),'co'),tracksAD(idxCloseToEdge & idxNoBoundary))
% iCurrent=2049; plot(tracksAD(iCurrent).xCoord(tracksAD(iCurrent).startingFrame),tracksAD(iCurrent).yCoord(tracksAD(iCurrent).startingFrame),'ro'); idxNoBoundary(iCurrent)
indNAatEdge = idxNoBoundary & idxCloseToEdge;
%% 2) maturing NAs to FAs
% For this, adhesions should have 1) negative overall velocity; 3)
% associated area should increase; 2) starts from NA state or within 1 um;
% 4) should be within 5 um from the edge for entire lifetime

[~,adhVel]=arrayfun(@(x) regression(x.iFrame(x.startingFrame:x.endingFrame), x.advanceDist(x.startingFrame:x.endingFrame)),tracksAD);
idxNegativeVel = adhVel<0;
firstNAidx = arrayfun(@(x) find(strcmp(x.state,'NA'),1,'first'),tracksAD,'UniformOutput',false);
firstFAidx = arrayfun(@(x) find(strcmp(x.state,'FA') | strcmp(x.state,'FC'),1,'first'),tracksAD,'UniformOutput',false);
idxStartingFromNAcell = cellfun(@(x,y) ~isempty(x) & ~isempty(y) & x<y, firstNAidx, firstFAidx,'UniformOutput',false);
idxStartingFromNAcell2 = cellfun(@(x) ~isempty(x), firstFAidx,'UniformOutput',false);
idxStartingFromNAcell3 = arrayfun(@(y) y.distToEdge(y.startingFrame) <1000/pixSize, tracksAD,'UniformOutput',false);
idxStartingFromNA = cellfun(@(x) ~isempty(x),idxStartingFromNAcell);
idxStartingFromNA2 = cellfun(@(x,y) ~isempty(x) & ~isempty(y),idxStartingFromNAcell2,idxStartingFromNAcell3);
idxCloseToEdge2 = arrayfun(@(x) mean(x.distToEdge)<5000/pixSize, tracksAD);
idxCloseToEdgeStarting = arrayfun(@(x) x.distToEdge(x.startingFrame) <1000/pixSize, tracksAD);
indMaturingNA = idxNegativeVel & (idxStartingFromNA | idxStartingFromNA2) & idxCloseToEdge2 & idxCloseToEdgeStarting;
%% 3) maturing FAs from existing FAs
% For this, adhesions should have 1) negative overall velocity; 2) should
% be always 'FA' or 'FC'; 3) associated area should increase; 
AlwaysFAcell = arrayfun(@(x) strcmp(x.state,'FA') | strcmp(x.state,'FC'),tracksAD,'UniformOutput',false);
lifeTimeCell = arrayfun(@(x)x.lifeTime+1,tracksAD,'UniformOutput',false);
idxAlwaysFAcell = cellfun(@(x,y) sum(x)==y,AlwaysFAcell,lifeTimeCell);
[~,faGrowth]=arrayfun(@(x) regression(x.iFrame(x.startingFrame:x.endingFrame), x.area(x.startingFrame:x.endingFrame)),tracksAD(idxAlwaysFAcell));
indexAlwaysFAcell=find(idxAlwaysFAcell);
positiveFAGrowth=false(size(idxAlwaysFAcell));
positiveFAGrowth(indexAlwaysFAcell(faGrowth>0))=true;
indMaturingFA = idxNegativeVel & positiveFAGrowth;
%% 4) decaying FAs from existing FAs
% For this, adhesions should have 1) negative overall velocity; 2) should
% be always 'FA' or 'FC'; 3) associated area should decrease; 
negativeFAGrowth=false(size(idxAlwaysFAcell));
negativeFAGrowth(indexAlwaysFAcell(faGrowth<0))=true;
indDecayingFA = idxNegativeVel & negativeFAGrowth;
%% 5) Others
indOthers = ~indDecayingFA & ~indMaturingFA & ~indMaturingNA & ~indNAatEdge;
%% Assign dynamicType
for ii=find(indNAatEdge)'
    tracksAD(ii).dynamicType=1;
end
for ii=find(indMaturingNA)'
    tracksAD(ii).dynamicType=2;
end
for ii=find(indMaturingFA)'
    tracksAD(ii).dynamicType=3;
end
for ii=find(indDecayingFA)'
    tracksAD(ii).dynamicType=4;
end
for ii=find(indOthers)'
    tracksAD(ii).dynamicType=5;
end
if ~isempty(outputPath)
    disp(['saving in ' outputPath '.'])
    save([outputPath filesep 'dynamicTypes.mat'], 'tracksAD', 'indNAatEdge', 'indMaturingNA','indMaturingFA','indDecayingFA','indOthers')
end
allFalse = false(size(indNAatEdge));
indAll={indNAatEdge,indMaturingNA,allFalse,indDecayingFA,indMaturingFA, indOthers,allFalse,allFalse,allFalse}; %To be compatible with 9-class-based
disp('Done!')
end