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
idxCloseToEdge = arrayfun(@(x) mean(x.distToEdge)<100000/pixSize, tracksAD);
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
iNA=2; iFC=3; iFA=4; iBA=1;
[~,adhVel]=arrayfun(@(x) regression(x.iFrame(x.startingFrame:x.endingFrame), x.advanceDist(x.startingFrame:x.endingFrame)),tracksAD);
idxNegativeVel = adhVel<0;
firstBAidxCell = arrayfun(@(x) find(x.state==iBA,1,'first'),tracksAD,'UniformOutput',false);
% firstBAidxCellEmpty = cellfun(@(x) isempty(x),firstBAidxCell);
% firstBAidx = ~firstBAidxCellEmpty .* cellfun(@
firstNAidx = arrayfun(@(x) find(x.state==iNA,1,'first'),tracksAD,'UniformOutput',false);
firstFAidx = arrayfun(@(x) find(x.state==iFC | x.state==iFA,1,'first'),tracksAD,'UniformOutput',false);
% firstNAidx = arrayfun(@(x) find(strcmp(x.state,'NA'),1,'first'),tracksAD,'UniformOutput',false);
% firstFAidx = arrayfun(@(x) find(strcmp(x.state,'FA') | strcmp(x.state,'FC'),1,'first'),tracksAD,'UniformOutput',false);
idxStartingFromNAToFAcell = cellfun(@(x,y) ~isempty(x) & ~isempty(y) & x<y, firstNAidx, firstFAidx,'UniformOutput',false);
idxEvenBeenFA = cellfun(@(x) ~isempty(x), firstFAidx,'UniformOutput',false);
idxCloseToEdge = arrayfun(@(y) y.distToEdge(y.startingFrame) <100000/pixSize, tracksAD,'UniformOutput',false);
idxStartingFromNAToFA = cellfun(@(x) ~isempty(x),idxStartingFromNAToFAcell);
idxEverBeenFAAndCloseToEdge = cellfun(@(x,y) ~isempty(x) & ~isempty(y),idxEvenBeenFA,idxCloseToEdge);
idxCloseToEdge2 = arrayfun(@(x) mean(x.distToEdge)<500000/pixSize, tracksAD);
idxCloseToEdgeStarting = arrayfun(@(x) x.distToEdge(x.startingFrame) <100000/pixSize, tracksAD);
indMaturingNA = (idxStartingFromNAToFA) & idxCloseToEdge2 & idxCloseToEdgeStarting; % &  | idxEverBeenFAAndCloseToEdge  & idxNegativeVel
idxStartingFromNACell = cellfun(@(x,y) ~isempty(x) & ~isempty(y) & x<y, firstBAidxCell, firstNAidx,'UniformOutput',false);
idxStartingFromNA= cellfun(@(x) ~isempty(x),idxStartingFromNACell);
indNonmaturingNA = (idxStartingFromNA) & indNAatEdge & ~indMaturingNA;
%% 3) maturing FAs from existing FAs: G5
% For this, adhesions should have 1) negative overall velocity; 2) should
% be always 'FA' or 'FC'; 3) associated area should increase; 
iFA=4; iFC=3;
AlwaysFAcell = arrayfun(@(x) x.state==iFA | x.state==iFC,tracksAD,'unif',false);
lifeTimeCell = arrayfun(@(x)x.lifeTime+1,tracksAD,'unif',false);
idxAlwaysFAcell = cellfun(@(x,y) sum(x)==y,AlwaysFAcell,lifeTimeCell);
[~,faGrowth]=arrayfun(@(x) regression(x.iFrame(x.startingFrame:x.endingFrame), x.area(x.startingFrame:x.endingFrame)),tracksAD(idxAlwaysFAcell));
indexAlwaysFAcell=find(idxAlwaysFAcell);
positiveFAGrowth=false(size(idxAlwaysFAcell));
positiveFAGrowth(indexAlwaysFAcell(faGrowth>0))=true;
indMaturingFA = positiveFAGrowth; %idxNegativeVel & 
%% 4) decaying FAs from existing FAs: G4
% For this, adhesions should have 1) negative overall velocity; 2) should
% be always 'FA' or 'FC'; 3) associated area should decrease; 
negativeFAGrowth=false(size(idxAlwaysFAcell));
negativeFAGrowth(indexAlwaysFAcell(faGrowth<0))=true;
indDecayingFA = negativeFAGrowth; %idxNegativeVel & 
%% 5) Others
indOthers = ~indDecayingFA & ~indMaturingFA & ~indMaturingNA & ~indNonmaturingNA;
%% Assign dynamicType
for ii=find(indNonmaturingNA)'
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
indAll={indNonmaturingNA,indMaturingNA,allFalse,indDecayingFA,indMaturingFA, indOthers,allFalse,allFalse,allFalse}; %To be compatible with 9-class-based
disp('Done!')
end