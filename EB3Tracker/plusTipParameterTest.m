function [allData]=plusTipParameterTest(projData,parametersToTest,secPerFrame,pixSizeNm)
% diagnostic test of tracking parameter settings

if nargin<4 
   error('--parametersToTest: not enough input parameters')
end



paramNames=parametersToTest(:,1);
paramDflts=parametersToTest(:,2);
paramRange=parametersToTest(:,3);


% check whether min/max radii values are ok - min can't be greater than max
% in any iteration
minRadDef=paramDflts{cellfun(@(x) isequal(x,'minRadius'),paramNames)};
maxRadDef=paramDflts{cellfun(@(x) isequal(x,'maxRadius'),paramNames)};
minRadAll=paramRange{cellfun(@(x) isequal(x,'minRadius'),paramNames)};
maxRadAll=paramRange{cellfun(@(x) isequal(x,'maxRadius'),paramNames)};

if any(minRadDef>maxRadAll) || any(minRadAll>maxRadDef)
    error('minRadius must be less than maxRadius')
end

% loop thru parameters and run through range of values to test
nParams=length(paramNames);
allData=cell(nParams,1);
c=1;
for i=1:nParams
    [data,dataCols]=paramTest(projData,paramNames,paramDflts,secPerFrame,pixSizeNm,paramNames{i},paramRange{i});
    allData{i,1}=data;
    nVal=length(paramRange{i});
    % keep track of which parameter we've tested
    [allDataList{c:c+nVal-1,1}]=deal(paramNames{i});
    c=c+nVal;
end


paramDir=[formatPath(projData.anDir) filesep 'paramTest'];

% combine data with parameter names and column descriptors
allData=cell2mat(allData);
allData=mat2cell(allData,ones(size(allData,1),1),ones(size(allData,2),1));
allData=[dataCols; [allDataList allData]];

save([paramDir filesep 'allData'],'allData')



function [data,dataCols]=paramTest(projData,paramNames,paramDflts,secPerFrame,pixSizeNm,pName,pRange)

% pName is the name of the parameter to test.  It can be:
% 'timeWindow'
% 'minRadius'
% 'maxRadius'
% 'maxFAngle'
% 'maxBAngle'
% 'maxShrinkFactor'
% 'fluctRad'

% pRange should be a vector containing the range of parameter pName to test

% data column descriptions
dataCols={' '...
    'maxTgapFrames'...
    'maxTgapSec'...
    'minRadius'...
    'maxRadius'...
    'maxFangle'...
    'maxBangle'...
    'maxShrinkFactor'...
    'fluctRad'...
    'nSubTracks'...
    'nFullTracks'...
    '%FullGrowthOnly'...
    '%FullWithFGap'...
    '%FullWithBGap'...
    '%GapsAtMaxTGap'...
    'dispToNNdistRatio'...
    'medFullTrackLftFrames'...
    'medFractionTimeGrowing'...
            };




% set the defaults
minRadius=paramDflts{cellfun(@(x) isequal(x,'minRadius'),paramNames)};
maxRadius=paramDflts{cellfun(@(x) isequal(x,'maxRadius'),paramNames)};
timeWindow=paramDflts{cellfun(@(x) isequal(x,'timeWindow'),paramNames)};
maxFAngle=paramDflts{cellfun(@(x) isequal(x,'maxFAngle'),paramNames)};
maxBAngle=paramDflts{cellfun(@(x) isequal(x,'maxBAngle'),paramNames)};
maxShrinkFactor=paramDflts{cellfun(@(x) isequal(x,'maxShrinkFactor'),paramNames)};
fluctRad=paramDflts{cellfun(@(x) isequal(x,'fluctRad'),paramNames)};

% test for 50 frames
timeRange=[1 50];
diagnostics=1;

nIter=length(pRange);

% define all parameters as vectors of repeated entires
minRadiusAll=minRadius.*ones(1,nIter);
maxRadiusAll=maxRadius.*ones(1,nIter);
timeWindowAll=timeWindow.*ones(1,nIter);
maxFAngleAll=maxFAngle.*ones(1,nIter);
maxBAngleAll=maxBAngle.*ones(1,nIter);
maxShrinkFactorAll=maxShrinkFactor.*ones(1,nIter);
fluctRadAll=fluctRad.*ones(1,nIter);

% reassign the one we're testing with the range in pRange
switch pName
    case 'minRadius'
        minRadiusAll=pRange;
    case 'maxRadius'
        maxRadiusAll=pRange;
    case 'timeWindow'
        timeWindowAll=pRange;
    case 'maxFAngle'
        maxFAngleAll=pRange;
    case 'maxBAngle'
        maxBAngleAll=pRange;
    case 'maxShrinkFactor'
        maxShrinkFactorAll=pRange;
    case 'fluctRad'
        fluctRadAll=pRange;
end



figDir=[formatPath(projData.anDir) filesep 'paramTest' filesep 'figs_' pName];
if isdir(figDir)
    rmdir(figDir,'s')
end
mkdir(figDir);

close all

s1=2;
strg1 = sprintf('%%.%dd',s1);
data=nan(nIter,14);
for i=1:nIter
    
    % assign specific parameter pRange
    minRadius=minRadiusAll(i);
    maxRadius=maxRadiusAll(i);
    timeWindow=timeWindowAll(i);
    maxFAngle=maxFAngleAll(i);
    maxBAngle=maxBAngleAll(i);
    maxShrinkFactor=maxShrinkFactorAll(i);
    fluctRad=fluctRadAll(i);

    % run tracking and post-processing
    plusTipTracker(projData,timeWindow,3,minRadius,maxRadius,maxFAngle,maxBAngle,maxShrinkFactor,fluctRad,timeRange,diagnostics)
    [projData]=plusTipPostTracking(projData,secPerFrame,pixSizeNm);
    
    % assign strings for figure names
    timeWindowStr = sprintf(strg1,timeWindow);
    minRadiusStr = sprintf(strg1,minRadius);
    maxRadiusStr = sprintf(strg1,maxRadius);
    maxFAngleStr = sprintf(strg1,maxFAngle);
    maxBAngleStr = sprintf(strg1,maxBAngle);
    maxShrinkFactorStr=sprintf('%3.1f',maxShrinkFactor);
    fluctRadStr=sprintf('%3.1f',fluctRad);
    
    % plots 1 and 2 are from forward/backward iterations of linking before
    % final forward, which is the one we save as plot 3.
    saveas(3,[figDir filesep 'linkingDistance_' minRadiusStr '_' maxRadiusStr '.tif']);
    saveas(4,[figDir filesep 'gapLifetimes_' timeWindowStr '_' maxFAngleStr '_' maxBAngleStr '_' maxShrinkFactorStr '_' fluctRadStr '.tif']);
    close all

    temp=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

    % add parameters to data array
    data(i, 1)=timeWindow;
    data(i, 2)=timeWindow*secPerFrame;
    data(i, 3)=minRadius;
    data(i, 4)=maxRadius;
    data(i, 5)=maxFAngle;
    data(i, 6)=maxBAngle;
    data(i, 7)=maxShrinkFactor;
    data(i, 8)=fluctRad;

    % number of subtracks
    data(i, 9)=size(temp,1);
    % number of full tracks
    data(i,10)=projData.numTracks;
    % percent of full tracks with growth only
    data(i,11)=100*(length(setdiff([1:projData.numTracks]',[projData.tracksWithFgap; projData.tracksWithFgap]))/projData.numTracks);
    % percent of full tracks with a forward gap
    data(i,12)=100*(length(projData.tracksWithFgap)/projData.numTracks);
    % percent of full tracks with a backward gap
    data(i,13)=100*(length(projData.tracksWithBgap)/projData.numTracks);
    % percent of all gaps that are tGapMax long
    idx=find(temp(:,5)==2 | temp(:,5)==3);
    gapLft=temp(idx,6); m=max(gapLft);
    data(i,14)=100*(sum(gapLft==m)/length(gapLft));
    
    % mean frame-frame displacement over within-frame NN distance
    data(i,15)=projData.meanDisp2medianNNDistRatio;
    
      
    % get total track lifetime and time spent growing
    tlft=zeros(projData.numTracks,1);
    glft=zeros(projData.numTracks,1);
    for iTrack=1:projData.numTracks
       idx=find(temp(:,1)==iTrack);
       type=temp(idx,5);
       
       tlft(iTrack)=sum(temp(idx,6)); % total time track exists
       glft(iTrack)=sum(temp(idx(type==1),6)); % total time growing
        
    end
    % median total tarck lifetime
    data(i,16)=median(tlft);
    % median fraction of time spent in growth
    data(i,17)=median(glft./tlft);
    
    
end
save([figDir filesep 'data'],'data');

