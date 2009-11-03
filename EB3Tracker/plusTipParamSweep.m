function [allData]=plusTipParamSweep(projData,parametersToTest,secPerFrame,pixSizeNm)
% plusTipParamSweep tests a range of tracking parameter settings to aid selection
%
% Input
%       for all inputs   : Call plusTipParamSweepGUI to set inputs using a GUI
%
% Output
%       allData          : cell array containing results of all tests
%       paramTest folder : subfolder of project folder, which contains a
%                          subfolder for each parameter tested.  these
%                          contain gapLifetime and linkingDistance figures
%                          and individual data (subset of allData) matrices
%
% Note: It is not great coding practice to put called functions
% within the GUI setup, but in this case it is ok since this function will
% never be called elsewhere.  They were previously separate but the similar
% names were confusing.


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


% where the parameter testing data will be stored
paramDir=[formatPath(projData.anDir) filesep 'paramTest'];
if isdir(paramDir)
    rmdir(paramDir,'s')
end

% loop thru parameters and run through range of values to test
nParams=length(paramNames);
c=1; cP=1;
for i=1:nParams
    if isempty(paramRange{i})
        continue
    end
    [data,dataCols]=paramTest(projData,paramNames,paramDflts,secPerFrame,pixSizeNm,paramNames{i},paramRange{i});
    allData{cP,1}=data;
    nVal=length(paramRange{i});
    % keep track of which parameter we've tested
    [allDataList{c:c+nVal-1,1}]=deal(paramNames{i});
    c=c+nVal;
    cP=cP+1;
end

% combine data with parameter names and column descriptors
allData=[dataCols; [allDataList vertcat(allData{:})]];
save([paramDir filesep 'allData'],'allData')



function [data,dataCols]=paramTest(projData,paramNames,paramDflts,secPerFrame,pixSizeNm,pName,pRange)

data=[];
dataCols=[];

% pName is the name of the parameter to test.  It can be:
% 'timeWindow'
% 'minRadius'
% 'maxRadius'
% 'maxFAngle'
% 'maxBAngle'
% 'maxShrinkFactor'
% 'fluctRad'

% pRange should be a vector containing the range of parameter pName to test


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
    plusTipCometTracker(projData,timeWindow,3,minRadius,maxRadius,maxFAngle,maxBAngle,maxShrinkFactor,fluctRad,timeRange,diagnostics)
    [projData]=plusTipPostTracking(projData,secPerFrame,pixSizeNm,[],0);

    % assign strings for figure names
    timeWindowStr = sprintf(strg1,timeWindow);
    minRadiusStr = sprintf(strg1,minRadius);
    maxRadiusStr = sprintf(strg1,maxRadius);
    maxFAngleStr = sprintf(strg1,maxFAngle);
    maxBAngleStr = sprintf(strg1,maxBAngle);
    maxShrinkFactorStr = sprintf('%3.1f',maxShrinkFactor);
    fluctRadStr=sprintf('%3.1f',fluctRad);

    % plots 1 and 2 are from forward/backward iterations of linking before
    % final forward, which is the one we save as plot 3.
    saveas(3,[figDir filesep 'linkingDistance_' minRadiusStr '_' maxRadiusStr '.tif']);
    saveas(4,[figDir filesep 'gapLifetimes_' timeWindowStr '_' maxFAngleStr '_' maxBAngleStr '_' maxShrinkFactorStr '_' fluctRadStr '.tif']);
    close all

    % pick which data to extract
    statNames1={'numTracks';'pair2pairDiffMicPerMinStd';'meanDisp2medianNNDistRatio';'percentFgapsReclass';'percentBgapsReclass'};
    statNames2=fieldnames(projData.stats);

    % GET PROJECT DATA
    count=1; iMov=1;
    % add data pulled from projData
    for iName=1:length(statNames1)
        plusTipData{iMov,count}=projData.(statNames1{iName});
        statNamesData{count}=statNames1{iName};
        count=count+1;
    end

    % add data pulled from projData.stats
    for iName=1:length(statNames2)
        values=projData.stats.(statNames2{iName});
        tempName=statNames2{iName};
        % some measurements have more than one value - here we put each in
        % a separate column and label with 2,3,...
        for v=1:length(values)
            plusTipData{iMov,count}=values(v);
            if v==1
                statNamesData{count}=tempName;
            else
                statNamesData{count}=[tempName '_' num2str(v)];
            end
            count=count+1;
        end
    end

    if i==1
        data=cell(nIter,8+length(plusTipData));
    end
    % add parameters to data array
    data{i,1}=minRadius;
    data{i,2}=maxRadius;
    data{i,3}=timeWindow;
    data{i,4}=timeWindow*secPerFrame;
    data{i,5}=maxFAngle;
    data{i,6}=maxBAngle;
    data{i,7}=maxShrinkFactor;
    data{i,8}=fluctRad;
    data(i,9:end)=plusTipData;

    if i==1
        % data column descriptions
        dataCols=[{projData.anDir...
            'minRadius'...
            'maxRadius'...
            'maxTgapFrames'...
            'maxTgapSec'...
            'maxFangle'...
            'maxBangle'...
            'maxShrinkFactor'...
            'fluctRad'}...
            statNamesData];
    end


end
[allDataList{1:size(data,1),1}]=deal(pName);

dataIndivTest=[dataCols; [allDataList data]];
save([figDir filesep 'dataIndivTest'],'dataIndivTest');