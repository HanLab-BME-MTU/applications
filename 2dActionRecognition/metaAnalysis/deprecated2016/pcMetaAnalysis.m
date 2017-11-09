%% query inlcudes the followinf fields: 
% baseDirname, 
% measures to process, 
% rule
% 
% Example:
% query.baseDirname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/Debug/';
% query.baseDirname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/';
% query.measures.lbp = true;
% query.measures.localMorphDynam = true;
% query.rule = 'ruleMetasVsNonMetas';%'ruleTumorInd';%
% query.metaDataFname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/PhaseContrast2DExperiments20140430.mat';
function [] = pcMetaAnalysis(query)

tmp = clock;

outputDirname = [query.baseDirname 'output/'];

if ~exist(outputDirname,'dir')
    unix(sprintf('mkdir %s',outputDirname));
end

resultsDirname = [...
    outputDirname, ...
    num2str(tmp(1)), pad(tmp(2),2), pad(tmp(3),2), pad(tmp(4),2), pad(tmp(5),2),...
    getMeasuresString(query.measures),...
    '_' query.rule '/'];

unix(sprintf('mkdir %s',resultsDirname));

[experimentsDirs,labels,labelsStr] = pcGetExperiments(query.baseDirname, query.rule,query.metaDataFname);

uniqueLabelsStr = unique(labelsStr,'stable');
uniqueLabel = unique(labels,'stable');

if isfield(query.measures,'lbp')
    [lbpData,lbpLabels,outLbp] = pcGetLbpData(experimentsDirs,labels,uniqueLabelsStr,uniqueLabel);
    save([resultsDirname 'lbpResults.mat'],'lbpData','lbpLabels','outLbp','labelsStr','resultsDirname');
    visualizeLBP(outLbp,resultsDirname);
end

if isfield(query.measures,'localMorphDynam')
    [sizes,eccentricities,speeds,matchingScores,morphDimLabels, outMorphDynam] = gcGetLocalMorphoDynamicsData(experimentsDirs,labels,uniqueLabelsStr,uniqueLabel);
    save([resultsDirname 'localMorphDynam.mat'],'sizes','eccentricities','speeds','matchingScores','outMorphDynam','labelsStr','resultsDirname');
%     visualizeLocalMorphDynamicsHistogramDiff(outMorphDynam,uniqueLabelsStr,resultsDirname);
end

end

function [ans] = getMeasuresString(measures)
ans = '';

if isfield(measures,'lbp')
    ans = [ans '_lbp'];
end

if isfield(measures,'localMorphDynam')
    ans = [ans '_localMorphDynam'];
end
end


function [] = visualizeLBP(out,resultsDirname)
nLabels = length(out.uniqueLabelsStr);
for curInd = 1 : nLabels
    curMap = out.maps{curInd};
    curN = out.Ns(curInd);
    curLableStr = out.uniqueLabelsStr{curInd};%uniqueLabel(curInd)
    h = figure;
    imagescnan(curMap);
    hold on;
    title(sprintf('%s (N = %d)',curLableStr,curN),'FontSize',28)
    caxis([0,0.03]); colorbar;        
    haxes = get(h,'CurrentAxes');    
    set(haxes,'XTick',[]);
    set(haxes,'YTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTickLabel',[]);
    set(haxes,'FontSize',28);    
    xlabel('1st PC','FontSize',28); ylabel('2nd PC','FontSize',28);
    hold off;
    idx=regexp([' ' curLableStr],'(?<=\s+)\S','start')-1;
    curLableStr(idx)=upper(curLableStr(idx));
    curLableStr(ismember(curLableStr,' ,.:;!')) = [];
    eval(sprintf('print -djpeg %s', [resultsDirname sprintf('lbpDistribution_%s.jpg',curLableStr)]));
    eval(sprintf('print -dbmp16m  %s', [resultsDirname sprintf('lbpDistribution_%s.bmp',curLableStr)]));
end
end
