function [] = pcMetaAnalysisMD(query)
% Given a query, performs meta analysis
%% query inlcudes the followinf fields: 
% dataDirname, 
% analysisDirname,
% measures to process, 
% rule
% 
%% Example:
% query.dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/PrimaryMelanoma/';
% query.analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/';
% query.measures.lbp = true;
% query.measures.localMorphDynam = true;
% query.rule = 'ruleJobTumor';%'ruleMetasVsNonMetas';%'ruleTumorInd';%
% query.metaDataFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/PhaseContrast2DExperiments_pilot.mat'
% query.metaDataFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/PhaseContrast2DExperiments20140430_MD.mat'
%
% Assaf Zaritsky, May 2014 
close all;
tmp = clock;

outputDirname = [query.analysisDirname 'output/'];

if ~exist(outputDirname,'dir')
    unix(sprintf('mkdir %s',outputDirname));
end

resultsDirname = [...
    outputDirname, ...
    num2str(tmp(1)), pad(tmp(2),2), pad(tmp(3),2), pad(tmp(4),2), pad(tmp(5),2),...
    getMeasuresString(query.measures),...
    '_' query.rule '/'];

unix(sprintf('mkdir %s',resultsDirname));

[experimentsDirs,labels,labelsStr] = pcGetExperimentsMD(query.dataDirname, query.analysisDirname, query.rule,query.metaDataFname);

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
    visualizeLocalMorphDynamicsHistograms(outMorphDynam,resultsDirname);
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

%% Visualization
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
    % file name and save figure
    idx=regexp([' ' curLableStr],'(?<=\s+)\S','start')-1;
    curLableStr(idx)=upper(curLableStr(idx));
    curLableStr(ismember(curLableStr,' ,.:;!')) = [];
    eval(sprintf('print -djpeg %s', [resultsDirname sprintf('lbpDistribution_%s.jpg',curLableStr)]));
    eval(sprintf('print -dbmp16m  %s', [resultsDirname sprintf('lbpDistribution_%s.bmp',curLableStr)]));
end
end

%% 
function [] = visualizeLocalMorphDynamicsHistograms(out,resultsDirname)
nLabels = length(out.sizes.uniqueLabelsStr);
for curInd = 1 : nLabels
    curN = out.sizes.Ns(curInd);
    curEccentircities = out.eccentricities.hists{curInd};
    curMatchingScores = out.matchingScores.hists{curInd};
    curSizes = out.sizes.hists{curInd};
    curSpeeds = out.speeds.hists{curInd};
    
    nBins = length(curEccentircities);
    
    curLableStr = out.sizes.uniqueLabelsStr{curInd};%uniqueLabel(curInd)
    h = figure;   
    hold on;
    title(sprintf('%s (N = %d)',curLableStr,curN),'FontSize',28);
    plot(1:nBins,curEccentircities,'-pr','MarkerFaceColor','r','MarkerSize',5);
    plot(1:nBins,curMatchingScores,'-og','MarkerFaceColor','g','MarkerSize',5);
    plot(1:nBins,curSizes,'-sc','MarkerFaceColor','c','MarkerSize',5);
    plot(1:nBins,curSpeeds,'-dk','MarkerFaceColor','k','MarkerSize',5);
    ylim([0,0.15])
    haxes = get(h,'CurrentAxes');    
    %     set(haxes,'XTick',[]);
    %     set(haxes,'YTick',[]);
    %     set(haxes,'XTickLabel',[]);
    %     set(haxes,'YTickLabel',[]);
    set(haxes,'FontSize',28);    
    xlabel('Bins','FontSize',28); ylabel('%','FontSize',28);
    %     legend('Eccentricity','Matching score','Size','Speed','FontSize',16,'Location','SouthOutside');
    hold off;
    % file name and save figure
    idx=regexp([' ' curLableStr],'(?<=\s+)\S','start')-1;
    curLableStr(idx)=upper(curLableStr(idx));
    curLableStr(ismember(curLableStr,' ,.:;!')) = [];
    eval(sprintf('print -djpeg %s', [resultsDirname sprintf('localMorphoDynamicDistributions_%s.jpg',curLableStr)]));
    eval(sprintf('print -dbmp16m  %s', [resultsDirname sprintf('localMorphoDynamicDistributions_%s.bmp',curLableStr)]));
end
end
