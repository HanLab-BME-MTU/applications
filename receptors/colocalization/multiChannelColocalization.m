function [occupHigh,lowColoc,highColoc] = multiChannelColocalization(bgEnrich,colocEnrich,randBgEnrich,visual) 
% MULTICHANNELCOLOCALIZATION measures colocalization characteristics of two channels based on their position relative to a third channel 
%
% Synopsis: [occupHigh, cellCtLow, cellCtHigh,highPt,lowPt] = scriptActinRes(bgEnrich,colocEnrich,randBgEnrich)
% Function finds an intensity threshold in a third channel to separate
% positions of high and low intensity. Function then divides colocalized
% pairs in the first two channels based on their location in the high and
% low intensity regimes. Point of function is to ask whether individual or
% paired properties of colocalized pairs change dependant on their
% relationship to a third channel
%
% Input:
% bgEnrich: enrichInd output from colocalization method colocalMeasurePt2Cnt.
% Channels used in previous method should be punctate channel and
% third/context channel. Cell structure containing enrichment percentage of
% objects in the context channel to detections in the punctate channel
%
% colocEnrich: ratioInd output from colocalization method
% colocalMeasurePt2Cnt. Channels used in previous method should be punctate
% channel and continuum channel (the initial colocalized pair of interest).
% Cell structure containing same information as bgEnrich but for punctate
% and continuum channel.
%
% randBgEnrich: randEnrichInd output from colocalMeasurePt2Cnt. Channels used
% in previous method should be punctate channel and background/context
% channel. Cell structure containing enrichment percentage for randomized
% positions for punctate and context channels.
%
% visual: Enter value of 1 if you would like boxplot visual of high versus
% low intensity enrichment values. Data taken from high/lowColoc outputs
%
% Output:
% occupHigh: 1x2 array for fraction of punctate detections found in high
% actin [actual random]
%
%
% lowColoc:punctate and continuum channel enrichment in low intensity regime of 
% third channel
%
% highColoc: punctate and continuum channel enrichment in high intensity regime of 
% third channel
%
% Anthony Vega 09/2014
%
%% Workspace
% Initialize output
    pTest = zeros(length(colocEnrich),1);
    rTest = zeros(length(colocEnrich),1);

    highColoc = cell(length(colocEnrich),1);
    lowColoc = cell(length(colocEnrich),1);
    plotHighCell = cell(length(colocEnrich),1);
    plotLowCell = cell(length(colocEnrich),1);

for k =1:length(colocEnrich)
    
    % Get enrichment from all randomized positions in background channel and
    % use these as a proxy to mean intensity ratio in background channel
    bgRand = randBgEnrich{k,1};
    bgThreshold = mean(bgRand(:,1));
    
    % Using this mean intensity, see which detections from punctate
    % channel fall above and below this threshold
    bgActual = bgEnrich{k,1};
    indexH = find(bgActual(:,1)>bgThreshold);
    indexL = find(bgActual(:,1)<bgThreshold);

    
    % Save number of high and low
    cellPtHigh = length(indexH);
    cellPtLow = length(indexL);
    
    %Use mean intensity to separate colocalization pairs
    colocPair = colocEnrich{k,1};
    highColoc{k,:} = colocPair(indexH,:);
    lowColoc{k,:} = colocPair(indexL,:);
    
    %For plotting purposes
    plotHigh = colocPair;
    plotHigh(indexL,:) = NaN;
    plotHighCell{k,:} = plotHigh;
    plotLow = colocPair;
    plotLow(indexH,:) = NaN;
    plotLowCell{k,:} = plotLow;

    % Random Test-how many random detections fall above or below threshold
    randHigh = length(find(bgRand(:,1)>bgThreshold));
    randLow = length(find(bgRand(:,1)<bgThreshold));
    
    %Compare fractions in high intensity for real and random data  
    pTest(k,1) = cellPtHigh/(cellPtHigh+cellPtLow) ;
    rTest(k,1) = randHigh/(randHigh+randLow); 

       
end
    occupHigh = [pTest rTest];

%% Visual
if visual ==1
    labelsCnt = {'Continuum Enrichment Low' ,'Continuum Enrichment High'};
    labelsPt = {'Punctate Enrichment Low','Punctate Enrichment High'};
    boxDataHigh = cell2mat(plotHighCell);
    boxDataLow = cell2mat(plotLowCell);
    h(1) = figure; 
    boxplot([boxDataLow(:,1),boxDataHigh(:,1)],'notch','on','labels',labelsCnt)
    ax = h(1).Children;
    ax.YLabel.String = 'Percent';
    h(2) = figure;
    boxplot([boxDataLow(:,2),boxDataHigh(:,2)],'notch','on','labels',labelsPt)
    ax = h(2).Children;
    ax.YLabel.String = 'Percent';
end  

end