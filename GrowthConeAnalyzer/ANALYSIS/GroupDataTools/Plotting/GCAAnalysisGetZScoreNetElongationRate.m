function [ output_args ] = GCAAnalysisGetZScoreNetElongationRate(toPlot,varargin)
% GCAAnalysisGetZSCoreNetElongationRate
% 
%% 
ip = inputParser;

ip.addRequired('toPlot');
ip.addParameter('cutOffPValue',0.05);
ip.addParameter('makePlots',true); 

ip.CaseSensitive = false;

ip.parse(toPlot,varargin{:});
%% 

deltasPerGroup = cell(size(toPlot.info.names,2),1);

for iGroup =1:size(toPlot.info.names,2)
    projListC = toPlot.info.projList{iGroup}(:,1);
    neuriteLengthStruct = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'MEASUREMENT_EXTRACTION' ...
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements' ...
        filesep 'neuriteLengthOutput.mat']),projListC);
    
    
    neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0);
    deltasPerGroup{iGroup,1} = cellfun(@(x) (x(end)-x(1)),neuriteLengths);
    
   
    
end

if ~isempty(ip.Results.cutOffPValue)
    [filter,pValues] = arrayfun(@(x) permTest(deltasPerGroup{x,1},deltasPerGroup{1,1},'alpha',ip.Results.cutOffPValue),1:numel(deltasPerGroup));
    
end

if ip.Results.makePlots
    for iGroup = 1:numel(deltasPerGroup)
        toPlot.NetOutgrowth.dataMat{iGroup} = deltasPerGroup{iGroup,1}';
    end
    meas{1} = 'NetOutgrowth';
    GCAGroupAnalysisGroupPlots(toPlot,'plotType','perCell','Interactive',false,'Measurements',meas); 
    
end 

% calculate the z-score


meanForZ = cellfun(@(x) nanmean(x),deltasPerGroup);
% Calculate the standard deviation in the values
stdForZ = cellfun(@(x) nanstd(x),deltasPerGroup);


% Standardize the Difference Metric by the Standard
% Deviation of the Control Distribution
valuesz = (meanForZ-meanForZ(1))./stdForZ(1);

if ~isempty(ip.Results.cutOffPValue)
 valuesz(~filter) =0 ; 
end 
% reorder based on hierarchical clustering (for now leave hard wired)
%% sortVector old 
%sortVector  = [1,4,3,6,2,5,7,8]; 
%valuesz = valuesz([1,4,3,6,2,5,7,8]); 
%% sortVector New 
sortVector = [1,4,3,6,2,5,8,7]; 
valuesz = valuesz(sortVector); 
% if ~isempty(ip.Results.cutOffPValue)
%     filter = filter(sortVector); 
% end 
    
sanityCheck = toPlot.info.names(1,(sortVector)); 
figure

valuesz = valuesz(2:end); 
valuesz = valuesz'; 
%                     toPlot.(params{iParam}).zScore = real(valuesz);
imagesc(valuesz,[-3,3])
set(gca,'XTickLabel',sanityCheck(2:end)); 
set(gca,'XTickLabelRotation',45); 
set(gca,'YTick',[1])
cMap = brewermap(11,'Rdbu'); 
cMap = flip(cMap,1); 
set(gcf,'ColorMap',cMap);
colorbar
saveas(gcf,'ZScoresNetNeuriteElongationRate.fig'); 
saveas(gcf,'ZScoresNetNeuriteElognationRate.eps','psc2'); 
save('valuesZ.mat','valuesz'); 
save('groupNames.mat','sanityCheck'); 
end

