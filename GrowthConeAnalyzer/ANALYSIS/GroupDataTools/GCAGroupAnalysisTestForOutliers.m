function [ output_args ] = GCAGroupAnalysisPlottingDistributionInTime(toPlotGroup,varargin)
%GCAGroupAnalysisTestForOutliers: view distributions per cell for different
%parameters and screen 
%% INPUT
% toPlotGroup: assumes data has been collected. 
% interactive: true

%% Check input
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlotGroup')
ip.addParameter('interactive',true,@(x) islogical(x))
ip.addParameter('outputDirectory',pwd,@(x) ischar(x)); 
ip.parse(toPlotGroup,varargin{:});
p = ip.Results;
%%
 
  saveDir =  [ip.Results.outputDirectory filesep 'Test Measurements Distributions']; 
  
  if  ~isdir(saveDir) ; 
      mkdir(saveDir) ; 
  end 

params = fieldnames(toPlotGroup);
params = params(cellfun(@(x) ~strcmpi(x,'info'),params)); 

if ip.Results.interactive == true
    paramSelect  = listSelectGUI(params,[],'move');
    params = params(paramSelect);
end


names = toPlotGroup.info.projList{1}(:,2);
measurementNames = cellfun(@(x) strrep(x,'_',' '),params,'uniformoutput',0);

for iParam = 1:numel(params);
    dataMat = toPlotGroup.(params{iParam}).dataMat{1}; 
    measurementNameC = measurementNames{iParam}; 
    gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,'sampleNames',names,'measurementName',measurementNameC);
    if ~isempty(ip.Results.outputDirectory)
        saveas(gcf,[saveDir filesep params{iParam} '.fig']) ;
    end
    close gcf
end





end

