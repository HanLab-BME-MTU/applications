
%% Colocalization
% The core function used for colocalization analysis is colocalMeasurePt2Cnt
MD = MovieData.load('/home2/avega/Documents/test_0001/test_0001.mat');%Indicate location of MD object file from scriptGeneralColocalization analysis
p = MD.getProcess(3).getParameters();
orgChannel = 2; %Original continuum channel used in scriptGeneralColocalization
p.ChannelObs = 3; %New context channel
p.SearchRadius = 2; %Radius around detection to use to read out intensity 
p.RandomRuns = 1;% Number of times randomized data is analyzed
MD.getProcess(3).setParameters(p);
MD.getProcess(3).run;


%% Comparison
% The core function used for colocalization analysis is multiChannelColocalization
% Load original colocalization results
load([p.OutputDirectory 'ColocalizationPt2Cnt/colocalInfo' num2str(p.ChannelRef) num2str(orgChannel) '.mat'],'enrichInd')
colocEnrich = enrichInd;
% Load new colocalization results
load([p.OutputDirectory 'ColocalizationPt2Cnt/colocalInfo' num2str(p.ChannelRef) num2str(p.ChannelObs) '.mat'],'enrichInd','randEnrichInd')
bgEnrich = enrichInd;
randBgEnrich = randEnrichInd;
% Run multi-channel colocalization analysis
visual = 1; %Shows boxplot visual of continuum enrichment at punctate detections separated by high and low context enrichment. Data taken from high/lowColoc outputs
[occupHigh,lowColoc,highColoc] = multiChannelColocalization(bgEnrich,colocEnrich,randBgEnrich,visual); 
% Save results
mkdir([p.OutputDirectory 'MultiChannelColocalization' ])
save([p.OutputDirectory 'MultiChannelColocalization/mccResults' ],'occupHigh','lowColoc','highColoc')
