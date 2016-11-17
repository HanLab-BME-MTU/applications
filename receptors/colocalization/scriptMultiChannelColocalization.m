function [] = scriptMultiChannelColocalization(MD)
    %% Colocalization
    % The core function used for colocalization analysis is
    % colocalMeasurePt2Cnt
        p = MD.getProcess(3).getParameters();
        orgChannel = 3; %Original continuum channel used in scriptGeneralColocalization
        p.ChannelObs = 1; %Context channel that you wish to study original colocalization in reference to
        p.SearchRadius = 2; %Radius around detection used for analysis
        p.RandomRuns = 1;% Number of times randomized data is analyzed
        MD.getProcess(3).setParameters(p);
        MD.getProcess(3).run;
    %% Comparison
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
end