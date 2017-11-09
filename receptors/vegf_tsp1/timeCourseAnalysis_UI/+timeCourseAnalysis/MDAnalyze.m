%% Time Course Analysis (MD-level)
function [MDExtra] = MDAnalyze(MD) %#ok<INUSD>
    %Need blank if not used
    MDExtra.blank = [];
    %{
    %loads 'partitionResult'
    if analysisPara.doPartition
        load(MD.processes_{channel, MD.getProcessIndex('PartitionAnalysisProcess')}.outFilePaths_{1});
    end
    %loads 'tracks' and 'diffAnalysisRes'
    if analysisPara.doPartition
        load(MD.processes_{channel, MD.getProcessIndex('MotionAnalysisProcess')}.outFilePaths_{1});
    end
    %}
end
