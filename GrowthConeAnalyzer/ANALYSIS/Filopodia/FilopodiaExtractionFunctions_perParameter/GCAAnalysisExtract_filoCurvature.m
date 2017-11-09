function [filoCurvMax] = GCAAnalysisExtract_filoCurvature(analInfo,filoFilterSet)
%% GCAAnalysisExtract_filoCurvature
% Collects Filopodia Max Curvature
%
%
%INPUT:
%
%analInfo :  REQUIRED:  large rx1 structure (currently takes some time to load)
%                       where r is the number of frames
%                       that includes all the information one would want about the
%                       segmentation of the movie (including the filopodia information)
%                       (eventually the analInfo needs to be
%                       saved per frame I think to make faster...too big of file)
%
%filoFilterSet: PARAM:  rx1 cell array where r (row) is the number of
%                       frames and each cell is a r1x1 logical vector
%                       that places a 1 for each filopodia to be included
%                       in the analysis.
%
% 
%
% OUTPUT:
%
% 
% filoCurvMax:   rx1 cell array where r (row) is the number of frames
%                    and each cell holds 
%                    an rx1 double of the max filopodia curvature 
%                    where r is the number of filo in a given frame  
%                    after the filo filter is applied
%%
frames = length(analInfo);

filoCurvMax = cell(frames,1);

for iFrame = 1:length(analInfo)-1
    
    filoInfo = analInfo(iFrame).filoInfo;
    
    if ~isempty(filoInfo)
        filterFrameC= filoFilterSet{iFrame};
        filoInfoFilt  = filoInfo(filterFrameC(:,1));
        
        
        % get the mean and max of the curvature values
        % curvMeans = arrayfun(@(x) mean(abs(filoInfoFilt(x).Ext_FiloCurvIndVals)),1:length(filoInfoFilt));
        curvMax = arrayfun(@(x) max(abs(filoInfoFilt(x).Ext_FiloCurvIndVals)),1:length(filoInfoFilt));
        
        %filoCurvMean{iFrame} = curvMeans';
        filoCurvMax{iFrame} =  curvMax';
    else
        filoCurvMax{iFrame} = [];
    end
end

    
    
end
% if ~isempty(outDir)
% save([outDir filesep 'param_filoLengthToVeil.mat'],'filoLengthToVeil');
% end
