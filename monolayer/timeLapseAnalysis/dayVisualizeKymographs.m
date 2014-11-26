function [] = dayVisualizeKymographs(geneDayDiff,mainDirname,propertyStr,metaData)
nGeneDay = length(geneDayDiff);

params.timePerFrame = metaData.timePerFrame;
params.patchSize = 15;

for iGeneDay = 1 : nGeneDay    
    params.pixelSize = geneDayDiff{iGeneDay}.pixelSize;
    dayGeneSeqStr = geneDayDiff{iGeneDay}.dayGeneSeqStr;
    
    geneKymograph = geneDayDiff{iGeneDay}.meanKDKymograph;
    controlKymograph = geneDayDiff{iGeneDay}.meanKControlKymograph;
    diffKymograph = geneDayDiff{iGeneDay}.diffKymograph;
        
    if strcmp(propertyStr,'Speed')
        params.caxis = [0 60];
    else if strcmp(propertyStr,'Directionality')
            params.caxis = [0 8];
        else
            if strcmp(propertyStr,'Coordination')
                params.caxis = [0 1];
            end
        end
    end
    
    params.fname = [mainDirname 'dayGeneControlKymograph/' dayGeneSeqStr '_' propertyStr '_KD.eps'];
    plotKymograph(geneKymograph,params);
    
    params.fname = [mainDirname 'dayGeneControlKymograph/' dayGeneSeqStr '_' propertyStr '_Ctrl.eps'];
    plotKymograph(controlKymograph,params);
    
    params.fname = [mainDirname 'dayGeneControlKymograph/' dayGeneSeqStr '_' propertyStr '_Diff.eps'];
    params.caxis = params.caxis - params.caxis(2)/2;
    plotKymograph(diffKymograph,params);
    
    close all;
end
end