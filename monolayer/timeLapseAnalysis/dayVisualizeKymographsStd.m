function [] = dayVisualizeKymographsStd(geneDayDiff,mainDirname,propertyStr,metaData)
nGeneDay = length(geneDayDiff);

params.timePerFrame = metaData.timePerFrame;
params.patchSize = 15;

matDname = [mainDirname 'dayGeneControlKymographStd/mat/'];
if ~exist(matDname,'dir')
    unix(sprintf('mkdir %s',matDname));
end

for iGeneDay = 1 : nGeneDay        
    params.pixelSize = geneDayDiff{iGeneDay}.pixelSize;
    dayGeneSeqStr = geneDayDiff{iGeneDay}.dayGeneSeqStr;
    
    %     if exist([mainDirname 'dayGeneControlKymographStd/' dayGeneSeqStr '_' propertyStr '_Diff.eps'],'file')
    if exist([matDname dayGeneSeqStr '_' propertyStr '_std_Ctrl.mat'],'file')
        continue;
    end
    
    geneStdKymograph = geneDayDiff{iGeneDay}.meanStdKDKymograph;
    controlStdKymograph = geneDayDiff{iGeneDay}.meanStdKControlKymograph;
    %     diffKymograph = geneDayDiff{iGeneDay}.diffKymograph;
        
    if strcmp(propertyStr,'Speed')
        params.caxis = [0 20];
    else if strcmp(propertyStr,'Directionality')
            params.caxis = [0 1];
        else
            if strcmp(propertyStr,'Coordination')
                error('not supports coordination');
            end
        end
    end
    
    params.fname = [mainDirname 'dayGeneControlKymographStd/' dayGeneSeqStr '_' propertyStr '_std_KD.eps'];
    plotKymograph(geneStdKymograph,params);
    save([matDname dayGeneSeqStr '_' propertyStr '_std_KD.mat']);
    
    params.fname = [mainDirname 'dayGeneControlKymographStd/' dayGeneSeqStr '_' propertyStr '_std_Ctrl.eps'];
    plotKymograph(controlStdKymograph,params);
    save([matDname dayGeneSeqStr '_' propertyStr '_std_Ctrl.mat']);
    
    %     params.fname = [mainDirname 'dayGeneControlKymograph/' dayGeneSeqStr '_' propertyStr '_Diff.eps'];
    %     params.caxis = params.caxis - params.caxis(2)/2;
    %     plotKymograph(diffKymograph,params);
    
    close all;
end
end