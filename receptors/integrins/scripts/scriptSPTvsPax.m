
for i = 1 : length(movieStructAlphaVPax);
    
    disp(num2str(i))
    
    if movieStructAlphaVPax(i).activityLevel > 0
        
        tmp = movieStructAlphaVPax(i).fileName{1};
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaV\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        load diffusionModeClassification.mat
        load numNeighborsPerTrack.mat
        
        seqOfEvents = vertcat(tracksFinal.seqOfEvents);
        maxFrames = max(seqOfEvents(:,1));
        
        lengthMinMax = [5 99];
        firstPaxImFile = [topDir '\imagesPax\' ls([topDir '\imagesPax\*0001.tif'])];
        firstPaxMaskFile = [topDir '\analysisFAs\masks\' ls([topDir '\analysisFAs\masks\*0001.tif'])];
        firstCellMaskFile = [topDir '\analysisCellEdgeSmall\SegmentationPackage\refined_masks\refined_masks_for_channel_1\' ...
            ls([topDir '\analysisCellEdgeSmall\SegmentationPackage\refined_masks\refined_masks_for_channel_1\*0001.tif'])];
        saveDir = [topDir '\analysisFAs'];

        trackDiffPaxInfo = sptDiffvsPaxInt(tracksFinal,diffModeAnalysisRes,numNeighbors,1:400:maxFrames+1,firstPaxImFile,firstPaxMaskFile,firstCellMaskFile,saveDir,1,lengthMinMax);
        close all

    end
    
end
