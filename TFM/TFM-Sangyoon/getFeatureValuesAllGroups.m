function [matrixFeatureValues] = getFeatureValuesAllGroups(tracksNA,feature,idGroups,startingOrEnding)
cellFeatureValues = cell(9,1);
for ii=1:9    
    if strcmp(startingOrEnding,'starting')
        cellFeatureValues{ii} =arrayfun(@(x) getfield(x,{1},feature,{(x.startingFrameExtra)}),tracksNA(idGroups{ii}));
    elseif strcmp(startingOrEnding,'ending')
        cellFeatureValues{ii} =arrayfun(@(x) getfield(x,{1},feature,{(x.endingFrameExtra)}),tracksNA(idGroups{ii}));
%         startingAmpTotal{ii} =arrayfun(@(x) (x.ampTotal(x.endingFrameExtra)),tracksNA(idGroups{ii}));
    else
        cellFeatureValues{ii} =arrayfun(@(x) mean(getfield(x,{1},feature,{':'})),tracksNA(idGroups{ii}));
%         startingAmpTotal{ii} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroups{ii}));
    end
end
[lengthLongestStartingAmpTotal]=max(cellfun(@(x) length(x),cellFeatureValues));

matrixFeatureValues = NaN(lengthLongestStartingAmpTotal,9);
for ii=1:9
    matrixFeatureValues(1:length(cellFeatureValues{ii}),ii) = cellFeatureValues{ii};
end
end