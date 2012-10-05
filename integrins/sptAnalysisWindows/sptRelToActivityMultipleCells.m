function [windowDistFromEdgeComb,sptPropInWindowComb0,...
    sptPropInWindowCombModDiff,sptPropInWindowCombModRatio,...
    sptPropInWindowCombModDiffRatio,sptPropActivityOnset] = ...
    sptRelToActivityMultipleCells(sptPropInWindowInd,windowDistFromEdgeInd)

%% Input

%get number of activity types and number of cells
[numType,numCell] = size(sptPropInWindowInd);

%get global analysis fields
globalField = fieldnames(sptPropInWindowInd);
numGlobalField = length(globalField);

%get fields indicating time relative to activity onset
%and position relative to cell edge
timePosField = fieldnames(windowDistFromEdgeInd);
numTimePosField = length(timePosField);

%% Window distance from edge

%reserve memory for output variable
for iTimePosField = 1 : numTimePosField
    windowDistFromEdgeComb.(timePosField{iTimePosField}) = [];
end
windowDistFromEdgeComb = repmat(windowDistFromEdgeComb,numType,1);

%go over activity types
for iType = 1 : numType
    
    %get information for current activity type
    windowDistFromEdgeType = windowDistFromEdgeInd(iType,:);
    
    %go over time and position subfields
    for iTimePosField = 1 : numTimePosField
        
        if iscell(windowDistFromEdgeType(1).(timePosField{iTimePosField}))
            
            windowDistFromEdgeComb(iType).(timePosField{iTimePosField}) = ...
                windowDistFromEdgeType(1).(timePosField{iTimePosField});
            
        else
            
            %get the size of the fields for each cell
            sizeSubSub = NaN(numCell,4);
            for iCell = 1 : numCell
                sizeSubSub(iCell,:) = size(windowDistFromEdgeType(iCell).(timePosField{iTimePosField}));
            end
            
            %get the largest size in each dimension and make that the
            %output size
            maxSize = max(sizeSubSub);
            [indCellMean0,indCellStd] = deal(NaN(maxSize(1),maxSize(2),2,numCell));
            
            %get the individual cell values
            for iCell = 1 : numCell
                tmp = windowDistFromEdgeType(iCell).(timePosField{iTimePosField});
                indCellMean0(1:sizeSubSub(iCell,1),1:sizeSubSub(iCell,2),:,iCell) = tmp(:,:,:,1);
                indCellStd(1:sizeSubSub(iCell,1),1:sizeSubSub(iCell,2),:,iCell) = tmp(:,:,:,2);
            end
            
            %average over all cells
            combCellMean = nanmean(indCellMean0,4);
            combCellStd  = nanstd(indCellMean0,[],4);
            
            %save results in output structure
            tmp = cat(4,combCellMean,combCellStd);
            windowDistFromEdgeComb(iType).(timePosField{iTimePosField}) = tmp;
            
        end
        
    end
    
end

%% Single particle measurements

%reserve memory for output variables
for iGlobalField = 1 : numGlobalField
    sptPropInWindowComb0.(globalField{iGlobalField}) = [];
end
[sptPropInWindowComb0,sptPropInWindowCombModDiff,sptPropInWindowCombModRatio,...
sptPropInWindowCombModDiffRatio,sptPropActivityOnset] = deal(repmat(sptPropInWindowComb0,numType,1));

%go over activity types
for iType = 1 : numType
    
    %get information for current activity type
    sptPropInWindowType = sptPropInWindowInd(iType,:);
    
    %go over global analysis fields
    for iGlobalField = 1 : numGlobalField
        
        %fetch current global analysis field
        propertiesGlobalField = vertcat(sptPropInWindowType.(globalField{iGlobalField}));
        
        if isempty(propertiesGlobalField)
            
            sptPropInWindowComb0(iType).(globalField{iGlobalField}) = [];
            sptPropInWindowCombModDiff(iType).(globalField{iGlobalField}) = [];
            sptPropInWindowCombModRatio(iType).(globalField{iGlobalField}) = [];
            sptPropInWindowCombModDiffRatio(iType).(globalField{iGlobalField}) = [];
            
        else
            
            %get local analysis sub-fields within global analysis field
            localField = fieldnames(propertiesGlobalField);
            numLocalField = length(localField);
            
            %go over local analysis sub-fields
            for iLocalField = 1 : numLocalField
                
                %get current property
                propertyLocalField = vertcat(propertiesGlobalField.(localField{iLocalField}));
                
                %go over the time and position subfields
                for iTimePosField = 1 : numTimePosField
                    
                    %get single particle measurements
                    propertyTimePos = vertcat(propertyLocalField.(timePosField{iTimePosField}));
                    
                    %get the size of the sub-sub-fields for each cell
                    sizeSubSub = ones(numCell,3);
                    size2nd = length(size(propertyTimePos(1).mean));
                    for iCell = 1 : numCell
                        sizeSubSub(iCell,1:size2nd) = size(propertyTimePos(iCell).mean);
                    end
                    
                    %get the largest size in each dimension and make that the
                    %output size
                    maxSize = max(sizeSubSub);
                    [indCellMean0,indCellStd,indCellNP] = ...
                        deal(NaN(maxSize(1),maxSize(2),maxSize(3),numCell));
                    
                    %get the individual cell mean, std and number of points
                    for iCell = 1 : numCell
                        indCellMean0(1:sizeSubSub(iCell,1),1:sizeSubSub(iCell,2),1:sizeSubSub(iCell,3),iCell) = propertyTimePos(iCell).mean;
                        indCellStd(1:sizeSubSub(iCell,1),1:sizeSubSub(iCell,2),1:sizeSubSub(iCell,3),iCell)  = propertyTimePos(iCell).std;
                        indCellNP(1:sizeSubSub(iCell,1),1:sizeSubSub(iCell,2),1:sizeSubSub(iCell,3),iCell)   = propertyTimePos(iCell).numPoints;
                    end
                    
                    %discard measurements with number of points < 20
                    minNP = 20;
                    indCellMean0(indCellNP<minNP) = NaN;
                    indCellStd(indCellNP<minNP) = NaN;
                    indCellNP(indCellNP<minNP) = NaN;
                    
                    %replace std = 0 with std = average(std)
                    tmp = indCellStd;
                    tmp(tmp==0) = NaN;
                    stdMean = nanmean(tmp,4);
                    for iCell = 1 : numCell
                        tmp = indCellStd(:,:,:,iCell);
                        indx0 = find(tmp==0);
                        tmp(indx0) = stdMean(indx0);
                        indCellStd(:,:,:,iCell) = tmp;
                    end
                    
                    %calculate cell sem
                    indCellSem = indCellStd ./ sqrt(indCellNP);
                    
                    %give each individual cell measurement a weight based on
                    %its sem
                    indCellWeight = 1 ./ (indCellSem.^2);
                    
                    %normalize the weights so that their sum = 1
                    sumIndCellWeight = nansum(indCellWeight,4);
                    indCellWeight = indCellWeight ./ repmat(sumIndCellWeight,[1 1 1 numCell]); %#ok<NASGU>
                    
                    %get property value at protrusion onset just behind the
                    %cell edge as the reference value for all measurements
                    if strcmp(timePosField{iTimePosField},'onset')
                        refValPerCell = indCellMean0(1,1,:,:);
                        tmp4output = squeeze(refValPerCell);
                        if size(tmp4output,2) > 1
                            tmp4output = tmp4output';
                        end
                        sptPropActivityOnset(iType).(globalField{iGlobalField}).(localField{iLocalField}) = tmp4output;
                    end
                    
                    %transform the mean in various ways
                    refValMat = repmat(refValPerCell,[maxSize(1) maxSize(2) 1 1]);
                    indCellMeanModDiff = indCellMean0 - refValMat; %#ok<NASGU> %difference
                    indCellMeanModRatio = indCellMean0 ./ refValMat; %ratio
                    indCellMeanModRatio(isinf(indCellMeanModRatio)) = NaN; %#ok<NASGU>
                    indCellMeanModDiffRatio = (indCellMean0 - refValMat) ./ refValMat; %difference+ratio
                    indCellMeanModDiffRatio(isinf(indCellMeanModDiffRatio)) = NaN; %#ok<NASGU>
                    
                    %                     %get the number of cells contributing to each measurement
                    %                     tmp = ~isnan(indCellNP);
                    %                     combCellNP = nansum(tmp,4);
                    %                     combCellNP(combCellNP==0) = NaN;
                    
                    %go over all modification types and combine
                    %measurements
                    modType = {'0','ModDiff','ModRatio','ModDiffRatio'};
                    for iMod = 1 : length(modType)
                        
                        %multiply each measurement by its weight
                        %this is used to calculate weighted average and to
                        %count number of cells contributing to weighted
                        %average measurement
                        eval(['vecMeanTimesWeight = indCellWeight.*indCellMean' modType{iMod} ';'])
                        
                        %calculate the weighted average of each measurement
                        combCellMean = nansum(vecMeanTimesWeight,4);
                        
                        %get number of cells contributing to the weighted
                        %average
                        combCellNP2 = ~isnan(vecMeanTimesWeight);
                        combCellNP2 = sum(combCellNP2,4);
                        combCellNP2(combCellNP2==0) = NaN;
                        
                        %calculate the weighted std of each measurement
                        eval(['sumSqDiff = nansum( indCellWeight .* (indCellMean' ...
                            modType{iMod} '-repmat(combCellMean,[1 1 1 numCell])).^2,4 );'])
                        denominator = (numCell-1) / numCell;
                        combCellStd = sqrt( sumSqDiff / denominator );
                        
                        %make all measurements with number of cells = 0 as NaN
                        combCellMean(isnan(combCellNP2)) = NaN;
                        combCellStd(isnan(combCellNP2)) = NaN;
                        
                        %save results in output structure
                        tmp2.mean = combCellMean;
                        tmp2.std = combCellStd;
                        tmp2.numPoints = combCellNP2;
                        eval(['sptPropInWindowComb' modType{iMod} ...
                            '(iType).(globalField{iGlobalField}).(localField{iLocalField}).(timePosField{iTimePosField}) = tmp2;'])
                        
                    end
                    
                end %(for iTimePosField = 1 : numTimePosField)
                
            end %(for iLocalField = 1 : numLocalField)
            
        end %(if isempty(propertiesGlobalField))
        
    end %(for iGlobalField = 1 : numGlobalField)
    
end %(for iType = 1 : numType)

%% ~~~ the end ~~~