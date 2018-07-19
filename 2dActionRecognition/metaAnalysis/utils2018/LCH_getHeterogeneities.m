%%
% data: n x m (# features)
% masterLabels, masterInds - partition criteria
% slaveLabels - labels withint each partition (figure)
function [] = LCH_getHeterogeneities(data,masterUniqueLabels,masterInds,slaveLabels)

nFeats = 10;
%
% if isempty(masterUniqueLabels)
%     LCH_hetero(data,slaveLabels); % [outDname filesep masterStr '_' suffix '.jpg']
% else

hetero = [];
source = [];
metEff = [];
for iMasterLabel = 1 : length(masterUniqueLabels)
    curInds = masterInds{iMasterLabel};
    curData = data(curInds,:);
    curMasterLabel = masterUniqueLabels{iMasterLabel};
    curSlaveLabels = slaveLabels(curInds);
    
    outStd = LCH_hetero(curData,curSlaveLabels,nFeats); % outFname = [outDname filesep masterStr '_' curMasterLabel '_' suffix '.jpg'];
    curSource = LCH_getSource(curMasterLabel);
    curMetEff = LCH_getMetEff(curMasterLabel);
    
    % output
    hetero = [hetero;outStd];
    
    tmpSource = cell(1,size(outStd,1));
    [tmpSource{:}] = deal(curSource);
    source = [source,tmpSource];
    tmpMetEff = cell(1,size(outStd,1));
    [tmpMetEff{:}] = deal(curMetEff);
    metEff = [metEff,tmpMetEff];
    
    % output
    fprintf(sprintf('***  %s  ***\n',curMasterLabel));
    fprintf([repmat(' %.2f ', 1, size(outStd,2)) '\n'], outStd');
    
end

indsCellLines = find(strcmp(source,'CellLines'));
indsTumors = find(strcmp(source,'Tumors'));
indsMelanocytes = find(strcmp(source,'Melanocytes'));

indsHigh = find(strcmp(metEff,'High'));
indsLow = find(strcmp(metEff,'Low'));

nFeats = size(hetero,2);

pvalsHighLow = nan(1,nFeats);
foldHighLow = nan(1,nFeats);
for ifeat = 1 : nFeats
    curStdHigh = hetero(indsHigh,ifeat);
    curStdLow = hetero(indsLow,ifeat);
    pvalsHighLow(ifeat) = ranksum(curStdHigh,curStdLow);
    foldHighLow(ifeat) = mean(curStdHigh) ./ mean(curStdLow);
end


if ~isempty(indsCellLines) && ~isempty(indsMelanocytes)
    pvalsTumorCellLine = nan(1,nFeats);
    pvalsTumorMelano = nan(1,nFeats);
    pvalsCellLineMelano = nan(1,nFeats);
    foldTumorCellLine = nan(1,nFeats);
    foldTumorMelano = nan(1,nFeats);
    foldCellLineMelano = nan(1,nFeats);
    for ifeat = 1 : nFeats
        curStdCellLine = hetero(indsCellLines,ifeat);
        curStdTumor = hetero(indsTumors,ifeat);
        curStdMelano = hetero(indsMelanocytes,ifeat);
        pvalsTumorCellLine(ifeat) = ranksum(curStdCellLine,curStdTumor);
        pvalsTumorMelano(ifeat) = ranksum(curStdTumor,curStdMelano);
        pvalsCellLineMelano(ifeat) = ranksum(curStdCellLine,curStdMelano);
        foldTumorCellLine(ifeat) = mean(curStdTumor)/mean(curStdCellLine);
        foldTumorMelano(ifeat) = mean(curStdTumor)/mean(curStdMelano);
        foldCellLineMelano(ifeat) = mean(curStdMelano)/mean(curStdCellLine);
    end
end
end

function [outStd] = LCH_hetero(mappedFeats,labels,nFeats)
inds = find(~cellfun(@isempty,labels));
mappedFeats = mappedFeats(inds,:);
labels = labels(inds);
uniqueLabels = unique(labels);
nUniqueLables = length(uniqueLabels);

outStd = nan(nUniqueLables,nFeats);

for ilabel = 1 : nUniqueLables
    indsLabel = find(strcmp(labels,uniqueLabels{ilabel}));
    curFeats = mappedFeats(indsLabel,:);
    for ifeat = 1 : nFeats
        outStd(ilabel,ifeat) = std(curFeats(:,ifeat));
    end
end

end