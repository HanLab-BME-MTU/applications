load('/work/bioinformatics/shared/dope/export/backup_dopeAnnotator_output_22May2017_0157-anevarez.mat')
doneAnnotations = cellAnnotations;
load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat')
doneKeys = {doneAnnotations.cellData.key};
doneKeys = unique(doneKeys);

cellDataSetAll = cellDataSet;
cellMat = cell2mat(cellDataSetAll);
allCell_Keys = {cellMat.key};
allCellKeys = unique(allCell_Keys);

remainingKeys = setxor(allCellKeys, doneKeys);
subSetData = cell(1, length(remainingKeys));
j = 1;
for i = 1:length(allCellKeys)
    
   if ismember({cellDataSetAll{i}.key}, remainingKeys)
    subSetData{j} = cellDataSetAll{i};
    j = j + 1;
   end
    
end




matFileName = '/work/bioinformatics/shared/dope/data/OMETIFF/dopeCheckPointGen2Gen3_24MAY2017_andres.mat';
cellDataSet = subSetData;
sessionAnnotations_24MAY2017 = doneAnnotations;
save(matFileName, 'cellDataSet', 'annotationSet','sessionAnnotations_24MAY2017','allCellKeys','doneKeys','remainingKeys');