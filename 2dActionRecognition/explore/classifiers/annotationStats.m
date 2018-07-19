function [] = annotationStats(annotationFname,uniqueLabelsFname)

close all; clc;

if nargin < 1
    %         annotationFname = 'C:\Assaf\LCH_CellExplorer\backup_dopeAnnotator_output_22May2017_0157-anevarez.mat';
    annotationFname = 'C:\Assaf\LCH_CellExplorer\dopeCheckPoint_Anno09Jun2017_1037.mat';
    uniqueLabelsFname = 'C:\Assaf\LCH_CellExplorer\uniqueLabels.mat';
end

load(annotationFname);

cellData = cellAnnotations.cellData;
keys = {cellData.key};
annos = {cellData.annotations};
repeats = cell2mat({cellData.repeatFlag});

indsRepeat = find(repeats);

n = length(annos);

flatAnnos = [];
for i = 1 : n 
    flatAnnos = [flatAnnos annos{i}];
end

uniqueLabels = unique(flatAnnos);

save(uniqueLabelsFname,'uniqueLabels');

nUniqueLabels = length(uniqueLabels);
valsInit = cell(1,nUniqueLabels);
for i = 1 : nUniqueLabels
    valsInit{i} = i;
end
label2indMap = containers.Map(uniqueLabels, valsInit);
annotationConfusionMatrix = zeros(nUniqueLabels,nUniqueLabels);

agreement = 0;
nrepeats = 0;
allRepeatKeys = {};

for i = 1 : length(indsRepeat)    
    repeatIndsBin = strfind(keys, keys{indsRepeat(i)});
    repeatInds = find(not(cellfun('isempty', repeatIndsBin)));

    validateRepeatInds(repeatInds,keys,repeats);
    
    curInd = repeatInds(1);
    curKey = keys{curInd};
    curAnnotation = annos{curInd};
    allRepeatKeys{length(allRepeatKeys) + 1} = curKey;
    for j = 2 : length(repeatInds)
        nrepeats = nrepeats + 1;
        repeatAnnotation = annos{repeatInds(j)};
        if sameAnnotation(curAnnotation,repeatAnnotation)
            agreement = agreement + 1;
        else
            printAnnotations(curKey,curAnnotation,repeatAnnotation);
        end
        
        annotationConfusionMatrix = updateConfusionMatrix(annotationConfusionMatrix,label2indMap,curAnnotation,repeatAnnotation);
        
    end
    
%     if length(repeatInds) ~= 2 
%         fprintf(sprintf('%s: ',keys{repeatInds(1)}));
%         for ir = 1 : length(repeatInds)
%             fprintf(sprintf('%d, ',repeatInds(ir)));
%         end
%         fprintf('\n\n');
%         nRepeatOverTwo = nRepeatOverTwo + 1;
%     end
end

figure; imagesc(annotationConfusionMatrix); colorbar;
tmp = repmat(sum(annotationConfusionMatrix,2),1,12);
annotationConfusionMatrixNorm = annotationConfusionMatrix./tmp;
figure; imagesc(annotationConfusionMatrixNorm); colorbar;%caxis([0,10]);

nUniqueKeyRepeat = unique(allRepeatKeys);



fprintf(sprintf('# cells = %d\n',n));
fprintf(sprintf('Average annotations per cell = %.1f\n',length(flatAnnos)/n));

for il = 1 : length(uniqueLabels)
    curLabel = uniqueLabels{il};
    countLabel = sum(strcmp(curLabel,flatAnnos));
    fprintf(sprintf('%s: %d cells\n',curLabel,countLabel));
end

%% Cell types
% TumorStr = {'m481','m214','m530','m610','um12','m514','m405','m528','m498','ut8','m634','m597'};
% MelanoStr = {'atcc','m116'};
% CellLineStr = {'wm3670','wm1361','wm1366','skmel2','a375','mv3'};
% 
% HighMetStr = {'m481','m214','um12','m514','m405','m634'};
% LowMetStr = {'m530','m610','m528','m498'};


TumorStr = {'m481','m610','m498','m634'};
MelanoStr = {'m116'};
CellLineStr = {'wm3670','wm1361','wm1366','skmel2','a375','mv3'};

HighMetStr = {'m481','m634'};
LowMetStr = {'m610','m498'};

cellTypes = cell(1,n);
nAll = 0;
nTumor = 0;
nMelano = 0;
nCellLine = 0;
nHighMet = 0;
nLowMet = 0;

phenotypeCount = zeros(6,10); % 6 groups (all, tumor, melano, cellLine,highMet,lowMet), 10 phenotypes

annotationCount = zeros(1,nUniqueLabels);

for i = 1 : n
    tmp = strsplit(keys{i},'_');
    cellTypes{i} = lower(tmp{2});
    
    curAnnotations = annos{i};
    curAnnotationInds = findAnnotationInds(curAnnotations,uniqueLabels);
    
    for ia = 1 : length(curAnnotations)
        annotationCount(label2indMap(curAnnotations{ia})) = annotationCount(label2indMap(curAnnotations{ia})) + 1;
    end
    
    if isempty(curAnnotationInds)
        continue;
    end
    
    % All
    nAll = nAll + 1;
    phenotypeCount = addOne(phenotypeCount,1,curAnnotationInds);
    
    % Tumor (+ high/low Met)
    if sum(strcmp(cellTypes{i},TumorStr)) > 0
        nTumor = nTumor + 1;
        phenotypeCount = addOne(phenotypeCount,2,curAnnotationInds);
        
        if sum(strcmp(cellTypes{i},HighMetStr)) > 0
            nHighMet = nHighMet + 1;
            phenotypeCount = addOne(phenotypeCount,5,curAnnotationInds);
        end
        if sum(strcmp(cellTypes{i},LowMetStr)) > 0
            nLowMet = nLowMet + 1;
            phenotypeCount = addOne(phenotypeCount,6,curAnnotationInds);
        end        
    end
    
    if sum(strcmp(cellTypes{i},MelanoStr)) > 0
        nMelano = nMelano + 1;
        phenotypeCount = addOne(phenotypeCount,3,curAnnotationInds);
    end
    
    if sum(strcmp(cellTypes{i},CellLineStr)) > 0
        nCellLine = nCellLine + 1;
        phenotypeCount = addOne(phenotypeCount,4,curAnnotationInds);
    end
end

fprintf(sprintf('\n Cell counts:\nall = %d, tumor = %d (high = %d, low = %d), melanocytes = %d, cell lines = %d\n',...
    nAll,nTumor,nHighMet,nLowMet,nMelano,nCellLine));


distributionCount = zeros(6,10);
distributionCount(1,:) = phenotypeCount(1,:) ./ nAll;
distributionCount(2,:) = phenotypeCount(2,:) ./ nTumor;
distributionCount(3,:) = phenotypeCount(3,:) ./ nMelano;
distributionCount(4,:) = phenotypeCount(4,:) ./ nCellLine;
distributionCount(5,:) = phenotypeCount(5,:) ./ nHighMet;
distributionCount(6,:) = phenotypeCount(6,:) ./ nLowMet;

figure; imagesc(phenotypeCount);
figure; imagesc(distributionCount);

end

%% uniqueLabels 1,2 - junk!
function curAnnotationInds = findAnnotationInds(curAnnotations,uniqueLabels)
curAnnotationInds = [];
for i = 1 : length(curAnnotations)
    curAnn = curAnnotations{i};
    curInd = find(strcmp(curAnn,uniqueLabels));
    curInd = curInd - 2;
    if curInd > 0
        curAnnotationInds = [curAnnotationInds curInd];
    end
end
end
%% 
function phenotypeCount = addOne(phenotypeCount,row,curAnnotationInds)
for i = 1 : length(curAnnotationInds)
    curInd = curAnnotationInds(i);
    phenotypeCount(row,curInd) = phenotypeCount(row,curInd) + 1;
end
end

%%
function bool = sameAnnotation(curAnnotation,repeatAnnotation)

bool = true;

if length(curAnnotation) ~= length(repeatAnnotation)
    bool = false;
    return
end

for i = 1 : length(curAnnotation)    
    IndexC = strfind(repeatAnnotation, curAnnotation{i});
    Index = find(not(cellfun('isempty', IndexC)));
    if length(Index) ~= 1 
        bool = false;
        return
    end
end
end

%% 
function [] = printAnnotations(curKey,curAnnotation,repeatAnnotation)
fprintf(sprintf('%s: \n',curKey));
if isempty(curAnnotation)
    fprintf('[]\n');
else
    for i = 1 : length(curAnnotation)
        fprintf(sprintf('%s, ',curAnnotation{i}));
    end
    fprintf('\n');
end

if isempty(repeatAnnotation)
    fprintf('[]\n');
else    
    for j = 1 : length(repeatAnnotation)
        fprintf(sprintf('%s, ',repeatAnnotation{j}));
    end
    fprintf('\n');
end
fprintf('\n');
end

%% 
function validateRepeatInds(repeatInds,keys,repeats)

if sum(repeats(repeatInds)) ~= (length(repeatInds)-1)
    fprintf(sprintf('%s: %d repeats, %d same key\n',keys{repeatInds(1)},sum(repeats(repeatInds)),length(repeatInds)));    
end

key = keys{repeatInds(1)};
for i = 2 : length(repeatInds)
    assert(strcmp(keys{repeatInds(i)},key));
end
end

%%
function annotationConfusionMatrix = updateConfusionMatrix(annotationConfusionMatrix,label2indMap,curAnnotation,repeatAnnotation)
repeat2cur = true(1,length(repeatAnnotation)); % to know which ones are not
cur2repeat = true(1,length(curAnnotation)); % to know which ones are not
for i = 1 : length(curAnnotation)
    curLabel = curAnnotation{i};
    ind = label2indMap(curLabel);
    IndexC = strfind(repeatAnnotation, curLabel);
    Index = find(not(cellfun('isempty', IndexC)));
    if ~isempty(Index)
        repeat2cur(Index) = false;      
        cur2repeat(i) = false;
        annotationConfusionMatrix(ind,ind) = annotationConfusionMatrix(ind,ind) + 1;            
    end
end

for i = 1 : length(curAnnotation)
    ind1 = label2indMap(curAnnotation{i});
    for j = 1 : length(repeatAnnotation)
        ind2 = label2indMap(repeatAnnotation{j});
        if repeat2cur(j) && cur2repeat(i) 
            annotationConfusionMatrix(min(ind1,ind2),max(ind1,ind2)) = annotationConfusionMatrix(min(ind1,ind2),max(ind1,ind2)) + 1;
        end
    end
end

end