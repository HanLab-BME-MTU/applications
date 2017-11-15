% simple script to pre-process the cells for Deep learning.
clc;
clear;
%% Gen Data
load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');

resizeOn = true;
cellSegmentation = true;
newDim = [217 217];

dataRootDir = '/work/bioinformatics/shared/dope/torch/test/217x217/segmented/';
dataRootDirVal = '/work/bioinformatics/shared/dope/torch/test/217x217/segmented/valSet';
dataRootDirTR = '/work/bioinformatics/shared/dope/torch/test/217x217/segmented/trainSet';

% storing over segmented examples
dataBlanksSeg = '/work/bioinformatics/shared/dope/torch/test/217x217/segmented/blanks'

randOrd = randperm(length(cellDataSet));
percentVal = .1;

parfor iR = 1:length(cellDataSet)
    i = randOrd(iR);
    MD = load(cellDataSet{i}.cellMD,'MD');
    %     I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
    MD = MD.MD;
    expStr = cellDataSet{i}.expStr;
    
    if rand > percentVal
        dataSetDir = 'trainSet';
    else
        dataSetDir = 'valSet';
    end

     
    for fidx = 1:MD.nFrames_
        I = mat2gray(MD.getChannel(1).loadImage(fidx));
        if resizeOn && size(I,1) ~= newDim(1)
            I = imresize(I, newDim);
        end
        
        if cellSegmentation
            [I mask] = segCellLCH(I,'preview',false);
            % check if 97% covering image, then
            rp = regionprops(mask);
            if rp.Area/size(I,1)^2  >= .9
                blank_mask = true;
            else
                blank_mask = false;
            end
            
        end
        
        frameNum = num2str(fidx);
        newFileOut = [MD.processes_{1}.funParams_.key '_f' frameNum '.png'];

        if MD.processes_{1}.funParams_.metEff == 0
            classDir = 'lowMet';
        elseif MD.processes_{1}.funParams_.metEff == 1
            classDir = 'highMet';
        else
            classDir = 'unMet';
        end
        
        dirOut = fullfile(dataRootDir, dataSetDir, classDir, expStr, newFileOut);
%         exrDir = fullfile(dirOut, newFileOut);
                           
        if blank_mask
            disp(['Blank mask: ' dirOut]);
            dirOut = fullfile(dataBlanksSeg, dataSetDir, classDir, expStr, newFileOut);
        end
        if exist(fileparts(dirOut),'dir') ~= 7
            mkdir(fileparts(dirOut));
        end
        imwrite(I, dirOut);
    end
    
end


%% Gather image data
% Need to figure out how to specify with subdirs correctly
% imdsTrain = imageDatastore(dataRootDirTR, 'IncludeSubfolders',true,'FileExtensions','.png','LabelSource','foldernames');
% imdsValid = imageDatastore(dataRootDirVal, 'IncludeSubfolders',true,'FileExtensions','.png','LabelSource','foldernames');

%% Gather image data
imdsTrain = imageDatastore(dataRootDirTR, 'IncludeSubfolders',true,'FileExtensions','.png');%'LabelSource','foldernames');
imdsValid = imageDatastore(dataRootDirVal, 'IncludeSubfolders',true,'FileExtensions','.png');%,'LabelSource','foldernames');

label_list_Train = cell(length(imdsTrain.Files),1);
label_list_Val = cell(length(imdsValid.Files),1);

label = 'highMet';
highMetLabel=cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false));
label_list_Train(highMetLabel) = {label};

label = 'unMet';
unMetLabel=cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false));
label_list_Train(unMetLabel) = {label};

label = 'lowMet';
lowMetLabel = cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false));
label_list_Train(lowMetLabel) = {label};
imdsTrain.Labels = categorical(label_list_Train);

label = 'highMet';
highMetLabel=cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false));
label_list_Val(highMetLabel) = {label};

label = 'unMet';
unMetLabel=cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false));
label_list_Val(unMetLabel) = {label};

label = 'lowMet';
lowMetLabel = cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false));
label_list_Val(lowMetLabel) = {label};
imdsValid.Labels = categorical(label_list_Val);


% label='lowMet';
% imdsTrain.Labels(cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false))) = label;
% label='unMet';
% imdsTrain.Labels(cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false))) = label;
% label='highMet';
% imdsTrain.Labels(cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false))) = label;
% 
% label='lowMet';
% imdsValid.Labels(cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false))) = label;
% label='unMet';
% imdsValid.Labels(cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false))) = label;
% label='highMet';
% imdsValid.Labels(cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false))) = label;

% get rid of unMet...
% [imds1] = splitEachLabel(imgs,1500,'Exclude','unMet');
% [trainData, valData] = splitEachLabel(imgs, 1500, 'randomize');

save(fullfile(dataRootDir, ['preTrainingDataStructure_217x217_segmented' datestr(datetime,'dd-mmm-yyyy-hhMM') '.mat']));