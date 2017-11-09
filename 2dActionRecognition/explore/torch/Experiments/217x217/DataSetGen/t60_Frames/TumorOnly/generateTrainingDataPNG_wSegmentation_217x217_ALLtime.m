% simple script to pre-process the cells for Deep learning.
clc;
clear;
%% Gen Data
% load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');
load('/work/bioinformatics/shared/dope/torch/test/cellData/cellMovieData_wAnno08May2017_0558.mat', 'cellDataSet')
% cellXploreDR('/work/bioinformatics/shared/dope/torch/test/cellData/cellMovieData_wAnno08May2017_0558.mat');
cellSegmentation = true;
newDim = [217 217];
resizeOn = false;

% First, build the tumor one....
%  --- > /work/bioinformatics/shared/dope/data/OMETIFF/all/sdcAll


dataRootDir = '/work/bioinformatics/shared/dope/torch/test/217x217/allTime';
dataRootDirVal = fullfile(dataRootDir, 'valSet');
dataRootDirTR = fullfile(dataRootDir, 'trainSet');

% storing over segmented examples
dataBlanksSeg = '/work/bioinformatics/shared/dope/torch/test/217x217/allTime/blanks';

randOrd = randperm(length(cellDataSet));
percentVal = .1;

parfor iR = 1:length(cellDataSet)
    i = randOrd(iR);
    javaaddpath('/home2/s170480/matlab/extern/bioformats/bioformats_package.jar','-end')
    try 
        MD = load(cellDataSet{i}.cellMD,'MD');
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
                if rp.Area/size(I,1)^2  >= .85
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

            if blank_mask
                disp(['Blank mask: ' dirOut]);
                dirOut = fullfile(dataBlanksSeg, dataSetDir, classDir, expStr, newFileOut);
            end
            if exist(fileparts(dirOut),'dir') ~= 7
                mkdir(fileparts(dirOut));
            end
            imwrite(I, dirOut);
        end
    catch ME
        disp(['FAILED MD : ' MD.movieDataPath_]);  
        disp(ME)
        disp(ME.stack)
    end
end


%% Gather image data
imdsTrain = imageDatastore(dataRootDirTR, 'IncludeSubfolders',true,'FileExtensions','.png');%'LabelSource','foldernames');
imdsValid = imageDatastore(dataRootDirVal, 'IncludeSubfolders',true,'FileExtensions','.png');%,'LabelSource','foldernames');


label_list_Train = cell(length(imdsTrain.Files),1);
label_list_Val = cell(length(imdsValid.Files),1);

label = 'highMet';
highMetLabel=cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false));
label_list_Train(highMetLabel) = {label};

label = 'lowMet';
lowMetLabel = cell2mat(cellfun(@(x) contains(x,label) ,imdsTrain.Files, 'Uniform', false));
label_list_Train(lowMetLabel) = {label};
imdsTrain.Labels = categorical(label_list_Train);

label = 'highMet';
highMetLabel=cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false));
label_list_Val(highMetLabel) = {label};

label = 'lowMet';
lowMetLabel = cell2mat(cellfun(@(x) contains(x,label) ,imdsValid.Files, 'Uniform', false));
label_list_Val(lowMetLabel) = {label};
imdsValid.Labels = categorical(label_list_Val);


uisave();

% save(fullfile(dataRootDir, ['preTrainingDataStructure_217x217_segmented' datestr(datetime,'dd-mmm-yyyy-hhMM') '.mat']));