
% simple script to pre-process the cells for Deep learning.
% Andrew R. Jamieson Oct 2017


% load data state 
clear;
clc;
% [filename pathname] = uigetfile('Which training set to use?');
% load(fullfile(pathname,filename));
[dataRootDirTR] = uigetdir('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/',...
                           'Which dir has train data?');
[dataRootDirVal] = uigetdir('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/',...
                            'Which dir has test data?');

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


uisave({'imdsTrain','imdsValid'},['preTrainingDataStructure_217x217_segmented' datestr(datetime,'dd-mmm-yyyy-hhMM') '.mat']);
% save(fullfile(dataRootDir, ['preTrainingDataStructure_217x217_segmented' datestr(datetime,'dd-mmm-yyyy-hhMM') '.mat']));
