% simple script to pre-process the cells for Deep learning.
clc;
clear;
%% Gen Data
load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');

resizeOn = true;
cellSegmentation = true;
newDim = [128 128];

dataRootDir = '/work/bioinformatics/shared/dope/torch/test/AAE/images/128x128/';
dataRootDirVal = fullfile(dataRootDir,'test');
dataRootDirTR = fullfile(dataRootDir,'train');

% storing over segmented examples
dataBlanksSeg = fullfile(dataRootDir,'blanks');

randOrd = randperm(length(cellDataSet));
percentVal = .1;

parfor iR = 1:length(cellDataSet)
    i = randOrd(iR);
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
            if rp.Area/size(I,1)^2  >= .9
                blank_mask = true;
            else
                blank_mask = false;
            end
            
        end
        
        frameNum = num2str(fidx);
        newFileOut = [cellDataSet{i}.key '_f' frameNum '.png'];

        if cellDataSet{i}.metEff == 0
            classDir = 'lowMet';
        elseif cellDataSet{i}.metEff == 1
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