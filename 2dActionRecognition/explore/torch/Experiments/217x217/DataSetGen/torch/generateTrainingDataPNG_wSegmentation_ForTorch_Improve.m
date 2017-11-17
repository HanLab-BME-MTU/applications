% simple script to pre-process the cells for Deep learning.
clc;
clear;
%% Gen Data
load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');

resizeOn = false;
padImage = false;%true;
cellSegmentation = true;
newDim = [256 256];
center_mass = true;
variance_filter = true;
cropImage = true;
cropSize = [150 150];

% dataRootDir = '/work/bioinformatics/shared/dope/torch/test/AAE/images/var/217x217_test/';
dataRootDir = '/work/bioinformatics/shared/dope/torch/test/AAE/images/var/150x150/';
% dataRootDir = '/work/bioinformatics/shared/dope/torch/test/AAE/images/64x64/';
dataRootDirVal = fullfile(dataRootDir,'test');
dataRootDirTR = fullfile(dataRootDir,'train');

% storing over segmented examples
dataBlanksSeg = fullfile(dataRootDir,'blanks');

randOrd = randperm(length(cellDataSet));
percentVal = .05;

parfor iR = 1:10%length(cellDataSet)
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
            
            [I mask] = segCellLCH(I,'preview', false, 'varFilterOut', variance_filter, ...
                                    'centerMass', center_mass);

            % check if 97% covering image, then
            rp = regionprops(mask);
            Area = rp.Area;
            if Area/size(I,1)^2  >= .9
                blank_mask = true;
            else
                blank_mask = false;
            end
            
        end
        
        if cropImage 
            x1 = ceil((size(I,1)-cropSize(1))/2);
            x2 = cropSize(1)+x1; 
            I =  I(x1:x2,x1:x2);
            I = I(1:cropSize(1),1:cropSize(2));
        end
        
        if padImage
            padSize=round((newDim(1) - size(I,1))/2)+1;
            I = padarray(I,[padSize padSize]);
            I = I(1:newDim(1), 1:newDim(2));
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