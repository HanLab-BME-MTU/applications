function segmentationParam = cellMasksSteerableFilterPlusManual(...
    firstImageFile,dirMask,segmentationParam)
%CELLMASKSSTEERABLEFILTERPLUSMANUAL is an interactive tool that allows the user to find cell edges based on intensity gradients
%
%SYNOPSIS segmentationParam = cellMasksSteerableFilterPlusManual(...
%    firstImageFile,dirMask,segmentationParam)
%
%INPUT  firstImageFile: Full path and name of first file in image stack.
%       dirMask       : Directory name for storing masks.
%       segmentationParam: Same as output, in case movie has been
%                       previously analyzed and only a subset is to be
%                       done/redone. The code looks at the field "success"
%                       (as defined below) and segments again images with
%                       success = 0.
%
%OUTPUT segmentationParam: Structure array with number of entries = number
%                       of images, containing the following fields:
%           .filterSigma     : Gaussian sigma used for line filter.
%           .closureRadius   : Radius used to close edge gaps.
%           .success         : 1 to indicate successful segmentation, 0
%                              otherwise.
%
%Khuloud Jaqaman, October 2011

%% Input

%get image information
[fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
dirImage=[fpath,filesep];
fNameImage=[fname,fno,fext];

%get mask information
dirMask=[dirMask,filesep];
fNameMask = ['mask_' fname];

%get all files in image stack
inFileList = getFileStackNames([dirImage,fNameImage]);
numFiles = length(inFileList);

%reserve memory for structure storing analysis parameters if it was not
%input
if nargin < 3 || isempty(segmentationParam)
    segmentationParam = repmat(struct('filterSigma',[],'closureRadius',[],'success',0),numFiles,1);
end

%find images which must be segmented
successVec = vertcat(segmentationParam.success);
indx2analyze = find(successVec==0)';

%% Masks

%analyze each image
for iFile = indx2analyze
    
    disp(' ');
    disp(['Image ' num2str(iFile) ' of ' num2str(numFiles)]);
    disp(' ');
    
    try
        
        %read image
        image = imread(inFileList{iFile});
        
        %generate mask
        genMask = 1;
        while genMask
            [mask,filterSigma,closureRadius] = steerableFilterPlusManual(image);
            disp(' ');
            userEntry = input('Success? y/n ','s');
            segmentationParam(iFile).success = strcmp(userEntry,'y');
            if strcmp(userEntry,'n')
                disp(' ');
                userEntry = input('Repeat segmentation for this image? "y" to repeat, "n" to go to next image. ','s');
                genMask = strcmp(userEntry,'y');
            else
                genMask = 0;
            end
            close all
        end
        
        %save mask
        [~,~,fno,fext]=getFilenameBody(inFileList{iFile});
        filename = [dirMask,fNameMask,fno,fext];
        imwrite(mask,filename);
        
        %store segmentation parameters
        segmentationParam(iFile).filterSigma = filterSigma;
        segmentationParam(iFile).closureRadius = closureRadius;
        
    catch %#ok<CTCH>
        
        disp(['image ' num2str(iFile) ' failed']);
        segmentationParam(iFile).success = 0;
        
    end
    
    close all
    
end

save(fullfile(dirMask,['segParam_' fname]),'segmentationParam');

%% ~~~ the end ~~~
