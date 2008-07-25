function [projList] = setupDirectories(selectROI,preserveROIs)
% SETUPDIRECTORIES creates/finds roi directories & returns searchable info
%
% DESCRIPTION: setupDirectories searches a user-selected, top-level
% directory for every sub-directory called /images, where presumably the
% original data for a movie resides.
%
% The function then creates a new directory at the same level called /rois,
% which will contain roi_1,...,roi_9 sub-directories.  Each roi_n directory
% will have an images and an analysis folder.  
% 
% If the user wants to select 1 or more rois (see 'selectROI'), (s)he can
% do so for up to nine per movie. Images will be cropped using a bounding
% rectangle and saved in roi_n/roi_images.  The selected ROI mask will be
% saved as a tif in roi_n/roi_analysis.  
%
% If the user does not want to select a ROI, the original images will
% simply be copied into the roi_n/roi_images folder.
%
% The user can either overwrite or retain /rois directories that have
% already been created by setupDirectories; in this way you can add more
% data and not undo previous analysis.  
%
% A cell array containing useful information will be the output (see
% 'projList' below).
%
%
% INPUT:
%
% selectROI:   1 if the user wants to select regions of interest (up to 9)
%              per image directory found in the user-selected top-level
%              directory; (0) if not (note: in this case, a rois/roi_1
%              directory will be created and the original images will be
%              copied into it.)
% preserveROI: (1) if projects that have already had rois selected should not
%              be overwritten; 0 if they should all be replaced
%
%
% OUTPUT:
%
% projList:         nSubprojects x 7 (or more) cell array of the form:
%                   [target oligo movie roi roi_imageDir roi_analysisDir
%                   origImDir (...)] 'projList' assumes that the data is in
%                   a directory hierarchy like this:
%                   target/oligo/movie/rois/roi_n, where roi_n is filled
%                   either by copying the original images or by querying
%                   the user to select a polygonal region of interest (see
%                   'selectROI' above).  The optional fields can contain
%                   data information or pointers to the results once the
%                   projects have been analyzed. 'projList' is saved in the
%                   top-level directory selected by the user.
% roiMask         : tif the size of the raw images containing the roiMask;
%                   this gets saved in the roi_n/roi_analysis folder
% roiMask_cropped : tif the size of the cropped images containing the roiMask;
%                   this gets saved in the roi_n/roi_analysis folder
%
% Created 20 July 2008 by Kathryn Applegate, Matlab R2008a

if nargin<2 || isempty(preserveROIs) || (preserveROIs~=0 && preserveROIs~=1)
preserveROIs=1;
end

if nargin<1 || isempty(selectROI) || (selectROI~=0 && selectROI~=1)
    selectROI=0;
end


topDir=uigetdir(pwd,'Please select top-level directory containing targets');

% imageDir is cell array containing all *\images directories in topDir
p=genpath(topDir);
imDir=strrep(p,';',' ');
[imageDir] = regexp(imDir,'\S*\\images\s','match')';

subProjCount=1; % counter for all rois from all projects
alreadyDoneList=0; % vector containing subProjCount index of projects that already have rois
for i=1:length(imageDir) % iterate through projects

    % this will be changed to 1 for projects that already have roi_n directories
    alreadyDoneFlag=0;

    % define image and rois directories
    runInfo(i,1).imDir=imageDir{i,1}(1:end-1);
    runInfo(i,1).roiDir=[runInfo(i,1).imDir(1:end-7) filesep 'rois'];

    % check for existence of rois directory; make dir if needed
    if isdir(runInfo(i,1).roiDir) && preserveROIs==1
        % there is a roi directory already and we want to preserve it
        alreadyDoneFlag=1;
        roiCount=[];
        % count how many roi_n directories there are
        for iRoi=1:9 % max rois per movie is 9
            if isdir([runInfo(i,1).roiDir filesep 'roi_' num2str(iRoi)])
                roiCount=[roiCount iRoi];
            end
        end

    else
        if isdir(runInfo(i,1).roiDir)
            rmdir(runInfo(i,1).roiDir,'s');
        end
        mkdir(runInfo(i,1).roiDir);
    end

    if alreadyDoneFlag==0 % need to create 1 or more rois

        % get list and number of images
        [listOfImages]=searchFiles('.tif',[],runInfo(i,1).imDir,0);

        % read first image
        img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));

        roiCount=1; % counter for rois for current project
        makeNewROI=1; % flag for making new roi

        % iterate til the user is finished or just copy if not choosing roi
        while makeNewROI==1 && roiCount<10
            % make new roi_n image/analysis directories
            currentRoiImDir=[runInfo(i,1).roiDir filesep 'roi_' num2str(roiCount) filesep 'roi_images'];
            mkdir(currentRoiImDir);
            currentRoiAnDir=[runInfo(i,1).roiDir filesep 'roi_' num2str(roiCount) filesep 'roi_analysis'];
            mkdir(currentRoiAnDir);

            if selectROI==1
                % draw polygon to make mask for bounding tracks
                [roiMask]=roipoly((img-min(img(:)))./(max(img(:))-min(img(:))));
                close

                % make bigger for edges of gaussian
                roiMaskDil=bwmorph(roiMask,'dilate',6);
                % get boundary box parameters for roiMask from dilated
                [rows,cols]=find(roiMaskDil);
                ymin=min(rows); height=max(rows)-min(rows)+1;
                xmin=min(cols); width=max(cols)-min(cols)+1;

                % crop and save ROI
                roiMask=logical(roiMask);
                croppedRoi=imcrop(roiMask,[xmin ymin width height]);
            else
                % make the ROI the whole image
                roiMask=ones(size(img));
                croppedRoi=roiMask;
            end

            % save original and croppsed roiMask
            imwrite(roiMask,[currentRoiAnDir filesep 'roiMask.tif']);
            imwrite(croppedRoi,[currentRoiAnDir filesep 'roiMask_cropped.tif']);

            % parse the path to get "words" used to identify target, oligo,
            % movie, and roi
            nChar=length(currentRoiImDir);
            filesepLoc=regexp(currentRoiImDir,'\\');
            wordStart=[1 filesepLoc+1]; wordEnd=[filesepLoc-1 nChar];
            words=cell(length(wordStart),1);
            for iWord=1:length(wordStart)
                words{iWord,1}=currentRoiImDir(wordStart(iWord):wordEnd(iWord));
            end

            % assign cell to be nProj x 7
            projList{subProjCount,1}=words{end-5,1};     % target
            projList{subProjCount,2}=words{end-4,1};     % oligo
            projList{subProjCount,3}=words{end-3,1};     % movie
            projList{subProjCount,4}=words{end-1,1};     % roi
            projList{subProjCount,5}=currentRoiImDir;    % roiImageDir
            projList{subProjCount,6}=currentRoiAnDir;    % roiAnalysisDir
            projList{subProjCount,7}=runInfo(i,1).imDir; % originalImageDir

            if selectROI==1
                reply=input('Do you want to select another ROI? y/n [n]: ','s');
            else
                reply='n';
            end
            if lower(reply)=='y'
                makeNewROI=1; % user said yes; make another one
                roiCount=roiCount+1; % counter for current condition rois
            else
                makeNewROI=0; % assume no; we're done
            end

            subProjCount=subProjCount+1; % counter for all conditions in top directory
        end

    elseif alreadyDoneFlag==1

        for iRoi=1:length(roiCount)
            % these are the current roi's im/an directories
            currentRoiImDir=[runInfo(i,1).roiDir filesep 'roi_' num2str(roiCount(iRoi)) filesep 'roi_images'];
            currentRoiAnDir=[runInfo(i,1).roiDir filesep 'roi_' num2str(roiCount(iRoi)) filesep 'roi_analysis'];

            % parse the path to get "words" used to identify target, oligo,
            % movie, and roi
            nChar=length(currentRoiImDir);
            filesepLoc=regexp(currentRoiImDir,'\\');
            wordStart=[1 filesepLoc+1]; wordEnd=[filesepLoc-1 nChar];
            words=cell(length(wordStart),1);
            for iWord=1:length(wordStart)
                words{iWord,1}=currentRoiImDir(wordStart(iWord):wordEnd(iWord));
            end

            % assign cell to be nProj x 7
            projList{subProjCount,1}=words{end-5,1};     % target
            projList{subProjCount,2}=words{end-4,1};     % oligo
            projList{subProjCount,3}=words{end-3,1};     % movie
            projList{subProjCount,4}=words{end-1,1};     % roi
            projList{subProjCount,5}=currentRoiImDir;    % roiImageDir
            projList{subProjCount,6}=currentRoiAnDir;    % roiAnalysisDir
            projList{subProjCount,7}=runInfo(i,1).imDir; % originalImageDir

            alreadyDoneList=[alreadyDoneList subProjCount];
            subProjCount=subProjCount+1; % counter for all conditions in top directory

        end

    end
end

for i=1:subProjCount-1
    
    % only do the ones that haven't been done
    if ~ismember(i,alreadyDoneList) 

        % get list and number of images from original data
        [listOfImages]=searchFiles('.tif',[],projList{i,7},0);
        nImages=size(listOfImages,1);
        s1=length(num2str(nImages));
        strg1=sprintf('%%.%dd',s1);

        % read first image, see if 8-bit or 16-bit
        img=imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]);
        if isa(img,'uint8')
            bitType=8;
        elseif isa(img,'uint16')
            bitType=16;
        else
            error('images should be 8-bit or 16-bit')
        end
        img=double(img);

        if selectROI==1 % CROP images from originals
            roiMask=double(imread([projList{i,6} filesep 'roiMask.tif']));
            % make bigger for edges of gaussian
            roiMaskDil=bwmorph(roiMask,'dilate',6);
            % get boundary box parameters for roiMask from dilated
            [rows,cols]=find(roiMaskDil);
            ymin=min(rows); height=max(rows)-min(rows)+1;
            xmin=min(cols); width=max(cols)-min(cols)+1;
            % crop and save images
            for iImage=1:size(listOfImages,1)
                img=double(imread([char(listOfImages(iImage,2)) filesep char(listOfImages(iImage,1))]));

                if bitType==8
                    croppedImg=uint8(imcrop(img,[xmin ymin width height]));
                else
                    croppedImg=uint16(imcrop(img,[xmin ymin width height]));
                end

                indxStr1=sprintf(strg1,iImage);
                imwrite(croppedImg,[projList{i,5} filesep 'image_cropped' indxStr1 '.tif']);
            end
        else % COPY images from originals
            copyfile(projList{i,7},projList{i,5});
        end
    end
end

save([topDir filesep 'projList'],'projList')



