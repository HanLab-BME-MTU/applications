function setupRoiDirectories(selectROI,overwriteROIs,doCrop)
% SETUPROIDIRECTORIES allows quick roi directory creation for tif series
%
% SYNOPSIS: setupRoiDirectories(selectROI,overwriteROIs,doCrop)
%
% DESCRIPTION: setupRoiDirectories searches a user-selected parent
% directory for every sub-directory called /images, where the
% original data for a movie resides.
%
% The function then creates up to 9 new directories at the
% same level, called roi_1, roi_2,...
%
% The user can either overwrite or retain roi directories that have
% already been created.
%
% INPUT  selectROI     : 0 if using the whole image as a ROI
%                        1 if the user wants to select regions of interest
%                        for each image directory found in the 
%                        user-selected top-level directory
%                        2 if the user wants to select the cell centers to 
%                        create background region of interest
%        overwriteROIs : 1 if projects that have already have rois should
%                        be overwritten; (0) if they should be preserved
%                        and more should be added
%        doCrop (opt)  : 1 if you want to generate cropped 16-bit images in
%                        a subdirectory of roi_x called 'cropped_imgs'.
%
%
% Created 20 July 2008 by Kathryn Applegate, Matlab R2008a

topDir=uigetdir(pwd,'Please select parent directory containing one or more TIFF series.');
if topDir==0
    return
end

% default - make roi_1 directory, roi is whole image
if nargin<1 || isempty(selectROI)
    selectROI=0;
end

% default - don't overwrite old data - just look for new projects under
% parent directory
if nargin<2 || isempty(overwriteROIs) || (overwriteROIs~=0 && overwriteROIs~=1)
    overwriteROIs=0;
end

reply='no';
if overwriteROIs==1
    reply = questdlg('All roi_x directories under the parent directory you just selected will be overwritten. Do you wish to continue?');
    if strcmpi(reply,'cancel') || strcmpi(reply,'no')
        return
    end
end
if strcmpi(reply,'yes')
    overwriteROIs=1;
end

% default - don't crop the images and store them
if nargin<3 || isempty(doCrop)
    doCrop=0;
end
% assume if we're cropping, we need to select a region to crop
if doCrop==1
    selectROI=1;
end

% find existing /images and /roi_x directories. "images" directories make
% up the list of all the projects; "roi_x" ones have already been analyzed
dirList=regexp(genpath(topDir),pathsep,'split');
% cell array of "images" directories
imageDirList=regexp(dirList,'(.+)images$','match','once');
imageDirList=imageDirList(~cellfun(@isempty,imageDirList));
% check whether there are any movies in the image directory list
if isempty(imageDirList)
    h=msgbox('No directories called "images" were found.');
    uiwait(h);
    return
end

% if we should overwrite ROI data, remove existing directories
if overwriteROIs==1
    % cell array of "roi_x" directories
    roiDirList  =regexp(dirList,'(.+)roi_(\d+)$','match','once')';
    roiDirList=roiDirList(~cellfun(@isempty,roiDirList));
    for iProj=1:length(roiDirList)
        temp=roiDirList{iProj}(1:end-1);
        rmdir(temp,'s');
    end
end

% flags for help messages
crpFlag=doCrop;
polyFlag=selectROI;
bgFlag=logical(selectROI);

% loop through projects
for iProj=1:length(imageDirList)

    % define image and roi directories
    imDir=imageDirList{iProj};
    roiDir=[imDir(1:end-7) filesep 'roi'];

    % count existing rois - roiCount will be number of first new one
    roiCount=1;
    flag=0;
    while flag==0
        if isdir([roiDir '_' num2str(roiCount)])
            roiCount=roiCount+1;
        else
            flag=1;
        end
    end

    % get list and number of images
    [listOfImages]=searchFiles('.tif',[],imDir,0);

    % make new ROIs until user is finished or we reach 10
    makeNewROI=1;
    while makeNewROI==1 && roiCount<10

        % read first image
        img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
        img=(img-min(img(:)))./(max(img(:))-min(img(:)));

        % make new roi_x analysis directory
        anDir=[roiDir '_' num2str(roiCount)];
        mkdir(anDir);

        if crpFlag==1
            h=msgbox('FYI: Cropped images and the cropping mask and coordinates will be stored in roi folder.');
            uiwait(h);
            crpFlag=0;
        end

        if selectROI==1

            if polyFlag==1
                h=msgbox('Select a polygon, right-click on last point, and select "Create Mask."');
                uiwait(h);
                polyFlag=0;
            end

            roiMask=[];
            while isempty(roiMask)
                try
                    % draw polygon to make mask
                    figure;
                    [roiMask,polyXcoord,polyYcoord]=roipoly(img);
                    roiYX=[polyYcoord polyXcoord; polyYcoord(1) polyXcoord(1)];
                    close(gcf)
                catch
                    h=msgbox('Please try again.','help');
                    uiwait(h);
                end
            end
        else
            % make the ROI the whole image
            [imL,imW]=size(img);
            roiMask=true(imL,imW);
            roiYX=[1 1; imL 1; imL imW; 1 imW; 1 1];
        end
        
        if selectROI==2
            % get within-cell point for each project
            figure; imagesc(img); colormap gray; axis equal
            if bgFlag
                h=msgbox('For each ROI, select a point in the cell center and press enter.');
                uiwait(h);
                bgFlag=0;
            end
            [x,y] = getpts;
            bgPtYX = [y(1) x(1)];

            save([anDir filesep 'bgPtYX'],'bgPtYX');
            close(gcf)
        end

        % save roiMask and coordinates
        imwrite(roiMask,[anDir filesep 'roiMask.tif']);
        save([anDir filesep 'roiYX'],'roiYX');
        disp(['Created: ' anDir])

        if doCrop==1

            minY=floor(min(roiYX(:,1)));
            maxY=ceil(max(roiYX(:,1)));
            minX=floor(min(roiYX(:,2)));
            maxX=ceil(max(roiYX(:,2)));

            mkdir([anDir filesep 'cropped_imgs']);
            for j=1:size(listOfImages,1)
                imgName=listOfImages{j,1};
                img=double(imread([listOfImages{j,2} filesep imgName]));
                img=uint16(img(minY:maxY,minX:maxX));
                imwrite(img,[anDir filesep 'cropped_imgs' filesep imgName]);
            end
        end

        if selectROI==1
            reply = questdlg('Do you want to select another ROI?');
        else
            reply='no';
        end
        if strcmpi(reply,'yes')
            makeNewROI=1; % user said yes; make another one
            roiCount=roiCount+1; % counter for current condition rois
        else
            makeNewROI=0; % assume no; we're done
        end
    end % while makeNewROI==1 && roiCount<10
end % iterate through projects

cd(topDir)
getProj(pwd)


