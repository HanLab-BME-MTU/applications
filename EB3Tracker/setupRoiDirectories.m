function setupRoiDirectories(selectROI,overwriteROIs,doCrop)
% SETUPROIDIRECTORIES allows quick roi directory creation for movie(s)
%
% DESCRIPTION: setupRoiDirectories searches a user-selected, top-level
% directory for every sub-directory called /images, where presumably the
% original data for a movie resides.
%
% The function then creates one or more (up to 9) new directories at the 
% same level, called roi_1,...,roi_9, depending on how many ROIs the user
% selects.
%
% The user can either overwrite or retain roi directories that have
% already been created; in this way you can add more data to the parent
% directory and select ROIs for just those new movies (leaving the old ones
% intact), or you may want to repeat ROI selection for all the movies.
%
% INPUT  selectROI     : 1 if the user wants to select regions of interest
%                        per image directory found in the user-selected top-level
%                        directory
%                        (0) if not (note: in this case, a roi_1 directory
%                        will be created and the "roi" will be the whole image)
%        overwriteROIs : 1 if projects that have already have rois should
%                        be overwritten 
%                        (0) if they should be preserved (in this case you
%                        cannot add more ROIs to the old projects - just to
%                        the new ones.  See note above.)
%        doCrop        : optional parameter. 1 if you want to generate
%                        cropped 16-bit images in a subdirectory of roi_x
%                        called images.
%
% OUTPUT roiMask       : tif the size of the raw images containing the roiMask;
%                        this gets saved in the roi_n folder
%        roiYX         : contains xy-coordinates of the roi
%
% Created 20 July 2008 by Kathryn Applegate, Matlab R2008a

topDir=uigetdir(pwd,'Please select top-level directory containing targets');
if topDir==0
    return
end

% default - make roi_1 directory, roi is whole image
if nargin<1 || isempty(selectROI) || (selectROI~=0 && selectROI~=1)
    selectROI=0;
end

% default - don't overwrite old data - just look for new projects under
% parent directory
if nargin<2 || isempty(overwriteROIs) || (overwriteROIs~=0 && overwriteROIs~=1)
    overwriteROIs=0;
end

reply='no';
if overwriteROIs==1
    reply = questdlg('This function will overwrite all existing ROI directories under the one you just selected. Do you wish to continue?');
    if strcmpi(reply,'cancel')
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

% find existing /images and /roi_x directories. "images" directories make
% up the list of all the projects; "roi_x" ones have already been analyzed

p=genpath(topDir);
if ispc
    tempDirList=strrep(p,';',' ');
else
    tempDirList=strrep(p,':',' ');
end
imageDirList=regexp(tempDirList,['\S*images\s'],'match')'; % cell array of "images" directories
roiDirList  =regexp(tempDirList,['\S*roi_\d\s'],'match')'; % cell array of "roi_x" directories


% if we should overwrite ROI data, remove it
if overwriteROIs==1 && ~isempty(roiDirList)
    for i=1:length(roiDirList)
        temp=roiDirList{i,1}; temp=temp(1:end-1);
        rmdir(temp,'s');
    end
end

if isequal(size(imageDirList),[0 0])
    h=msgbox('No directories called "images" were found.');
    uiwait(h);
end

roiCount=0;
for i=1:length(imageDirList) % iterate through projects

    % define image and roi directories
    imDir=imageDirList{i,1}(1:end-1);
    roiDir=[imDir(1:end-7) filesep 'roi'];

    % check for existence of rois directory and make dir if needed
    if ~isdir([roiDir '_1'])
        % get list and number of images
        [listOfImages]=searchFiles('.tif',[],imDir,0);

        % read first image
        img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
        img=(img-min(img(:)))./(max(img(:))-min(img(:)));
        roiCount=1; % counter for rois for current project
        makeNewROI=1; % flag for making new roi

        % iterate til the user is finished
        while makeNewROI==1 && roiCount<10
                % make new roi_n image/analysis directories
                currentRoiAnDir=[roiDir '_' num2str(roiCount)];
                mkdir(currentRoiAnDir);
                
                if selectROI==1
                    roiMask=[];
                    while isempty(roiMask)
                        try
                            % draw polygon to make mask
                            figure
                            [roiMask,polyXcoord,polyYcoord]=roipoly(img);
                        catch
                            h=msgbox('Please try again.','help');
                            uiwait(h);
                        end
                    end
                    close
                    roiYX=[polyYcoord polyXcoord; polyYcoord(1) polyXcoord(1)];
                else
                    % make the ROI the whole image
                    [imL,imW]=size(img);
                    roiMask=logical(ones(imL,imW));
                    roiYX=[1 1; imL 1; imL imW; 1 imW; 1 1];
                    
                    % get within-cell point for each project
                    figure
                    imagesc(img); colormap gray;
                    if i==1
                        h=msgbox('For each ROI, select a point in the cell center and press enter.');
                        uiwait(h);
                    end
                    [x,y] = getpts;
                    bgPtYX = [y(1) x(1)];
                    
                    save([currentRoiAnDir filesep 'bgPtYX'],'bgPtYX');
                    close
                end

                % save roiMask and coordinates
                imwrite(roiMask,[currentRoiAnDir filesep 'roiMask.tif']);
                save([currentRoiAnDir filesep 'roiYX'],'roiYX');
                disp(currentRoiAnDir)
                
                if doCrop==1
                    minY=floor(min(roiYX(:,1)));
                    maxY=ceil(max(roiYX(:,1)));
                    minX=floor(min(roiYX(:,2)));
                    maxX=ceil(max(roiYX(:,2)));
                    
                    mkdir([currentRoiAnDir filesep 'images']);
                    for j=1:size(listOfImages,1)
                        imgName=[char(listOfImages(j,2)) filesep char(listOfImages(j,1))];
                        img=double(imread(imgName));
                        img=uint16(img(minY:maxY,minX:maxX));
                        imwrite(img,[currentRoiAnDir filesep 'images' filesep char(listOfImages(j,1))]);
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
    else % if there's no roi_1 directory
        
    end
end % iterate through projects