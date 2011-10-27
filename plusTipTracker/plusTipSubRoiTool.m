function plusTipSubRoiTool(projList,selectType,distUnits,distVal,timeUnits,timeVal,cellRoiYX,pickExclude)
% plusTipSubRoiTool allows sub-ROI selection and extracts MT growth tracks
%
% SYNOPSIS: plusTipSubRoiTool(projList,selectType,distUnits,distVal,timeUnits,timeVal,cellRoiYX,pickExclude)
%
% selectType:
%   0: manual, freehand regions
%   1: auto, cell body and peripheral region
%   2: auto, cell body and four peripheral regions
% distUnits/Val control how the division is made between the cell boundary
%               and the cell body
%       distUnits: {'fraction','microns'}
%       distVal  : if distUnits is fraction: between 0 and 1
%                  if distUnits is microns : >0
% timeUnits/Val control how within-roi tracks are defined
%       timeUnits: {'fraction','frames'}
%       timeVal  : if distUnits is fraction: between 0 and 1
%                  if distUnits is seconds : >0
%% OPTION TO TURN ON MICROPATTERN
micropattern = 0;  % if set to 1 call micropattern function and ignore rest
%selectType = 'cellPeriphSingle';
useSegMask = 0; % option to use a previously defined segmented Mask 
% could change to make it read in the the name of the mask
subRoiFilename = 'subRois';


%%
if micropattern == 1
   plusTipSubRoiToolMicropatterns(projList,selectType,distUnits,distVal,timeUnits,timeVal,cellRoiYX,pickExclude);
else % proceed with Kathryn's function
%% Check Input

homeDir=pwd;
warningState = warning;
warning('off','Images:initSize:adjustingMag')
warning('off','MATLAB:divideByZero')



if ispc
    fileExt='.emf';
else
    fileExt='.tif';
end

if nargin<1 || isempty(projList)
    projList=combineProjListFiles(0);
    if isempty(projList)
        return
    end

    a=projList2Cell(projList);
    a=a(:,1);

    % allow multiple projects to be selected
    if isempty(a)
        selection=[];
    else
        [selection,selectionList]=listSelectGUI(a,[],'move',1);
    end

    if isempty(selection)
        msgbox('plusTipSubRoiTool: No projects selected or tracking has not been completed.')
        return
    else
        projList=projList(selection,1);
    end
end


if nargin<2 || isempty(selectType)
    selectType=0; % default is manual
end

if nargin<3 || isempty(distUnits)
    distUnits='fraction';
end
if nargin<4 || isempty(distVal)
    distVal=.5;
end
if nargin<5 || isempty(timeUnits)
    timeUnits='fraction';
end
if nargin<6 || isempty(timeVal)
    timeVal=[];
end
if nargin<7 || isempty(cellRoiYX)
    cellRoiYX=[];
end
if nargin<8 || isempty(pickExclude)
    pickExcludeInput=0;
else
    pickExcludeInput=pickExclude;
end


switch selectType
    case 0
        selectType='manual';
    case 1
        selectType='cellPeriphSingle';
    case 2 
        selectType='cellPeriphQuad';
end

if ~ismember(lower(distUnits),{'fraction','microns'})
    error('plusTipSubRoiTool: distUnits must be fraction or microns')
end
if ~isempty(strmatch(lower(distUnits),'fraction')) && ~(distVal>0 && distVal<1)
    error('plusTipSubRoiTool: distUnits is fraction, distVal must be in 0-1')
end

if ~ismember(lower(timeUnits),{'fraction','seconds'})
    error('plusTipSubRoiTool: timeUnits must be fraction or seconds')
end
if ~isempty(strmatch(lower(timeUnits),'fraction')) && ~(timeVal>0 && timeVal<=1)
    error('plusTipSubRoiTool: timeUnits is fraction, timeVal must be in 0-1')
end
%% Body
collectPlots = 1; 
nProj=length(projList);
if collectPlots == 1
up1 = getFileNameBody(projList(1,1).anDir);
collectedDataPath = getFileNameBody(up1); 
mkdir([collectedDataPath filesep 'collectedSubRoiPlots'])
end 

for iProj=1:nProj

    anDir=projList(iProj,1).anDir;

    % check if current project is a sub-roi.  if it is, we extract tracks
    % only, since selecting sub-rois recursively is not allowed.  if it is
    % a full roi, we go on with roi selection.
    if ~isempty(strfind(anDir,'sub'))
        continue
    end

    
    subanDir=[anDir filesep subRoiFilename];
    if isdir(subanDir)  
        rmdir(subanDir,'s');  
        
    end
    mkdir(subanDir);
    cd(anDir)

    % load projData
    projData=load([anDir filesep 'meta' filesep 'projData']);
    projData=projData.projData;
    imDir=projData.imDir;

    % get list of images, read first image, and make RGB (gray)
    [listOfImages]=searchFiles('.tif',[],formatPath(imDir),0);
    img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
    [imL,imW]=size(img(:,:,1));
    img1=(img-min(img(:)))./(max(img(:))-min(img(:)));
    img=repmat(img1,[1,1,3]);

    cMap=hsv(10);



    % make sub-roi directory under roi_x directory and go there
    flag=0;
    roiCount=1; % counter for rois for current project
    img2show=img;
    while flag==0
        currentRoiAnDir=[subanDir filesep 'sub_' num2str(roiCount)];
        if isdir(currentRoiAnDir)
            % add the current roi to the composite image
            tempRoi=imread([currentRoiAnDir filesep 'roiMask.tif']);
            [img2show]=addMaskInColor(img2show,tempRoi,[1 1 1]);
            roiCount=roiCount+1;
        else
            flag=1;
        end
    end
    roiStart=roiCount;

    % if single project and roi coords loaded, get roiMask
    roiMask=[];
    if nProj==1 && ~isempty(cellRoiYX)
        if ~isequal(size(cellRoiYX),size(img(:,:,1)))
            roiYX=cellRoiYX;
            roiMask=roipoly(img,roiYX(:,2),roiYX(:,1));
        end
    end
    % otherwise let user load a mask or create a new one if this is the
    % first subroi, or load the already-saved one from subanDir
    if isempty(roiMask)
        if roiCount>1 % for sub_2 or later
            % load the previously-used roiMask
            roiMask=imread([subanDir filesep 'roiMask.tif']);
            roiYX=load([subanDir filesep 'roiYX.mat']); roiYX=roiYX.roiYX;
        else % for sub_1, establish what cell region will be
            if useSegMask ~= 1 
            choice=questdlg('Before creating Sub-ROIs, you need to define the cell boundary.','Cell ROI option','Draw new','Load roiYX.mat','Draw new');
            if ~isempty(strmatch(choice,'Draw new'))
                c=1;
                while isempty(roiMask)
                    try
                        figure
                        [roiMask,polyXcoord,polyYcoord]=roipoly(img);
                        close
                        roiYX=[polyYcoord polyXcoord; polyYcoord(1) polyXcoord(1)];
                    catch
                        h=msgbox('Please try again.','help');
                        uiwait(h);
                        c=c+1;
                        if c<=3 % give them 3 chances, then abort
                            return
                        end
                    end
                end
            else % load from file
                c=1;
                while isempty(roiMask)
                    try
                        [FileName,PathName] = uigetfile({'*.*'},'Select roiYX.mat');
                        p=load([PathName FileName]);
                        roiYX=p.roiYX;
                        roiMask=roipoly(img,roiYX(:,2),roiYX(:,1));
                    catch
                        h=msgbox('Please try again.','help');
                        uiwait(h);
                        c=c+1;
                        if c<=3 % give them 3 chances, then abort
                            return
                        end
                    end
                end
            end
            else 
                 p = load([anDir filesep 'roiYXSeg.mat']); % load the roiYX file in the anDir
                        roiYX=p.roiYX;
                        roiMask=roipoly(img,roiYX(:,2),roiYX(:,1)); 
                
                
        
            end
        end 
    end 

    % set cell boundary in white to composite image
    [img2show]=addMaskInColor(img2show,roiMask,[1 1 1]);

    % make distance transform before getting areas to exclude
    % this gives distance in MICRONS for every cell pixel to nearest edge
    pixSizMic=projData.pixSizeNm/1000; % side of a pixel in microns
    weightedRoi=bwdist(swapMaskValues(roiMask)).*pixSizMic;

    % inner mask will be used unless manual option is used
    if strcmpi('fraction',distUnits)
        distCutoff=distVal.*max(weightedRoi(:));
    else
        distCutoff=distVal;
    end
    % innerMask's boundary is fraction/microns inwards from the cellRoi boundary
    innerMask=weightedRoi>distCutoff;

    if innerMask==roiMask
        fractionFromEdge=1;
    else
        fractionFromEdge=0;
    end

    excludeMask=swapMaskValues(roiMask);

    % EXCLUDEMASK IS SUBMASK (OR INVERSE OF CELL MASK) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     pickExcludeInput=0;
    %     fileList=searchFiles('subMask',[],subanDir,0);
    %     if ~isempty(fileList)
    %         excludeMask=imread([fileList{1,2} filesep fileList{1,1}]);
    %     else
    %         excludeMask=~roiMask;
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pickExclude=pickExcludeInput;
    while pickExclude==1
        choice=questdlg({'Exclusion mask will be applied to all sub-rois of current project',...
            'created in this round. Press cancel to use full cell.'},...
            'Exclude region option','Draw new','Load a mask from file','Cancel','Draw new');
        if ~isempty(strmatch(choice,'Draw new'))
            figure
            imshow(roiMask.*img1,[]);
            hold on
            axis equal
            h=msgbox('Draw an ellipse over the area to exclude and double-click when finished','help');
            uiwait(h);
            h=imellipse;
            vert=wait(h);
            close(gcf)
            excludeMask=excludeMask | roipoly(img,vert(:,1),vert(:,2));
            close(gcf)

            reply = questdlg('Do you want to exclude another region?');
            if strcmpi(reply,'yes')
                pickExclude=1; % user said yes; make another one
            else
                pickExclude=0; % assume no; we're done
            end

        elseif ~isempty(strmatch(choice,'Load a mask from file'))
            c=1;
            while pickExclude==1
                try
                    [FileName,PathName] = uigetfile({'*.*'},'Select excludeMask.tif');
                    excludeMask=imread([PathName FileName]);
                    pickExclude=0;
                catch
                    h=msgbox('Please try again.','help');
                    uiwait(h);
                    c=c+1;
                    if c<=3 % give them 3 chances, then abort
                        return
                    end
                end
            end
        else % user hit cancel
            % don't use an exclusion mask
            pickExclude=0;
        end
    end

    % save whole cell roiMask and coordinates
    imwrite(roiMask,[subanDir filesep 'roiMask.tif']);
    save([subanDir filesep 'roiYX.mat'],'roiYX');

    switch selectType
        case 'manual'
            % this will be handled in the roiSelection part below

        case 'cellPeriphSingle'
            % split cell into two parts - central and periphery - based on some number
            % of microns or fraction from the periphery inwards
            roiSet=zeros(imL,imW,2);
            roiSet(:,:,1)=innerMask;
            roiSet(:,:,2)=roiMask-innerMask;

        case 'cellPeriphQuad'
            % divide the cell into 5 sub-rois consisting of a central polygon and four
            % quadrants around the periphery

            % get list of all the pixels on innerRoi's boundary
            if fractionFromEdge==1
                innerMaskYX=[nan nan];
            else
                [y1,x1]=ind2sub([imL,imW],find(innerMask,1));
                innerMaskYX = bwtraceboundary(innerMask,[y1,x1],'N');
            end

            figure
            imshow(roiMask.*img1,[])
            hold on
            axis equal
            plot(innerMaskYX(:,2),innerMaskYX(:,1))
            plot(roiYX(:,2),roiYX(:,1))
            h=msgbox('Draw a line across the cell and double-click when finished','help');
            uiwait(h);
            h=imline;
            position = wait(h);
            close(gcf)

            % position of the ends of the user-chosen line
            lineEndsYX=position(:,2:-1:1);

            % get roi centroid
            stats=regionprops(bwlabel(roiMask),'centroid');
            centerRoiYX=stats.Centroid(2:-1:1);

            [xAll,yAll]=meshgrid(1:imW,1:imL);
            if abs(lineEndsYX(1,1)-lineEndsYX(2,1))<1 % y's are the same, user chose horizontal line
                r11=(yAll<=lineEndsYX(1,1));
                r12=(yAll> lineEndsYX(1,1));
                r21=(xAll<=centerRoiYX(1,2));
                r22=(xAll> centerRoiYX(1,2));
            elseif abs(lineEndsYX(1,2)-lineEndsYX(2,2))<1 % x's are the same, user chose vertical line
                r11=(xAll<=lineEndsYX(1,2));
                r12=(xAll> lineEndsYX(1,2));
                r21=(yAll<=centerRoiYX(1,1));
                r22=(yAll> centerRoiYX(1,1));
            else
                % user-chosen line slope and y-intercept
                m1=diff(lineEndsYX(:,1))/diff(lineEndsYX(:,2));
                b1=lineEndsYX(1,1)-m1*lineEndsYX(1,2);
                % perpendicular line going through roi's centroid
                m2=-1/m1;
                b2=centerRoiYX(1,1)-m2*centerRoiYX(1,2);
                % y-coordinates of both lines across all the x-pixels
                yLine1=repmat(m1.*[1:imW]+b1,[imL,1]);
                yLine2=repmat(m2.*[1:imW]+b2,[imL,1]);
                % divide the image into two parts on either side of line 1
                r11=(yAll<=yLine1);
                r12=(yAll> yLine1);
                % divide the image into two parts on either side of line 2
                r21=(yAll<=yLine2);
                r22=(yAll> yLine2);
            end

            if fractionFromEdge==1
                roiSet=zeros(imL,imW,4);
                roiSet(:,:,1)=r11 & r21 & (roiMask-innerMask);
                roiSet(:,:,2)=r12 & r21 & (roiMask-innerMask);
                roiSet(:,:,3)=r12 & r22 & (roiMask-innerMask);
                roiSet(:,:,4)=r11 & r22 & (roiMask-innerMask);
            else
                roiSet=zeros(imL,imW,5);
                roiSet(:,:,1)=innerMask;
                roiSet(:,:,2)=r11 & r21 & (roiMask-innerMask);
                roiSet(:,:,3)=r12 & r21 & (roiMask-innerMask);
                roiSet(:,:,4)=r12 & r22 & (roiMask-innerMask);
                roiSet(:,:,5)=r11 & r22 & (roiMask-innerMask);
            end
        otherwise
            disp('Option not supported')
    end

    % iterate til the user is finished
    makeNewROI=1; % flag for making new roi
    while makeNewROI==1
        % make new roi_n image/analysis directories
        % string for number of files

        currentRoiAnDir=[subanDir filesep 'sub_' num2str(roiCount)];
        mkdir(currentRoiAnDir);
        mkdir([currentRoiAnDir filesep 'meta']);
        mkdir([currentRoiAnDir filesep 'feat']);

        %copyfile([anDir filesep 'feat' filesep 'movieInfo.mat'],[currentRoiAnDir filesep 'feat' filesep 'movieInfo.mat']);
        [a,b,c]=copyfile([anDir filesep 'feat' filesep 'movieInfo.mat'],[currentRoiAnDir filesep 'feat' filesep 'movieInfo1.mat']);
        [a,b,c]=movefile([currentRoiAnDir filesep 'feat' filesep 'movieInfo1.mat'],[currentRoiAnDir filesep 'feat' filesep 'movieInfo.mat']);

        switch selectType
            case 'manual'
                tempRoi=[];
                while isempty(tempRoi)
                    try
                        if makeNewROI==1
                            h=msgbox('Draw a sub-ROI','help');
                            uiwait(h);
                        end
                        figure
                        [tempRoi,polyXcoord,polyYcoord]=roipoly(img2show);
                    catch
                        disp('Please try again.')
                    end
                end
                close(gcf)
            case 'cellPeriphQuad'
                tempRoi=roiSet(:,:,roiCount-roiStart+1);
            case 'cellPeriphSingle'
                tempRoi=roiSet(:,:,roiCount-roiStart+1);
            otherwise
                disp('Option not supported')
        end

        % get intersection with max region in the cell
        tempRoi=tempRoi & roiMask & ~excludeMask;

        % add the current roi to the composite image
        [img2show]=addMaskInColor(img2show,tempRoi,cMap(max(1,mod(roiCount,10)),:));

        % shrink max region in the cell for next round by excluding current roi
        roiMask=roiMask-tempRoi;

        % get coordinates of vertices (ie all pixels of polygon boundary)
        if sum(tempRoi(:))~=0
            [y1,x1]=ind2sub([imL,imW],find(tempRoi,1)); % first pixel on boundary
            roiYX = bwtraceboundary(tempRoi,[y1,x1],'N'); % get all pixels on boundary
        else
            roiYX=[nan nan];
        end

        % save sub-roi mask and coordinates
        save([currentRoiAnDir filesep 'roiYX.mat'],'roiYX');
        imwrite(tempRoi,[currentRoiAnDir filesep 'roiMask.tif']);
        imwrite(excludeMask,[currentRoiAnDir filesep 'excludeMask.tif']);

        switch selectType
            case 'manual'
                reply = questdlg('Do you want to select another ROI?');
            case 'cellPeriphSingle'
                if roiCount<size(roiSet,3)+roiStart-1
                    reply='yes';
                else
                    reply='no';
                end
            case 'cellPeriphQuad'
                if roiCount<size(roiSet,3)+roiStart-1
                    reply='yes';
                else
                    reply='no';
                end
            otherwise
                disp('Option not supported')
        end
        if strcmpi(reply,'yes')
            makeNewROI=1; % user said yes; make another one
            roiCount=roiCount+1; % counter for current condition rois
        else
            makeNewROI=0; % assume no; we're done
        end

        % add new sub-roi to the list of projects for extraction
        n=length(projList);
        projList(n+1,1).imDir=imDir;
        projList(n+1,1).anDir=currentRoiAnDir;

    end % while makeNewROI==1


    % plot old sub-roi boundaries in white and new ones in color
    % add a number to center of each sub-roi to show sub-roi number
    for iRoi=1:roiCount
        if iRoi==1
            % get whole cell mask
            roiMask=imread([subanDir filesep 'roiMask.tif']);
            roiYX=load([subanDir filesep 'roiYX.mat'],'roiYX'); roiYX=roiYX.roiYX;

            % initialize mask to store accumulated subs for this
            % iteration
            subMask=zeros(size(roiMask));

            % plot cell outline on image
            h=figure;
            imshow(repmat(roiMask,[1 1 3]).*img)
            hold on
            plot(roiYX(:,2),roiYX(:,1),'w');
        end

        % load iRoi outline coordinates and mask
        currentRoiAnDir=[subanDir filesep 'sub_' num2str(iRoi)];
        tempRoiMask=imread([currentRoiAnDir filesep 'roiMask.tif']);
        roiYX=load([currentRoiAnDir filesep 'roiYX.mat']); roiYX=roiYX.roiYX;

        plot(roiYX(:,2),roiYX(:,1),'Color',cMap(max(1,mod(iRoi,10)),:));

        % make weighted mask using distance transform to find position
        % where text should go
        weightedRoi=bwdist(~tempRoiMask);
        [r,c]=find(weightedRoi==max(weightedRoi(:)));
        text(c(1),r(1),num2str(iRoi),'color','r','fontsize',14);

        % make mask containing all sub-rois selected for this project
        % during this run
        if iRoi>=roiStart
            subMask=subMask | tempRoiMask;
        end
    end


    imwrite(subMask,[subanDir filesep 'subMask_' num2str(roiStart) '_' num2str(roiCount) '.tif']);

    % save composite image
    indxStr1 = num2str(roiStart); indxStr2 = num2str(roiCount);
    saveas(h,[subanDir filesep 'subROIs_' indxStr1 '_' indxStr2 '.fig'])
    saveas(h,[subanDir filesep 'subROIs_' indxStr1 '_' indxStr2 fileExt])
    close(h)

    % create updated projList for the roi_x folder containing all the sub-projects
    cd('..')
    getProj(pwd);
end


    
    
% look for repeats and only extract from unique sub-rois
subDirList=projList2Cell(projList);
projCell=unique(subDirList(:,1));
nProj=length(projCell);
progressText(0,'Extracting tracks from Sub-ROIs');

for iProj=1:nProj
    % create new projData from original data and save it in new meta folder
    currentRoiAnDir=projCell{iProj,1};
    plusTipSubRoiExtractTracks(currentRoiAnDir,timeUnits,timeVal);
    progressText(iProj/nProj,'Extracting tracks from Sub-ROIs');
end

cd(homeDir)
warning(warningState);
disp('Sub-ROIs...finished')

end % if micropattern

function [img2show]=addMaskInColor(img,roiMask,c)
%subfunction to add new polygon outline to composite image - this is needed
%because you can't pass vector graphics info to roipoly function, and we
%want to be able to visualize the regions that have already been selected.
temp=double(bwmorph(roiMask,'remove'));
borderIdx=find(temp);
nPix=numel(roiMask);

img2show=img;
img2show(borderIdx)=c(1);
img2show(borderIdx+nPix)=c(2);
img2show(borderIdx+2*nPix)=c(3);
