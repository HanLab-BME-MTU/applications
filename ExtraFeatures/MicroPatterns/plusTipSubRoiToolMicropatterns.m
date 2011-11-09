function plusTipSubRoiToolMicropatterns(projList,selectType,distUnits,distVal,timeUnits,timeVal,cellRoiYX,pickExclude)
% plusTipSubRoiTool allows sub-ROI selection and extracts MT growth tracks
%
% SYNOPSIS: plusTipSubRoiTool(projList,selectType,distUnits,distVal,timeUnits,timeVal,cellRoiYX,pickExclude)
%
% Crazy Dirty Code for Generating SubRois, Segregating Stats, and Performing Analysis for Mijungs Micropatterns 
% MB 04-2011 (NOTE super inefficient but gets the job done) 
%% PARAMETERS
%
% Copyright (C) 2011 LCCB 
%
% This file is part of plusTipTracker.
% 
% plusTipTracker is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% plusTipTracker is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with plusTipTracker.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Choose H or Y Pattern 
% HPattern = 0 is YPattern (ie Default is YPattern)
% HPattern = 1 is HPattern (ie turn on the HPattern)

% Choice of H or Y Pattern
    
HPattern = 1;
        
    %%% CELL REGION PARAMETERS %%%
    numWindows = 2; % Number Windows Moving Toward Cell Center total number
    % of windows will equal n numWindows with defined region from cell edge
    % plus a leftover region in the interior (if numWindows = 2, total
    % windows will equal 3)
    windowSize = distVal; % Width of Outer Windows in Microns (Relative to Cell Edge) IN MICRONS
    %%%%%%%%%%%%%%%%%
 
% Choose Name for folders: input as strings 

subRoiFolderName = 'SUBROIS_3um_TargetedToRegion_MB_11_08';
analysisFolderName = 'ANALYSIS_3um_TargetedToRegion_MB_11_08';

%%% Choose if want to also perform analysis 
doAnalysis = 1; % set this to 1 to "turn-on" anlaysis 

%Parameters for Pooled Analysis 
doBtw = 0; % set to 1 to perform between Group Comparisons 
% Here we are comparing different subRegions within the same condition.
doWtn = 0; % set to 1 to perform within Group Comparisons: set to 1 if you want to take a look 
% at the variability between the individual subregions before pooling
doPlot = 1; % set to 1 to make boxplots
removeBeginEnd = 1; % 1 if you want to remove tracks from beginning and end for 
                    % analyis -- note need to check for subRois that this
                    % is not buggy. 0 to analyze all tracks 
justExtractTracks = 0;% ignore this for now
%% Initialize: CheckInput Parameters
homeDir=pwd;
warningState = warning;
warning('off','Images:initSize:adjustingMag')
warning('off','MATLAB:divideByZero')

fileExt='.tif';

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

nProj=length(projList);
originalList = projList;

if justExtractTracks ~= 1 
%%
for iProj=1:nProj

    anDir=projList(iProj,1).anDir;

    % check if current project is a sub-roi.  if it is, we extract tracks
    % only, since selecting sub-rois recursively is not allowed.  if it is
    % a full roi, we go on with roi selection.
    if ~isempty(strfind(anDir,'sub'))
        continue
    end
  
    subanDir=[anDir filesep subRoiFolderName];
   
    if isdir(subanDir) ==1 
        rmdir(subanDir,'s')
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

    cMap=hsv(9);



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
           % choice=questdlg('Before creating Sub-ROIs, you need to define the cell boundary.','Cell ROI option','Draw new','Load roiYX.mat','Draw new');
            %if ~isempty(strmatch(choice,'Draw new'))
             %   c=1;
              %  while isempty(roiMask)
               %     try
                %       figure
                 %       [roiMask,polyXcoord,polyYcoord]=roipoly(img);
                  %      close
                   %     roiYX=[polyYcoord polyXcoord; polyYcoord(1) polyXcoord(1)];
                   %catch
                    %    h=msgbox('Please try again.','help');
                     %   uiwait(h);
                      %  c=c+1;
                       % if c<=3 % give them 3 chances, then abort
                        %    return
                        %end
                   % end
                %end
                %else % load from file
                %c=1;
                 %while isempty(roiMask)
                  %   try
                   %      [FileName,PathName] = uigetfile({'*.*'},'Select roiYX.mat');
                         %p =load([PathName FileName]);
                 p = load([anDir filesep 'roiYXSeg.mat']); % load the roiYX file in the anDir
                        roiYX=p.roiYX;
                        roiMask=roipoly(img,roiYX(:,2),roiYX(:,1));
                   % roiMask = imread([anDir filesep 'roiMaskSeg.tif']); 
                  % catch
                    %    h=msgbox('Please try again.','help');
                      %  uiwait(h);
                     %   c=c+1;
                       % if c<=3 % give them 3 chances, then abort
                        %   return
                        %end
                  %  end
                %end
            %end
        end
    end


    % set cell boundary in white to composite image
    [img2show]=addMaskInColor(img2show,roiMask,[1 1 1]);

    % make distance transform before getting areas to exclude
    % this gives distance in MICRONS for every cell pixel to nearest edge
    pixSizMic=projData.pixSizeNm/1000; % size of a pixel in microns
    weightedRoi=bwdist(swapMaskValues(roiMask)).*pixSizMic; %Distance Transform

    % inner mask will be used unless manual option is used
    if ~isempty(strmatch('fraction',distUnits))
        distCutoff=distVal.*max(weightedRoi(:));
    else
        distCutoff=distVal;
    end
%%  Code Added to Calculate Multiple Inner Masks
%   innerMask's boundary is fraction/microns inwards from the cellRoi boundary
    
    
    if windowSize ~=0 
        
        % initialize matrix to hold inner masks 
        innerMasks = zeros(imL,imW,numWindows);
    
            % fill matrix holding distance transform derived masks
            for iWindow = 1:numWindows
                innerMasks(:,:,iWindow) = weightedRoi>windowSize*iWindow;
            end 
        else % no need to calculate inner masks as user does not want to do windowing
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
            
            roiSet(:,:,1)=roiMask-innerMaskLarge;
            roiSet(:,:,2)=roiMask-innerMaskSmall;
            roiSet(:,:,3)=innerMaskSmall;
        case 'cellPeriphQuad'
            
            % divide the cell into 5 sub-rois consisting of a central polygon and four
            % quadrants around the periphery
            
            % get list of all the pixels on innerRoi's boundary
            %if fractionFromEdge==1
             %   innerMaskYX=[nan nan];
            %else
            %    [y1,x1]=ind2sub([imL,imW],find(innerMask,1));
            %   % 
%if justExtractSubTracks == 1
   %  innerMaskYX = bwtraceboundary(innerMask,[y1,x1],'N');
            %end

            %figure
            %imshow(roiMask.*img1,[])
            %hold on
            %axis equal
            %plot(innerMaskYX(:,2),innerMaskYX(:,1))
            %plot(roiYX(:,2),roiYX(:,1))
            %h=msgbox('Draw a line across the cell and double-click when finished','help');
            %uiwait(h);
            %h=imline;
            %position = wait(h);
            %close(gcf)

            % position of the ends of the user-chosen line
            %Note:(:,2:-1:1) just a way of inverting x and y coordinates
            %lineEndsYX=position(:,2:-1:1);
            
            
%% Start Modifications for Micropatterns: Parameters
    
    
    if HPattern== 1 
        angle1 = 60;
        angle2 = 30;
    else 
        
     %   angle1 = 45;
      angle2 = 15;
      %  angle3 = 75;
     
    end 
        
    if iProj > 1 
         save([subanDir filesep 'sub_ROIs_all'], 'roiSet');

        clear roiSet
    else 
    end 
%% Modifications for Micropatterns: Body (current set up H-Pattern)
            % get roi centroid
            stats=regionprops(bwlabel(roiMask),'centroid');
            centerRoiYX=stats(1).Centroid(2:-1:1);
             %centerRoiYX=stats.Centroid;
             
            [xAll,yAll]=meshgrid(1:imW,1:imL);
           
                %Divide Image Along X and Y axes
                r11=(xAll<=centerRoiYX(1,2));
                r12=(xAll> centerRoiYX(1,2));
                r21=(yAll<=centerRoiYX(1,1));
                r22=(yAll> centerRoiYX(1,1));
                
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if HPattern== 1 % H-Pattern
               
                    m1 = tand(angle1); %slope of line (NOTE: will look like has a neg   
                    m2 =-m1; % negative slope of line1 
                    
                
                % Find the line that is angle1 from the x-axis and 
                % passes throught the image center
              
                % slope in image b/c y-axis reveresed for images)
                b1 = centerRoiYX(1,1)-m1*centerRoiYX(1,2); % y-intercept
                yLine1 = repmat(m1.*[1:imW]+b1,[imL,1]); % Line coordinates 
                
                %Divide Image Based on This Line
                d11= (yAll<=yLine1); % 1st angle negative slope: below (again reversal based on image coordinate system)
                d12 = (yAll>yLine1); % 1st angle negative slope: above
                
                % Find Corresponding Line -Angle1 From 
               
                b2 = centerRoiYX(1,1)-m2*centerRoiYX(1,2); % find y-intercept
                yLine2 = repmat(m2.*[1:imW]+b2,[imL,1]); % Line2 Coordinates
                d21 = (yAll<=yLine2); % 1st angle positive slope below 
                d22 = (yAll>yLine2); % 1st angle positive slope above  
                
                
                else % y-pattern
                    
                
                    sortedMaskX = sortRows(roiYX,2);
                    minXCoord = sortedMaskX(1,:);
                    maxXCoord = sortedMaskX(end,:); 
               
                    angleWidth = 15;
                    
                    % find slope 
                    xMinMidPt = (minXCoord(1,1) - centerRoiYX(1,1))/(minXCoord(1,2)-centerRoiYX(1,2));
                    angleXMinMidPt = atand(xMinMidPt);
                    angleAd1 = angleXMinMidPt - angleWidth;
                    angleAd2 = angleXMinMidPt + angleWidth;
                    
                    xMaxMidPt = (maxXCoord(1,1) - centerRoiYX(1,1))/(maxXCoord(1,2)-centerRoiYX(1,2));
                    angleXMaxMidPt = atand(xMaxMidPt);
                    angleAd3 = angleXMaxMidPt - angleWidth;
                    angleAd4 = angleXMaxMidPt + angleWidth;
                    
                    m1 = tand(angleAd1); 
                    m2 = tand(angleAd2);
                    m3 = tand(angleAd3);
                    m4 = tand(angleAd4);
                    
                    
                % Find the line that is angle1 from the x-axis and 
                % passes throught the image center
              
                % slope in image b/c y-axis reveresed for images)
                b1 = centerRoiYX(1,1)-m1*centerRoiYX(1,2); % y-intercept
                yLine1 = repmat(m1.*[1:imW]+b1,[imL,1]); % Line coordinates 
                
                %Divide Image Based on This Line
                d11= (yAll<=yLine1); % 1st angle negative slope: below (again reversal based on image coordinate system)
                d12 = (yAll>yLine1); % 1st angle negative slope: above
                
                % Find Corresponding Line -Angle1 From 
               
                b2 = centerRoiYX(1,1)-m2*centerRoiYX(1,2); % find y-intercept
                yLine2 = repmat(m2.*[1:imW]+b2,[imL,1]); % Line2 Coordinates
                d21 = (yAll<=yLine2); % 1st angle positive slope below 
                d22 = (yAll>yLine2); % 1st angle positive slope above
                
                
                b3 = centerRoiYX(1,1)-m3*centerRoiYX(1,2);
                yLine3 = repmat(m3.*[1:imW]+b3,[imL,1]);
                d31 = (yAll<=yLine3);
                d32 = (yAll>yLine3);
                
                b4 = centerRoiYX(1,1)-m4*centerRoiYX(1,2);
                yLine4 = repmat(m4.*[1:imW]+b4,[imL,1]);
                d41 = (yAll<=yLine4);
                d42 = (yAll>yLine4);
                
                
                %Make Top Adhesion Regions
                
                AdTopLeft = d22 & d11 & roiMask;
                
                AdTopRight = d41 & d32 & roiMask;
                    
                NonAdTop =  d21 & d31 & roiMask;
                end   % if HPattern     
                
                % Window Ad, NonAd if user desires
                if exist('innerMasks','var') == 1
                
                    % Initialize matrix to hold the distance cut-off rois (For
                    % all) 
                    cellEdgeRoi= zeros(imL,imW,numWindows+1);
                
                    % Make Cell Edge Specific Regions Based on specified distance from the cell edge 
                    cellEdgeRoi(:,:,1) = roiMask - innerMasks(:,:,1);
                    for i = 1:(size(innerMasks,3)-1)
                        cellEdgeRoi(:,:,i+1) = innerMasks(:,:,i)-innerMasks(:,:,i+1);
                    end 
                    cellEdgeRoi(:,:,numWindows+1) = innerMasks(:,:,numWindows);
                
                end 
                
           
              
                if HPattern == 1
                 
                roiNonAdTop = r21 & d11 & d21 & roiMask;
                
                % if the user wants to do windowing need to further divide
                % region
                 if exist('cellEdgeRoi','var') == 1
                    for i= 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i) = roiNonAdTop & cellEdgeRoi(:,:,i);
                    end 
                 else
                     roiSet(:,:,1) = roiNonAdTop; % no further subdivision of regions
                 end 
                 
               
                roiNonAdBot = r22 & d22 & d12 & roiMask;
                
                
               
                 % Divide Non-Ad region inot Cell-Edge Based Windows
                if exist('cellEdgeRoi','var') == 1
               
                
                %Find the number of those already set
                numSet = size(roiSet,3);
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = roiNonAdBot & cellEdgeRoi(:,:,i);
                    end 
                else 
                    roiSet(:,:,2) = roiNonAdBot; % no further subdivision of regions
                    
                end
                
                
                %Large regions new name
                roiAdLeft1Big = r11 & r21 & d12 & roiMask;
                roiAdLeft2Big = r11 & r22 & d21 & roiMask;
                
               
                
                roiAdRight1Big = r12 & r21 & d22 & roiMask; 
                roiAdRight2Big = r12 & r22 & d11 & roiMask;
               
               
                % if angle2 is not empty find the lines with an absolute angle2  
                % from the x-axis and intersects the cell center
                %. Use these lines to subdivide the image further
               
                  
                   
                       
                end 
                
                if ~isempty(angle2)
               
                    if HPattern== 1
               
                % Find Line3 (angle2 from x-axis) 
                m3 = tand(angle2); % find slope 
                m4 =-m3;
               
                
                b3 = centerRoiYX(1,1)-m3*centerRoiYX(1,2); % find y-intercept
                yLine3 = repmat(m3.*[1:imW]+b3,[imL,1]); % line coordinates
                
                %Subdivide Image Using Line3 
                d31 = (yAll<=yLine3);
                d32 = (yAll>yLine3);
                
                % Find Line4 (-angle2 from x-axis)
             
                b4 = centerRoiYX(1,1)-m4*centerRoiYX(1,2);
                yLine4 = repmat(m4.*[1:imW]+b4,[imL,1]);
                d41 = (yAll<=yLine4);
                d42 = (yAll>yLine4);
                
                    
                
                
                 % cluster subregions based on H-pattern
                    
                % Adhesion regions to combine
                AdTopLeft = roiAdLeft1Big & d32;
                AdTopRight =  roiAdRight1Big & d42  ;
                AdBotLeft = roiAdLeft2Big &  d41  ;
                AdBotRight = roiAdRight2Big &  d31;
               

                AdLeft = AdBotLeft + AdTopLeft;
                AdRight = AdBotRight + AdTopRight;
                
                if exist('cellEdgeRoi','var') == 1 
                    numSet = size(roiSet,3);
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = AdLeft & cellEdgeRoi(:,:,i);
                    end 
                
               
                
                    numSet = size(roiSet,3);
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = AdRight & cellEdgeRoi(:,:,i);
                    end 
                else 
                    roiSet(:,:,3) = AdLeft;
                    roiSet(:,:,4) = AdRight;
                    
                end  
             
                
                    %hints r11 = left hand y-axis r12 = right hand of
                    %y=axis r22 below x-axis
               
                    
                % 1 above 2 below
                AdCorners = zeros(imL,imW,4);
                AdCorners(:,:,1) = r21 & d12 & d31 & roiMask; % top left 
                AdCorners(:,:,2) = r21 & d22 & d41 & roiMask; % top right
                AdCorners(:,:,3) = r22 & d21 & d42 & roiMask; % bottom 
                AdCorners(:,:,4) = r22 & d11 & d32 &  roiMask; % bottom
                
                if exist('cellEdgeRoi','var') == 1 
                numSet = size(roiSet,3);
              
                    for iCorn = 1:4 
                        for i = 1:size(cellEdgeRoi,3)
                            roiSet(:,:,i+numSet) = AdCorners(:,:,iCorn) & cellEdgeRoi(:,:,i);
                         
                        end 
                        numSet = size(roiSet,3);
                    end 
                else 
                    
                    roiSet(:,:,5) = AdCorners(:,:,1);
                    roiSet(:,:,6) = AdCorners(:,:,2);
                    roiSet(:,:,7) = AdCorners(:,:,3);
                    roiSet(:,:,8) = AdCorners(:,:,4); 
                    
                end   
                  
              
               
                
                else % y-pattern 
                    
                  
                    sortedMask = sortRows(roiYX);
                    minYCoord = sortedMask(end,:); 
                    
                    angleWidth = 15;
                    % find slope 
                    mMidPt = (minYCoord(1,1) - centerRoiYX(1,1))/(minYCoord(1,2)-centerRoiYX(1,2));
                    angleMidPt = atand(mMidPt);
                    angleAdLt = angleMidPt - angleWidth;
                    angleAdRt = angleMidPt + angleWidth;
                    
                    m5 = tand(angleAdLt); 
                    m6 = tand(angleAdRt);
                    
                    
                     b5 = centerRoiYX(1,1)-m5*centerRoiYX(1,2); % y-intercept
                     yLine5 = repmat(m5.*[1:imW]+b5,[imL,1]); % Line coordinates 
                
                    %Divide Image Based on This Line
                    d51= (yAll<=yLine5); % 1st angle negative slope: below (again reversal based on image coordinate system)
                    d52 = (yAll>yLine5); % 1st angle negative slope: above
                
                    % Find Corresponding Line -Angle1 From 
               
                    b6 = centerRoiYX(1,1)-m6*centerRoiYX(1,2); % find y-intercept
                    yLine6 = repmat(m6.*[1:imW]+b6,[imL,1]); % Line2 Coordinates
                    d61 = (yAll<=yLine6); % 1st angle positive slope below 
                    d62 = (yAll>yLine6); % 1st angle positive slope above
                
                    
                    AdBot = d52 & d62 & roiMask;
                    
                    NonAdBotLt = d61 & d12 & roiMask;
                   
                    NonAdBotRt = d51 & d42 & roiMask;
                    
                    
                    
                    
                   if exist('cellEdgeRoi','var') == 1
                   %Divide Ad Left Region Into Cell-Edge Based Windows 
                     for i= 1:size(cellEdgeRoi,3)
                     roiSet(:,:,i) = AdTopLeft & cellEdgeRoi(:,:,i);
                     end 
                    
                    
                    
      
                    %Find the number of those already set
                    numSet = size(roiSet,3);
                
                    % Divide Ad Right Region into Cell-Edge Based Windows
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = AdTopRight & cellEdgeRoi(:,:,i);
                    end 
                   
                   %Find the number of those already set
                    numSet = size(roiSet,3);
                
                    % Divide Non-Ad region inot Cell-Edge Based Windows
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = AdBot & cellEdgeRoi(:,:,i);
                    end
                    
                  
                   %Find the number of those already set
                    numSet = size(roiSet,3);
                
                    % Divide Non-Ad region inot Cell-Edge Based Windows
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = NonAdBotLt & cellEdgeRoi(:,:,i);
                    end
                   
                   
                   
                 
                   %Find the number of those already set
                    numSet = size(roiSet,3);
                
                    % Divide Non-Ad region inot Cell-Edge Based Windows
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = NonAdTop & cellEdgeRoi(:,:,i);
                    end
                   
                   
                  
                   
                   
                   
                   %Find the number of those already set
                    numSet = size(roiSet,3);
                
                    % Divide Non-Ad region inot Cell-Edge Based Windows
                    for i = 1:size(cellEdgeRoi,3)
                        roiSet(:,:,i+numSet) = NonAdBotRt & cellEdgeRoi(:,:,i);
                    end
                   
                   end  % if exist 
                   
                end 
                    
                    
                end 
              
               
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
        dir{roiCount} = currentRoiAnDir; 
        

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
    subDirList{iProj} = dir';

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
    %cd('..')
    %getProj(pwd);
    
    % save([subanDir filesep 'projListSubRois'],'projList');
end % for iProj

else % skip above if just extract tracks is turned on
    
  
   
   
end % if justExtractTracks
    

%% Partition The Data for Each Sub-Roi

projCell = vertcat(subDirList{:}); 


% look for repeats and only extract from unique sub-rois
%subDirList=projList2Cell(projList);
%projCell=unique(subDirList(:,1));

nProj=length(projCell);
progressText(0,'Extracting tracks from Sub-ROIs');
for iProj=1:nProj
    % create new projData from original data and save it in new meta folder
    currentRoiAnDir=projCell{iProj,1};
    plusTipSubRoiExtractTracks(currentRoiAnDir,timeUnits,timeVal);
    progressText(iProj/nProj,'Extracting tracks from Sub-ROIs');
end

%% Start Make GroupLists


    pathUp1 = getFilenameBody(projList(1).imDir);
    [pathUp2,groupName] = getFilenameBody(pathUp1);    
    
    statDir = [pathUp2 filesep analysisFolderName];
    
    if isdir(statDir)
        rmdir(statDir,'s')
    else 
    end 
      mkdir(statDir);
   save([statDir filesep 'projCellSubRois'],'projCell'); 
 
    %Make GroupLists for Each Type of Window   
    %for all types of windows make a groupList combining all projects
    
    nProjOrig = length(originalList);
   
    groupListDir = [statDir filesep 'groupLists'];
    mkdir(groupListDir);

if windowSize ~= 0 

if HPattern == 1    
   
        
    %%%%% GroupLists Non-Adhesion %%%%%%%%%%
    
    %Individual groupLists for Small Windows From Cell Edge
    for iWindow = 1:numWindows
    
     %NonAdDir = [statDir filesep 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
     %mkdir(NonAdDir);
   
    
    groupListCount = 1;
    for iProj = 1:nProjOrig
        
        groupList{groupListCount,1} = [ groupName 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(iWindow)];
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(iWindow+numWindows+1)];
    
        groupListCount = groupListCount + 2;
    end % if Proj
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListNonAdCompare = groupList;
    else 
        groupListNonAdCompare = [groupListNonAdCompare; groupList];
    end
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListNonAd' num2str(iWindow*windowSize) '_uM'],'groupList');
    
    % Hold Window Groups for Between Region Comparisons
    groupListNonAd(:,:,iWindow) = groupList;
    
    end  % end iWindow   
    
    % Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:nProjOrig
            groupList{groupListCount,1} = [ groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(numWindows+1)];
            groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(2*(numWindows+1))];
            groupListCount = groupListCount + 2;
        end
        
     %NonAdDir = [statDir filesep 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
     %mkdir(NonAdDir);
     
     %Save Individual GroupList For Central Region
     save([groupListDir filesep 'groupListNonAd_GreaterThan' num2str(numWindows*windowSize) 'uM'],'groupList');   
     
     % Hold Central Region Group for Between Region Comparisons
     groupListNonAd(:,:,numWindows+1) = groupList;
   
     %Add Central Region to Window Comparison GroupList 
     groupListNonAdCompare = [groupListNonAdCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListNonAdCompare; 
     save([groupListDir filesep 'groupListCompareWindows_NonAd'],'groupList');
     
     %%%%%% GroupLists Adhesion %%%%%% 
     
       
        clear groupList
       
        %Individual groupLists for Small Windows From Cell Edge
        for iWindow = 1:numWindows
    
        %AdDir = [statDir filesep 'Adhesion_' num2str(iWindow*windowSize) 'uM'];
        %mkdir(AdDir);
   
    
        groupListCount = 1;
        for iProj = 1:nProjOrig
            groupList{groupListCount,1} = [ groupName 'Adhesion_' num2str(iWindow*windowSize) 'uM']; 
            groupList{groupListCount+1,1} = [ groupName 'Adhesion_' num2str(iWindow*windowSize) 'uM'];
            windowNumRight = 2*(numWindows+1) + iWindow;
            windowNumLeft = windowNumRight+numWindows + 1;
            groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNumRight)];
            groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNumLeft)];
             
            groupListCount = groupListCount +2;
        end  
        
         %Create GroupList to Compare Among Windows in Same Region Type
        if iWindow == 1
            groupListAdCompare = groupList;
        else 
            groupListAdCompare = [groupListAdCompare; groupList];
        end
        
         % Hold Window Groups for Between Region Comparisons
            groupListAd(:,:,iWindow) = groupList;
        
        %Save GroupList Files for Each Window Size/Region Type 
        save([groupListDir filesep 'groupListAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
        
        
       
        
        end %for iWindows
        
       
        
        % Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:nProjOrig
            
            groupList{groupListCount,1} = [groupName 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            windowNumRight = 2*(numWindows+1) + numWindows+1;
            windowNumLeft = windowNumRight+numWindows+1;
            groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNumRight)]; % 
            groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNumLeft)];
            groupListCount = groupListCount + 2;
        end
        
     %AdDir = [statDir filesep 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
     %mkdir(AdDir);
     
     
     %Save Individual GroupList For Central Region
     save([groupListDir filesep 'groupListAd_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
     
     % Hold Central Region Group for Between Region Comparisons
     groupListAd(:,:,numWindows+1) = groupList;
     
     %Add Central Region to Window Comparison GroupList 
     groupListAdCompare = [groupListAdCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListAdCompare;
     save([groupListDir filesep 'groupListCompareWindows_Ad'],'groupList');
     
 %%    
     
      %%%%%% GroupLists Adhesion Corners %%%%%% 
      clear groupList
       

    for iWindow = 1:numWindows
    
    % AdCornDir = [statDir filesep 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
    % mkdir(AdCornDir);
   
   
    groupListCount = 1;
    
    for iProj = 1:nProjOrig
        
        groupList{groupListCount,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+2,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+3,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        
        windowNum1 = 4*(numWindows+1) + iWindow;
        windowNum2 = windowNum1 + numWindows +1;
        windowNum3 = windowNum2 +numWindows +1;
        windowNum4 = windowNum3 + numWindows +1;
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)];
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
        groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
        groupList{groupListCount+3,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum4)];
    
        groupListCount = groupListCount + 4;
           
    end % for iProj
    
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListAdCornCompare = groupList;
    else 
        groupListAdCornCompare = [groupListAdCornCompare; groupList];
    end % end if
     
    % Hold Window Groups for Between Region Comparisons
            groupListAdCorn(:,:,iWindow) = groupList;
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListAdCorn_' num2str(iWindow*windowSize) '_uM'],'groupList');
    
    end   % for iWindows 
 
    
     %Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:nProjOrig
            
            groupList{groupListCount,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+2,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+3,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            
            windowNum1 = 4*(numWindows+1) + numWindows +1;
            windowNum2 = windowNum1 + numWindows + 1;
            windowNum3 = windowNum2 + numWindows +1;
            windowNum4 = windowNum3 + numWindows +1;
            
            groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)]; 
            groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
            groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
            groupList{groupListCount+3,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum4)];
            
            groupListCount = groupListCount + 4;
        end % for iProj
        
    % AdCornDir = [statDir filesep 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
    % mkdir(AdCornDir);
    
    %Save Individual Group List for Central Region 
    save([groupListDir filesep 'groupListAdCornGreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
    
    % Hold Central Region Group for Between Region Comparisons
    groupListAdCorn(:,:,numWindows+1) = groupList;
    
    %Add Central Region to Window Comparison GroupList 
     groupListAdCornCompare = [groupListAdCornCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListAdCornCompare;
     save([groupListDir filesep 'groupListCompareWindows_AdCorn'],'groupList');
     
     
    
    
     clear groupList
     
 %%%%%%%%%%%%%%%%%% GroupList Same Window Different Region Type %%%%%%%
    for iWindow = 1:numWindows
        groupList = [groupListAd(:,:,iWindow) ; groupListAdCorn(:,:,iWindow); groupListNonAd(:,:,iWindow)];
        save([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM'],'groupList');
    end % for iWindow
    
    groupList = [groupListAd(:,:,numWindows+1) ; groupListAdCorn(:,:,numWindows+1) ; groupListNonAd(:,:,numWindows+1)];
    save([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM'],'groupList');
    
   %[groupData]  = plusTipPoolGroupData(groupList,1,0,1,1);
   
  %[discrimMat] =  plusTipTestDistr(groupData,[]); 
   %[ ] = ('stats',groupData,[],1,1,1,0);
%% Y-Pattern 
else 
    %%%%%% Create Y-Pattern GroupLists NonAdhesion Regions %%%%%% 
      clear groupList
    for iWindow = 1:numWindows
    
    % AdCornDir = [statDir filesep 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
    % mkdir(AdCornDir);
   
   
    groupListCount = 1;
    
    for iProj = 1:nProjOrig
        
        groupList{groupListCount,1} = [ groupName '_NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName '_NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+2,1} = [ groupName '_NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
      
        windowNum1 = 3*(numWindows + 1) + iWindow;
        windowNum2 = windowNum1 + numWindows +1;
        windowNum3 = windowNum2 +numWindows +1;
      
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)];
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
        groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
       
    
        groupListCount = groupListCount + 3;
           
    end % for iProj
    
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListNonAdhesionCompare = groupList;
    else 
        groupListNonAdhesionCompare = [groupListNonAdhesionCompare; groupList];
    end % end if
     
    % Hold Window Groups for Between Region Comparisons
            groupListNonAd(:,:,iWindow) = groupList;
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListNonAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
    
    end   % for iWindows 
 
    
     %Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:nProjOrig
            
            groupList{groupListCount,1} = ['grp_' groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = ['grp_' groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+2,1} = ['grp_' groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
      
            
            windowNum1 = 3*(numWindows +1) + numWindows +1 ;
            windowNum2 = windowNum1 + numWindows + 1;
            windowNum3 = windowNum2 + numWindows +1;
    
            
            groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)]; 
            groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
            groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];  
            
            groupListCount = groupListCount + 3;
        end % for iProj
        
    % AdCornDir = [statDir filesep 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
    % mkdir(AdCornDir);
    
    %Save Individual Group List for Central Region 
    save([groupListDir filesep 'groupListNonAdGreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
    
    % Hold Central Region Group for Between Region Comparisons
    groupListNonAd(:,:,numWindows+1) = groupList;
    
    %Add Central Region to Window Comparison GroupList 
     groupListNonAdhesionCompare = [groupListNonAdhesionCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListNonAdhesionCompare;
     save([groupListDir filesep 'groupListCompareWindows_NonAd'],'groupList');
    
%% 
 %%%%%% Create Y-Pattern GroupLists: Adhesion Regions %%%%%% 
      clear groupList
       

    for iWindow = 1:numWindows
    
   
   
   
    groupListCount = 1;
    
    for iProj = 1:nProjOrig
        
        groupList{groupListCount,1} = [ groupName '_Adhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName '_Adhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+2,1} = [ groupName '_Adhesion_' num2str(iWindow*windowSize) 'uM'];
      
        windowNum1 = iWindow;
        windowNum2 = windowNum1 + numWindows +1;
        windowNum3 = windowNum2 +numWindows +1;
      
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)];
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
        groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
       
    
        groupListCount = groupListCount + 3;
           
    end % for iProj
    
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListAdhesionCompare = groupList;
    else 
        groupListAdhesionCompare = [groupListAdhesionCompare; groupList];
    end % end if
     
    % Hold Window Groups for Between Region Comparisons
            groupListAd(:,:,iWindow) = groupList;
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
    
    end   % for iWindows 
 
    
     %Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:nProjOrig
            
            groupList{groupListCount,1} = [ groupName '_Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName '_Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+2,1} = [ groupName '_Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
      
            
            windowNum1 = numWindows +1;
            windowNum2 = windowNum1 + numWindows + 1;
            windowNum3 = windowNum2 + numWindows +1;
    
            
            groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)]; 
            groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
            groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];  
            
            groupListCount = groupListCount + 3;
        end 
    
  %Save Individual Group List for Central Region 
    save([groupListDir filesep 'groupListAdGreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
    
    % Hold Central Region Group for Between Region Comparisons
    groupListAd(:,:,numWindows+1) = groupList;
    
    %Add Central Region to Window Comparison GroupList 
     groupListAdhesionCompare = [groupListNonAdhesionCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListNonAdhesionCompare;
     save([groupListDir filesep 'groupListCompareWindows_Ad'],'groupList');
    
   %%%%%%%%%%%%%%%%%% GroupList Same Window Different Region Type %%%%%%%
    for iWindow = 1:numWindows
        groupList = [groupListAd(:,:,iWindow) ; groupListNonAd(:,:,iWindow)];
        save([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM'],'groupList');
    end % for iWindow
    
    groupList = [groupListAd(:,:,numWindows+1) ; groupListNonAd(:,:,numWindows+1)];
    save([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM'],'groupList'); 
     
end % End if H-Pattern   
end 
clear groupList;
if windowSize == 0 % window size = 0 so don't need all the above
    groupListCount = 1;
    for iProj = 1:nProjOrig
       
        groupList{groupListCount,1} = [groupName 'NonAdhesion'];
        groupList{groupListCount+1,1} = [groupName 'NonAdhesion'];
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_1'];
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_2'];
        
        groupListCount = groupListCount + 2;
       
    end     
    
   
    save([groupListDir filesep 'groupListNonAd'],'groupList');
    groupListNonAd = groupList; 
    
    
    clear groupList 
    
    groupListCount = 1;
    for iProj = 1:nProjOrig
        groupList{groupListCount,1} = [groupName 'Adhesion'];
        groupList{groupListCount+1,1} = [groupName 'Adhesion'];
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_3']; 
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_4'];
    
        groupListCount = groupListCount + 2;
    end 
    
    
    save([groupListDir filesep 'groupListAd'],'groupList');
    groupListAd = groupList; 
    
    
    clear groupList
    groupListCount = 1;
    for iProj = 1:nProjOrig 
        groupList{groupListCount,1} = [groupName 'AdhesionCorner']; 
        groupList{groupListCount+1,1} = [groupName 'AdhesionCorner'];
        groupList{groupListCount+2,1} = [groupName,'AdhesionCorner'];
        groupList{groupListCount+3,1} = [groupName,'AdhesionCorner']; 
        
        groupList{groupListCount,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_5'];
        groupList{groupListCount+1,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_6'];
        groupList{groupListCount+2,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_7'];
        groupList{groupListCount+3,2} = [projList(iProj).anDir filesep subRoiFolderName filesep 'sub_8'];
        
        groupListCount = groupListCount + 4;
    end
    
    save([groupListDir filesep 'groupListAdCorn'],'groupList');
    groupListAdCorn= groupList;
    clear groupList 
    groupList = [groupListNonAd; groupListAd; groupListAdCorn];
    
    save([groupListDir filesep 'groupListCompareRegions'],'groupList');
end  

%% Perform Analysis



if doAnalysis == 1

% Make Directories for Analysis

analSaveDir1 = [ statDir filesep 'CompareBtwRegionTypes'];
mkdir(analSaveDir1);

analSaveDir2 = [statDir filesep 'CompareWindowsWithinRegionType'];
mkdir(analSaveDir2);



if HPattern== 1
    regionTypes{1,1} = 'Ad';
    regionTypes{2,1} = 'AdCorn';
    regionTypes{3,1} = 'NonAd';
    
else 
    regionTypes{1,1} = 'Ad';
    regionTypes{2,1} = 'NonAd';
end

% Pool the data from the groupList

for iWindow = 1:numWindows

    saveDir = [analSaveDir1 filesep num2str(iWindow*windowSize) '_uM'];
    mkdir(saveDir);
    
    
    
    %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM.mat']);
    
    groupList = groupList.groupList;
    %Perform Pooling of Data
    groupData = plusTipExtractGroupData(groupList,removeBeginEnd); 
    save([saveDir filesep 'groupData'], 'groupData'); 
    plusTipPoolGroupData(groupData, saveDir, doWtn, doPlot); 
    
    % Collect Data
   
    [statsCellGS,statsCellFG, statsCellBG] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    plusTipTestDistrib(groupData,saveDir,0.005,20,21);
    
    meanValueGS = statsCellGS(2:(length(regionTypes)+1),4);
    meanValueGL = statsCellGS(2:(length(regionTypes)+1),6);
    meanValueGD = statsCellGS(2:(length(regionTypes)+1),8);
    
    if iWindow == 1
       forPlot.meanValueGS = [regionTypes meanValueGS];
       forPlot.meanValueGL = [regionTypes meanValueGL];
       forPlot.meanValueGD = [regionTypes meanValueGD]; 
       
    else
        
       forPlot.meanValueGS = [forPlot.meanValueGS meanValueGS];
       forPlot.meanValueGL = [forPlot.meanValueGL meanValueGL]; 
       forPlot.meanValueGD = [forPlot.meanValueGD meanValueGD]; 
    
    end
    
    
    
    
    
    
    

end % end for iWindows

saveDir = [analSaveDir1 filesep 'GreaterThan_' num2str(numWindows*windowSize) 'uM' ];
mkdir(saveDir);

 %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM.mat']);
    
    groupList = groupList.groupList;
    %Perform Pooling of Data
    [groupData] = plusTipExtractGroupData(groupList,removeBeginEnd); 
    save([saveDir filesep 'groupData'],'groupData'); 
    
    plusTipPoolGroupData(groupData, saveDir, doWtn, doPlot); 
    
   
    
    % Collect Data
   
    
    
   [ statsCellGS, statsCellFG, statsCellBG, stats ] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    
    %Perform Statistical Tests
    
   plusTipTestDistrib(groupData,saveDir,0.005,20,21); 
    
    
    meanValueGS = statsCellGS(2:(length(regionTypes)+1),4);
    meanValueGL = statsCellGS(2:(length(regionTypes)+1),6);
    meanValueGD = statsCellGS(2:(length(regionTypes)+1),8);
    
       forPlot.meanValueGS = [forPlot.meanValueGS meanValueGS];
       forPlot.meanValueGL = [forPlot.meanValueGL meanValueGL]; 
       forPlot.meanValueGD = [forPlot.meanValueGD meanValueGD]; 
    
    
    
    save([statDir filesep 'forPlotting'],'forPlot');
    
    
    
    
    
    
%Do Same For Spatial Comparisons Within Region Type
    
    % Pool the data from the groupList


for iRegionType= 1:length(regionTypes)
    
    % Collect Data
   
    saveDir = [analSaveDir2 filesep char(regionTypes{iRegionType})];
    mkdir(saveDir);
    
    


    %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareWindows_' char(regionTypes{iRegionType,1}) '.mat']);
    groupList = groupList.groupList;
    
    
    %Perform Pooling of Data
    groupData = plusTipExtractGroupData(groupList,removeBeginEnd); 
    save([saveDir filesep 'groupData'],'groupData')
    plusTipPoolGroupData(groupData, saveDir, doWtn, doPlot); 
    
    
    [ statsCellGS, statsCellFG, statsCellBG] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    
   %Perform Statistical Tests
   plusTipTestDistrib(groupData,saveDir,0.005,20,21); 
    
    
end % end for iRegionType

        
end % end if doAnalysis
    


     



cd(homeDir)
warning(warningState);
disp('Sub-ROIs...finished')

end 




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
end