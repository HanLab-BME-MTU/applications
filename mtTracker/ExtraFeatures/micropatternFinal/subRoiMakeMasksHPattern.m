function [ roiSet,roiYX] = subRoiMakeMasksHPattern( currentProjDir, maskParams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% loops through a number of different types of H-pattern masks 
% specificied in the above 
% assume H pattern 
% 
% 
%
% parameters to loop through
%numWindows = [2 1]; 
%windowSize = [1 3];% vector of window sizes to loop through  

 
    %if HPattern== 1 
        angle1 = 60;
        angle2 = 30;
    %else 
        
     %   angle1 = 45;
     % angle2 = 15;
      %  angle3 = 75;
     
    %end 
if nargin < 1 || isempty(currentProjDir)
  currentProjDir =  uigetdir(pwd,'Please Select a Project Directory'); 
end 

%% loop over all mask types the user would like to make  

    
%%Initiate 
        upOneDir = getFilenameBody(currentProjDir); 
        imDir = [upOneDir filesep 'images']; 
        
        % get list of images, read first image, and make RGB (gray)
        [listOfImages]=searchFiles('.tif',[],imDir,0);
        img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
        [imL,imW]=size(img(:,:,1));
        
        
      
        
        
        %%% Start Get subRoiSets %%% 
        
        % load mask for current proj 
        p = load([currentProjDir filesep 'roiYXSeg.mat']); % load the roiYX file in the anDir
        roiYX=p.roiYX;
        roiMask=roipoly(img,roiYX(:,2),roiYX(:,1));
      
        
        % perform distanceTransform
        s = load([currentProjDir filesep 'meta' filesep 'projData.mat']); 
        pixSizMic = s.projData.pixSizeNm/1000; 
        weightedRoi=bwdist(swapMaskValues(roiMask)).*pixSizMic; %Distance Transform
        %distCutoff = maskParams.windowSize(iWin); 
        
        % Calculate the mulitple inner masks 
         
        % initialize matrix to hold inner masks 
        innerMasks = zeros(imL,imW,maskParams.numWindows);
    
            % fill matrix holding distance transform derived masks
            for iWindow = 1:maskParams.numWindows
                innerMasks(:,:,iWindow) = weightedRoi>maskParams.windowSize*iWindow;
            end 
      
           
            if maskParams.numWindows ~=0 
                
                % Initialize matrix to hold the distance cut-off rois (For
                    % all) 
                    cellEdgeRoi= zeros(imL,imW,maskParams.numWindows+1);
                
                    % Make Cell Edge Specific Regions Based on specified distance from the cell edge 
                    cellEdgeRoi(:,:,1) = roiMask - innerMasks(:,:,1);
                    for i = 1:(size(innerMasks,3)-1)
                        cellEdgeRoi(:,:,i+1) = innerMasks(:,:,i)-innerMasks(:,:,i+1);
                    end 
                    cellEdgeRoi(:,:,maskParams.numWindows+1) = innerMasks(:,:,maskParams.numWindows);
             end 
            
           if maskParams.subRegions ~= 1 
               roiSet = cellEdgeRoi ; % no subregions
           else % further divide these regions to make the appropriate masks 
          
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
                 
                % cluster maskParams.subRegions based on H-pattern
                    
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
           end   % end if  
           
    
          
%            if maskParams.subRegions == 1;
%                subRegions = 'subRegions';
%            else
%                subRegions = ''; 
%            end
%            
%          filename = ['roiSet_' num2str(maskParams.numWindows) '_'...
%              num2str(maskParams.windowSize) '_um_windows_' subRegions]; 


 %save([currentProjDir filesep dirName filesep filename ],'roiSet');
end 

