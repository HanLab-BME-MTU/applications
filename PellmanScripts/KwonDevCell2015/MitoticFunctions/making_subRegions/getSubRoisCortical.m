function [ projListMod] = getSubRoisCortical(projList,edgeDist,subRoiType,subRoiFilename)
%getSubRoisCortical: makes cortical subRois from two pole mitotic points: records pole to cortex distance.
%useful for Pellman collaboration
% INPUT: 
% projList: a cell array of of directories where the whole cell analysis
% that you want to partition is stored. 
% edgeDist: the distance from the edge (in um) that you would like to make
% the cortical subRoi 
% subRoiType: 1 if bipolar, 2 if full cortical 
% OUTPUT:            
% projListMod: if subRoiType = 1 a modified cell array of directories that
%             correspond only to bipolar cells 
%             if subRoiType =2 (full cortical) it will return the same of
%             the original project list

%% Make subRois is input requires
if subRoiType ==1;
    
    for iProj = 1:numel(projList)
        anDir = char(projList(iProj));
       % anDir = projList(iProj).anDir; 
     
        %     anDir = strrep(anDir,'Mijung','Mijung_RPE_MitoticSegmented');
        % imDir = projList(iProj).imDir;
        %     imDir = strrep(imDir,'Mijung','Mijung_RPE_MitoticSegmented');
        %anDir = projList{iProj};
        % check for number of poles
        %s  = load([anDir filesep 'poleInfo.mat']);
        s.poleInfo.numPoles =2; % just set = to 2 for 11-24 data as already screened the set for bipolar
        
        % if polar
        %if subRoiType == 1 % polar
        % test to make sure bipolar
        %%
        if s.poleInfo.numPoles ==2
            
          if exist([ anDir filesep subRoiFilename filesep 'sub_2' filesep 'masks'],'file')==7; 
           % check inside 
              maskTest= searchFiles('.tif',[],[ anDir filesep subRoiFilename filesep 'sub_1' filesep 'masks'],0); 
            [test,~] = size(maskTest); 
         
            if test == 100
                display( [projList{iProj} ' Run Previously SKIPPING']); 
            end % if test ; 
          else 
              
          
              
          
             
                
          
         
                display( ['Getting SubRois for ' projList{iProj}]); 
                
          
            % proceed
            projListFlag(iProj,1) = 1;
            
            x = getFilenameBody(anDir);
            imDir = [x filesep 'images'];
            [up2,con,cellNum] = getFilenameBody(x);
            [up3,~] = getFilenameBody(up2);
            
            % load mask internal
            % maskExternal = logical(imread([anDir filesep 'roiMaskSeg.tif']));
           
            maskDir = [anDir filesep 'masks'];
            
            % load masks
            [masks] = searchFiles('.tif',[],maskDir,0);
            mask1 = [masks{1,2} filesep masks{1,1}];
            [path body no ext ] = getFilenameBody(mask1);
            
            % if length(no)>1
            %     padded = 1;
            % else
            %     padded = 0;
            % end
            %
            % if padded == 0
            %    [ sortedMaskList, sortednum] = sortUnpaddedList(masks);
            % end
            
            [listOfImages] =  searchFiles('.tif',[],imDir,0);
            img1 = [listOfImages{1,2} filesep listOfImages{2,2}];
            [~,~,no,~] = getFilenameBody(img1);
            if length(no)>1
                padded2 = 1;
            else padded2 = 0 ;
            end
            if padded2 == 0;
                [sortedImages , sortednum2]= sortUnpaddedList(listOfImages);
            end
            
            nImTot = length(sortedImages);
            fmt = '%03d';
            
             if exist([anDir filesep 'roiMaskKinetochore.tif'])==0 
          
                mask = logical(imread(mask1)); 
                
                maskInternal = zeros(size(mask)); 
             else 
                   maskInternal = logical(imread([anDir filesep 'roiMaskKinetochore.tif']));
            
            end 
            %% Run through all frames and make masks
            
            for iFrame = 1:length(sortedImages)
%                 fileNameIm = [char(sortedImages(iFrame,1)) filesep char(sortedImages(iFrame,2)),...;
%                     num2str(sortednum2(iFrame)) char(sortedImages(iFrame,4))];
                  fileNameIm = [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}]; 
                img =  double(imread(fileNameIm));
                
                %     if isempty(maskInternal)
                %         maskInternal = roipoly;
                %     end
                fileNameMask = [masks{iFrame,2} filesep masks{iFrame,1} ];
                
                %     fileNameMask  = [char(sortedMaskList(iFrame,1)) filesep char(sortedMaskList(iFrame,2)),...;
                %         num2str(sortednum(iFrame)) char(sortedMaskList(iFrame,4))];
                maskExternal = double(imread(fileNameMask));
                maskExtraCell = ~maskExternal;
                maskExternal = logical(getLargestCC(maskExternal));
                
                
                
                mask = maskExternal-maskInternal;
                if sum(maskInternal(:) ) ~= 0 
                % get long centroid
                figure('Visible','off')
                
                imshow(img,[]);
                props = regionprops(maskInternal,'centroid','orientation');
                centroid = props.Centroid;
                end 
                %     angle =- props.Orientation;
                %
                % for iPt = 1:2
                % h=msgbox('click the centrosome','help');
                % uiwait(h);
                % h=impoint;
                % position = wait(h);
                % ptMask = zeros(size(img));
                % idx = sub2ind(size(img),round(position(2)),round(position(1)));
                % point(iPt,2) = position(2);
                % point(iPt,1) = position(1);
                % ptMask(idx) = 1;
                % poles(:,:,iPt) = ptMask;
                %
                %
                % end
                %
                % close gcf
                
                s = load([anDir filesep 'poleInfo.mat']);
                % s = load([anDir filesep 'subRois_Poles2' filesep 'poles.mat']);
                %   poles = s.poles;
                %
                % pole1= poles(:,:,1);
                %  pole2 = poles(:,:,2);
                %  [y1,x1] = ind2sub(size(img),find(pole1==1));
                %  [y2,x2] = ind2sub(size(img),find(pole2==1));
                %  point = [x1,y1;x2,y2];
                
                point = s.poleInfo.coords;
                poles = s.poleInfo.poleMasksInd;
               %poles = s.poleInfo.poleMasks; % field in 2/11 data  
               mAxis = (point(1,2)-point(2,2))/(point(1,1)-point(2,1));
                angle = atand(mAxis);
                bAxis = point(1,2) - mAxis*point(1,1);
                
                mAxisPer = -1/mAxis;
                anglePer = atand(mAxisPer);
                bAxisPer = centroid(2)-mAxisPer*centroid(1);
                
                % full perpendicular ROIS for Assymetry
                [imL, imW ] = size(mask);
                x = [1:imW];
                y = mAxis +bAxis;
                
                [~,yAll]=meshgrid(1:imW,1:imL);
                yLineAxis=repmat(mAxis.*[1:imW]+bAxis,[imL,1]);
                yPer = repmat(mAxisPer.*[1:imW]+bAxisPer,[imL,1]);
                
                r11 = yAll<=yPer;
                r12 = yAll>yPer;
                
                
                
                % this will be used to find intersection and test for long axis label
                x = [1:imW] ;
                y = mAxis.*x+bAxis;
                x = x(y<imL);
                y = y(y<imL);
                x = x(y>0);
                y= y(y>0);
                pixelTest = sub2ind(size(img),ceil(y),x);
                lineMask = zeros(size(mask));
                lineMask(pixelTest) = 1;
                
                % sort the poles to correspond to one pole or the other.
                % pole test
                wholeMask = imfill(mask);
                test1 = r11& wholeMask;
                test2 = r12& wholeMask;
                idxtest1 =  find(test1==1);
                idxpole1 = find(poles(:,:,1)==1);
                overlap= intersect(idxtest1,idxpole1);
                % switch the regions
                if isempty(overlap)
                    r11 = yAll>yPer;
                    r12 = yAll<=yPer;
                end
                %% Asymmetry SubRois
                %     m1 =  tand(anglePlusMinus + angle);
                %     m2 = tand(angle -anglePlusMinus);
                %     b1 = centroid(2) - m1*centroid(1);
                %     b2 = centroid(2) - m2*centroid(1);
                %
                %
                %
                
                %
                %      x = [1:imW];
                %      y = mAxis +bAxis;
                %     [~,yAll]=meshgrid(1:imW,1:imL);
                %
                %     yLineAxis=repmat(mAxis.*[1:imW]+bAxis,[imL,1]);
                % %     y1 = repmat(m1.*[1:imW]+b1,[imL,1]);
                % %     y2 = repmat(m2.*[1:imW]+b2,[imL,1]);
                %     line = yAll==yLineAxis;
                %
                %     r11 = yAll<=y1;
                %     r12 = yAll>y1;
                %
                %     r21 = yAll<= y2;
                %     r22 = yAll>y2;
                %
                %     roiSetTest(:,:,1)= r11 & mask & r22;
                %     roiSetTest(:,:,2) = r12 & mask & r21;
                %     roiSetTest(:,:,3) = r11 & mask & r21;
                %     roiSetTest(:,:,4)  = r12 & mask & r22;
                %
                %     %First test for long axis label;
                
                %
                %     countLong = 1;
                %     countShort = 3;
                %     for i = 1:4
                %         CCTest = bwconncomp(roiSetTest(:,:,i));
                %       First  pixels = CCTest.PixelIdxList{1};
                %         if ~isempty(intersect(pixelTest,pixels));
                %
                %             roiSet(:,:,countLong) = roiSetTest(:,:,i);
                %         countLong = countLong+1;
                %         else
                %             roiSet(:,:,countShort) = roiSetTest(:,:,i);
                %             countShort = countShort+1;
                %         end
                %      end
                %     figure('Visible','on');
                %
                %     imshow(img,[])
                %     hold on
                %     plot(x,y,'b');
                %      colors = ['r','b','g','y'];
                %     for iRoi = 1:4
                %
                %         currentSubRoiDir = [anDir filesep 'subRois_Poles' filesep 'sub_' num2str(iRoi)];
                %         if isdir(currentSubRoiDir)
                %             rmdir(currentSubRoiDir,'s');
                %         end
                %         mkdir([currentSubRoiDir filesep 'meta'])
                %         mkdir([currentSubRoiDir filesep 'feat'])
                %         movieInfo1  = [anDir filesep 'feat' filesep 'movieInfo.mat'];
                %         movieInfo2 = [currentSubRoiDir filesep 'feat' filesep 'movieInfo.mat'];
                %         copyfile(movieInfo1,movieInfo2);
                %         roiMask = roiSet(:,:,iRoi);
                %         imwrite(roiMask,[currentSubRoiDir filesep 'roiMask.tif']);
                %
                %         [y1,x1]=ind2sub([imL,imW],find(roiSet(:,:,iRoi),1)); % first pixel on boundary
                %         roiYX = bwtraceboundary(roiSet(:,:,iRoi),[y1,x1],'N');
                %         save([currentSubRoiDir filesep 'roiYX.mat'],'roiYX');
                %
                %
                %         plot(roiYX(:,2),roiYX(:,1),colors(iRoi));
                %
                %         % make weighted mask using distance transform to find position
                %         % where text should go
                %         weightedRoi=bwdist(~roiSet(:,:,iRoi));
                %         [r,c]=find(weightedRoi==max(weightedRoi(:)));
                %         text(c(1),r(1),num2str(iRoi),'color','r','fontsize',14);
                %
                %
                %    %   plusTipSubRoiExtractTracksFORPACKAGE(currentSubRoiDir,'turnFiguresOn',0);
                %
                %     end
                %
                %
                %        spy(poles(:,:,1)|poles(:,:,2),25,'y');
                %        scatter(centroid(1),centroid(2),10,'c');
                %     if ~isdir([saveDir filesep 'CollectedSubRois' filesep 'Poles']);
                %         mkdir([saveDir filesep 'CollectedSubRois' filesep 'Poles']);
                %
                %     end
                %     saveas(gcf,[saveDir filesep 'CollectedSubRois' filesep 'Poles' filesep ['subRois_Cell_' num2str(cellNum) '.tif']]);
                %      saveas(gcf,[anDir filesep 'subRois_Poles' filesep 'subRoiPoles.eps'],'psc2');
                %      save([anDir filesep 'subRois_Poles' filesep 'poles.mat'],'poles');
                %
                %
                %     close(gcf)
                % %   name = ['HSET_M10P_Db' num2str(cellNum) '_w1DIC_t1.TIF'];
                % %   name = strrep(name,'b0','b');
                % %    DIC =  double(imread([up2 filesep 'HM_Db_DIC_other_files' filesep 'images' filesep name]));
                %    figure('visible','off');
                % %    imshow(DIC,[]);
                %    hold on
                %    for iRoi = 1:4
                %          [y1,x1]=ind2sub([imL,imW],find(roiSet(:,:,iRoi),1)); % first pixel on boundary
                %         roiYX = bwtraceboundary(roiSet(:,:,iRoi),[y1,x1],'N');
                %          plot(roiYX(:,2),roiYX(:,1),colors(iRoi));
                %
                %         % make weighted mask using distance transform to find position
                %         % where text should go
                %         weightedRoi=bwdist(~roiSet(:,:,iRoi));
                %         [r,c]=find(weightedRoi==max(weightedRoi(:)));
                %         text(c(1),r(1),num2str(iRoi),'color','r','fontsize',14);
                %    end
                %    if ~isdir([saveDir filesep 'CollectedSubRois' filesep 'PolesDIC'])
                %        mkdir([saveDir filesep 'CollectedSubRois' filesep 'PolesDIC'])
                %    end
                %
                %    saveas(gcf,[saveDir filesep 'CollectedSubRois' filesep 'PolesDIC' filesep ['subRois_Cell_DIC' num2str(cellNum) '.tif']]);
                %
                %    close(gcf)
                %
                %     clear roiSet roiSetTest poles;
                %% SET-UP CORTICAL SUBREGION FOLDERS
                
                %%%%% MAKE SUBROIS %%%%%
                % find initial distance transform from cell edge (1 um)
                subRoiEdgeMask = imfill(mask);
                distTrans = bwdist(~subRoiEdgeMask);
                
                forIntersectAll = subRoiEdgeMask;
                forIntersectAll(distTrans>1) = 0;
                
                
                distTrans = distTrans.*.108; %  convert to microns
                subRoiEdgeMask(distTrans>edgeDist )=0;
                
                
                
                roiSet(:,:,1) = subRoiEdgeMask & r11;
                roiSet(:,:,2) = subRoiEdgeMask & r12;
                forIntersect(:,:,1) = forIntersectAll & r11;
                forIntersect(:,:,2) = forIntersectAll & r12;
                
                
                for iSub = 1:2
                    
                    subRoiEdgeDir = [anDir filesep subRoiFilename filesep 'sub_' num2str(iSub)];
                    if iFrame == 1 % set up new subRoiDirectory
                        if ~isdir(subRoiEdgeDir)
                            mkdir(subRoiEdgeDir)
                        end
                        mkdir([subRoiEdgeDir filesep 'meta'])
                        mkdir([subRoiEdgeDir filesep 'feat'])
                        movieInfo1  = [anDir filesep 'feat' filesep 'movieInfo.mat'];
                        movieInfo2 = [subRoiEdgeDir filesep 'feat' filesep 'movieInfo.mat'];
                        copyfile(movieInfo1,movieInfo2);
                        
                        if ~isdir([subRoiEdgeDir filesep 'masks'])
                            mkdir([subRoiEdgeDir filesep 'masks'])
                        end
                        
                        if ~isdir([subRoiEdgeDir filesep 'extraCell']);
                            mkdir([subRoiEdgeDir filesep 'extraCell']);
                        end
                        
                    end % iFrame ==1
                    % save the mask
                    fmt = '%03d';
                    imwrite(roiSet(:,:,iSub),[subRoiEdgeDir filesep 'masks' filesep 'roiMask' num2str(iFrame, fmt) '.tif']);
                    imwrite(maskExtraCell,[subRoiEdgeDir filesep 'extraCell' filesep 'roiMaskEC' num2str(iFrame,fmt) '.tif']);
                    
                    
                    
                end % iSub
                
                %%% REGION-OVER-VIEW SANITY PLOT %%%
                % search for DIC file using list images
%                 [DICs] = searchFiles('DIC',[],anDir,0);
%                 
%                 imgDIC = double(imread([char(DICs(1,2)) filesep char(DICs(1,1))]));
%                 
%                 
               toPlot(:,:,1) = img;
%                 toPlot(:,:,2) = imgDIC;
             names{1} = 'EB_Comet';
%                 names{2} = 'DIC';
              for iPlot = 1:1
                    figure('Visible','off');
                    imshow(toPlot(:,:,iPlot),[])
                    hold on
                    plot(x,y,'r');
                    plotScaleBar(100,5,'Label','10um');
                    colors = ['b','g'];
                    for iSub = 1:2
                        roiYX = bwboundaries(roiSet(:,:,iSub));
                        roiY = roiYX{1}(:,1);
                        roiX = roiYX{1}(:,2);
                        plot(roiX,roiY,colors(iSub));
                        weightedRoi=bwdist(~roiSet(:,:,iSub));
                        [r,c]=find(weightedRoi==max(weightedRoi(:)));
                        text(c(1),r(1),num2str(iSub),'color','r','fontsize',14);
                        
%                        calculate p (distance from pole to cortex)
                        pixBound = find(forIntersect(:,:,iSub) ==1);
                        edgePt = intersect(pixBound,pixelTest);
                        while isempty(edgePt)
                            dil = zeros(size(mask));
                            dil(pixBound) = 1;
                            dil =  imdilate(dil,strel('disk',1));
                            pixBound = find(dil==1);
                            edgePt = intersect(pixBound,pixelTest);
                        end
                        if length(edgePt) >1
                            edgePt = edgePt(1); % just take the first
                        end
                        [yIntersect(iSub), xIntersect(iSub)] = ind2sub(size(img),edgePt);
                        scatter(xIntersect(iSub),yIntersect(iSub),100,colors(iSub),'filled');
                        d = sqrt((point(iSub,1)-xIntersect(iSub)).^2 + (point(iSub,2)-yIntersect(iSub)).^2);
                        p_iFrame = d.*0.108;
%                        put in a cell to save for all of them 
                  
                        if iFrame ==1
                            p{iProj}(iSub) = p_iFrame; % save the first frame for sorting % either that or can save the average
                            
                            s = load([anDir filesep 'meta' filesep 'projData.mat' ]); 
                            projData = s.projData; 
                            projData.poleCortexDist = p_iFrame; % save it here.. 
                            save([anDir filesep 'meta' filesep 'projData.mat'],'projData'); 
                            
                        end
                        text(xIntersect(iSub)+1,yIntersect(iSub)+1,['P = ' num2str(p_iFrame,2) ' um'],'color','y');
                        
                        poleCortexDist(iSub,iFrame) = p_iFrame; 
                        
                       
                    end % for iSub
%                    plot the poles
                    for iPole= 1:2
                        spy(poles(:,:,iPole),colors(iPole),30)
                    end
                    spindleL = sqrt((point(1,1)-point(2,1)).^2 + (point(1,2)-point(2,2)).^2);
                    spindleL = spindleL*0.108;
                    text(centroid(1)+1,centroid(2)+1,['L = ' num2str(spindleL,2),' um'],'color','y');
                    
                    
                    
%                     % save two places
                   saveas(gcf,[anDir filesep subRoiFilename filesep names{iPlot} '.eps'],'psc2');
%                     collectedDir = [up3 filesep 'CollectedSubRois' filesep subRoiFilename filesep 'RegionOverView' filesep  con];
                    
%                     
%                     if ~isdir([collectedDir filesep 'Cell' cellNum])
%                         mkdir([collectedDir filesep 'Cell' cellNum]);
%                     end
%                     saveas(gcf,[collectedDir filesep 'Cell' cellNum  filesep 'subRois_Cortical_' names{iPlot} '_' con '_' cellNum '_' num2str(iFrame,fmt) '.tif']);
%                     close(gcf)
%                     % save DIC as well
%                     
            end % iPlot
%                 % save the pole masks for later use likely better way to store but who
                % cares
                save([anDir filesep subRoiFilename  filesep 'poles.mat'],'poles');
                
                clear roiSet forIntersectAll r11 r12 forIntersect toPlot
            end
          %save([anDir filesep 'poleCortexDist.mat'],'poleCortexDist'); % make sure to save the calculation for pole cortex dist for each frame
            save([anDir filesep subRoiFilename filesep 'poleCortexDist'],'poleCortexDist'); 
          end % ifexist 
        else
            projListFlag(iProj,1) = 0; % don't perform analysis
            
        end   % if pole info
        clear toPlot poleCortexDist
           
    end % for iProj
    % make new projList
    projListMod = projList(logical(projListFlag));
    %p = p(logical(projListFlag));
    save('projListBipolar.mat','projListMod');
    
    
    
elseif subRoiType == 2 % just do full  % a little repetitive but good enough for now
    
    %
    for iProj = 1:numel(projList)
        anDir = char(projList(iProj));
        
        x = getFilenameBody(anDir);
        imDir = [x filesep 'images'];
        [up2,con,cellNum] = getFilenameBody(x);
        [up3,~] = getFilenameBody(up2);
        
        subRoiEdgeDir = [anDir filesep subRoiFilename filesep 'sub_1'];
        if exist([subRoiEdgeDir filesep 'meta' filesep 'projData.mat'])~=0; 
            display(['Skipping ' anDir 'Already Run']);
        else 
        if ~isdir(subRoiEdgeDir)
            mkdir(subRoiEdgeDir)
        end
        mkdir([subRoiEdgeDir filesep 'meta'])
        mkdir([subRoiEdgeDir filesep 'feat'])
        movieInfo1  = [anDir filesep 'feat' filesep 'movieInfo.mat'];
        movieInfo2 = [subRoiEdgeDir filesep 'feat' filesep 'movieInfo.mat'];
        copyfile(movieInfo1,movieInfo2);
        
        
        
        
        
        
        % SEARCH FOR MASKS AND LOAD IF NECESSARY
        if exist([anDir filesep 'masks'],'dir') ==7; 
            listOfMasks = searchFiles('tif',[],[anDir filesep 'masks'],0,'all',1); % . 
        start = 1; 
        stop = size(listOfMasks,1); 
        maskOutDir = [subRoiEdgeDir filesep 'masks']; % will make a mask file 
        if isdir(maskOutDir)
           rmdir(maskOutDir,'s')       
        end 
            mkdir(maskOutDir)
        
        
        multMask = 1;     
            %maskExternal = logical(imread(listOfMasks{iFrame})); 
        else 
            projListCS.anDir = anDir; 
            imDir = upDirectory(anDir,1); 
            projListCS.imDir = [imDir filesep 'images']; 
            wrapperThresholding(projListCS,'Rosin'); 
            %make sure to read in the masks 
             listOfMasks = searchFiles('tif',[],[anDir filesep 'masks'],0,'all',1); % 
          start = 1; 
         % stop = 1;
         stop = size(listOfMasks,1); 
%          listOfMasks{1} = [anDir filesep 'roiMaskSeg.tif'];
         multMask = 1;  
        %maskExternal = logical(imread([anDir filesep 'roiMaskSeg.tif']));
         maskOutDir = [subRoiEdgeDir filesep 'masks']; % will make a mask file 
        if isdir(maskOutDir)
           rmdir(maskOutDir,'s')       
        end 
            mkdir(maskOutDir)
        end
        
        % load the kinetochore files if applicable
        if exist([anDir filesep 'roiMaskKinetochore.tif'])==0 
            maskInternal = zeros(size(maskExternal)); 
        else 
           
       maskInternal = logical(imread([anDir filesep 'roiMaskKinetochore.tif']));
        end 
        
        for iFrame = start:stop
        %[listOfImages] = searchFiles('.tif',[],imDir,0);
        %imageName1 = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
        %img =  double(imread(imageName1));
        maskExternal = logical(imread(listOfMasks{iFrame})); 
        mask = maskExternal-maskInternal;
        
        subRoiEdgeMask = imfill(mask);
        distTrans = bwdist(~subRoiEdgeMask);
        distTrans = distTrans.*.108; %  convert to microns
        subRoiEdgeMask(distTrans>edgeDist )=0;
        
        roiMask = subRoiEdgeMask;
        
        if multMask ==1; 
         imwrite(roiMask,[maskOutDir filesep 'roiMask' num2str(iFrame,'%03d') '.tif']);    
        else
            imwrite(roiMask,[maskOutDir filesep 'roiMask.tif']); 
        end 
        
        %imwrite(roiMask,[subRoiEdgeD 'roiMask.tif']);
        end    % iFrame 
        
        plusTipSubRoiExtractTracksFORPACKAGE_multMasks(subRoiEdgeDir,'fraction',0.01,'extractType','multMasks1'); 
        end  
    end % iProj
    
    
    
    projListMod = projList;
    
else % don't need to make subRois
    
    projListMod = projList;
end % if subRoiType
% %%
% 