function trackCells(hObject,projNum)
% trackCells finds and links coordinates in a serie of phase contrast images 
%
% SYNOPSIS       trackCells(hObject,projNum)
%
% INPUT          hObject : a handle to the Gui which called the function
%                projNum : which job is currently being dealt with
%
% OUTPUT         All outputs are written directly to disk 
%                M : described in trackLinker
%                cellProps : cellProps(:,1)=coord(:,1);
%	 	             cellProps(:,2)=coord(:,2);
%		             cellProps(:,3)=belongsto(:);  (number of cluster - label)
%		             cellProps(:,4)=numberOfOccurences(:);  (how many cells in the cluster this cell is in)
%		             cellProps(:,5)=bodycount(:);  (area of the cluster with the number given in belongsto)
%		             cellProps(:,6)=perimdivare(:);  (cluster)
%                binaryImage is the binary image of the areas occupied by cells                
%
% DEPENDENCIES   trackCells uses { imClusterSeg
%				   trackLinker
%				   checkMinimalCellCell
%				   findNucloiTrack
%				   body
%				   ptFindHalos
%				   templfindertrack }
%                                  
%                trackCells is used by { PolyTrack }
%
% REMARK         trackCells fetches directly in GUI PolyTrack:
%                   handles : structure with information used within GUI
% 				  from handles.jobs(projNum):
% 					imagedirectory : where are the images 
% 					imagename : what are the images called
% 					imagenameslist : list of images within imagedirectory with imagename
% 					firstimage : which images shall we start with (refers to imagenameslist)
% 					lastimage : which images shall be the last one (refers to imagenameslist)
% 					intensityMax : highest value image can have (calc from bitdepth)
% 					fi_nucleus / la_nucleus : nucloi intensity first / last image
% 					fi_background / la_background : background intensity first / last image
% 					fi_halolevel / la_halolevel : halo intensity first image
% 					leveladjust : factor to adjust intensity difference nucloi/ background
% 					minsize : minimal size of nucloi 
% 					maxsize : maximal size of nucloi
% 					minsdist : minimal distance between two cells
% 					minmaxthresh : onoff - should minima and segmentation be used 
% 					clustering : onoff - should clustering be used
% 					increment : image to image increment (imagenameslist)
% 					noiseparameter : used to calculate threshold whithin templfindertrack for ignoring certain pixels
% 					savedirectory : where shall calculated data be saved to
% 					sizetemplate : what size should a template have
% 					boxsize : what size should the searcharea(template matching) have
% 					timestepslide : hao many timesteps should retrospectively be analysed during programm execution
% 					minedge : minimal distance to edge to still use templatesearch
% 					mincorrqualtempl : minimal quality of correlation(template) 
% 					mintrackcorrqual : minimal quality of correlation (track comparison)
%
% Colin Glass, Feb 04

%%%% NOTE: it might be a good idea to change the input from a guihandle and
%%%% a number to a structure activejob, which includes the jobvalues

% Fetch and extract values from the 
handles = guidata(hObject);

imageDirectory = handles.jobs(projNum).imagedirectory;
imageName = handles.jobs(projNum).imagename;

increment = handles.jobs(projNum).increment;
firstImaNum = handles.jobs(projNum).firstimage;
lastImaNum = handles.jobs(projNum).lastimage;
imageNamesList = handles.jobs(projNum).imagenameslist;

levNucFirst = handles.jobs(projNum).fi_nucleus;
levBackFirst = handles.jobs(projNum).fi_background;
levHaloFirst = handles.jobs(projNum).fi_halolevel;      

levNucLast = handles.jobs(projNum).la_nucleus;
levBackLast = handles.jobs(projNum).la_background;
levHaloLast = handles.jobs(projNum).la_halolevel;      

% Range in which direct assignement ,of found coordinates, from frame to
% frame is accepted
maxSearch = handles.jobs(projNum).maxsearch;

saveDirectory = handles.jobs(projNum).savedirectory;

percentBackground = handles.jobs(projNum).noiseparameter;
sizeTemplate = handles.jobs(projNum).sizetemplate;
boxSizeImage = handles.jobs(projNum).boxsize;

% Minimal/maximal size of the black spot in the cells
minSizeNuc = handles.jobs(projNum).minsize;
maxSizeNuc = handles.jobs(projNum).maxsize;

% Get the maximum image intensity of the images in the job
maxImageIntensity = handles.jobs(projNum).intensityMax;

% An educated guess of the minimal distance between neighbouring cells
% (better to big than to small, for this value is used for the static
% search. Searches with templates are more tolerant)
minDistCellToCell = handles.jobs(projNum).minsdist; 
% This value influences the level between the nucloi and the background, as
% calculated from input (clicking in the pictures, which pop up shortly
% after the programm starts rolling
levelChanger = handles.jobs(projNum).leveladjust;

% Minimum 4
howManyTimeStepSlide = handles.jobs(projNum).timestepslide;

% Minimal distance to edge for tracking with template
minEdge = handles.jobs(projNum).minedge;

minimalQualityCorr = handles.jobs(projNum).mincorrqualtempl;
minTrackCorr = handles.jobs(projNum).mintrackcorrqual;

segmentation = handles.jobs(projNum).minmaxthresh;
clustering = handles.jobs(projNum).clustering;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following parameters should not be changed!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How much the blobs found in ptFindHalos shall be eroded. This is an indirect
% size criteria for ptFindHalos. Increase - minimal size of halos will be
% increased, decrease - ... decreased
erodeDiskSize = round((sqrt(minSizeNuc))/2) ;

% Range within witch to look for correlations of tracks over several pics. 
distanceCorrel = maxSearch * 1.5; 

% Within which distance shall the program body look for an area, to which
% it can allocate a cell (coordinates of a cell)
approxDistance = round(minDistCellToCell / 2);

% The following markers identify the origins of coordinates. They get added to the
% values of coordinates. Depending on how the coordinate was found, a certain marker 
% will be added. In this way, the coordinates themselves always carry all information.

% If you wish, you can add other markers, to distinguish between other classes of cells.
% DO NOT USE 0.9 AS A MARKER!!! (x <= 0.8 is ok), because the rounding will change the 
% coordinates themselves (which is not good)
% If you wish to mark cells with more then one property, markers placed at the second 
% digit (0.01 , 0.02 ...) could be used.
goodCellMarker = 0.1;
goodCell = round(goodCellMarker * 10);

templateCellMarker = 0.2;
templateCell = round(templateCellMarker * 10);

newCellMarker = 0.3;
newCell = round(newCellMarker * 10);

newCellTemplateMarker = 0.4;
newCellTemplate = round(newCellTemplateMarker * 10);

fileExtension = imageName(end-3:end);
bodyFileName = imageName(1:end-7);

howManyImg = lastImaNum - firstImaNum + 1;

levDiffFirst = abs(levNucFirst - levBackFirst) * levelChanger;
levDiffLast = abs(levNucLast - levBackLast) * levelChanger;

incrDiff = (levDiffLast - levDiffFirst) / (howManyImg - 1);
incrBack = (levBackLast - levBackFirst) / (howManyImg - 1);
incrHalo = (levHaloLast - levHaloFirst) / (howManyImg - 1);
incrNuc = (levNucLast - levNucFirst) / (howManyImg - 1);

i = 0;
k = 0;
a = 0;
lostOnEdge = 0;

newImg = [];
MPM = [];
cellProps = [];
M = [];
coord = [];

emptyM = zeros(1,4);

% Go to the directory where the images are
cd (imageDirectory);

% Check that the first and last image numbers are actually the right way around
if ~ (lastImaNum > firstImaNum)
   fprintf (1, 'Last image # is smaller than first image # in job %s\n', projNum);
   return
else    
   % Index is equal to the number of the image, countLoops keeps
   % track of the loops
   countLoops = 0;
    
   % Loop through all the images of the job using the specified increment
   for jImageNum = firstImaNum:increment:lastImaNum
        
      % Let the user know where he/she is in the sequence by printing the number
      countLoops = countLoops + 1;
      fprintf (1, 'Image number = %d, job number = %d\n', jImageNum, projNum);
        
      % Make sure we start in the image directory
      cd (imageDirectory);

      % Get the filename of the image with number jImageNum
      fileName = char(imageNamesList(jImageNum));
        
      % Read the current image and normalize the intensity values to [0..1]
      newImg = imreadnd2(fileName, 0, maxImageIntensity);
              
      % Calculate the size of the normalized image
      [img_h,img_w] = size(newImg);
        
      % If image is the first one in the sequence do the following:
      if jImageNum == firstImaNum
         % Prepare the space for the binary image
         binaryImage = zeros(img_h, img_w);
                    
         % Determine the intensity level of the halos
         haloLevel = (levHaloFirst - levBackFirst) * 2 / 3 + levBackFirst;
                      
         % In case the clustering algorithm has been selected use the kmeans method
         % to segment the image
         if clustering
            % Here's the call to the image segmentation function
            [seg_img, obj_val, nothing, mu0] = imClusterSeg(newImg, 0, 'method', 'kmeans', ...
               'k_cluster', 3, 'mu0', [levNucFirst;levBackFirst;levHaloFirst]);
                        
            % Find areas that are really dark and match cells into them
            [coordNuc, regmax] = findNucloiTrack(seg_img, levDiffFirst, minSizeNuc, maxSizeNuc, 1);
                        
            % Find cells that look like the third eye (round, big spots of
            % pure light). We do this because the pictures are of poor
            % quality and display huge halos around certain cells
            [haloCoord, logihalo] = ptFindHalos(seg_img, erodeDiskSize, haloLevel, 1);
	                   
         % Do the following when the segmentation algorithm has been selected (only
         % difference is the last parameter in the function call: 2
         elseif segmentation
            % Find areas that are really dark and match cells into them
            [coordNuc, regmax] = findNucloiTrack(newImg, levDiffFirst, minSizeNuc, maxSizeNuc, 2);

            % And again find the cells that look like big bright spots
            [haloCoord, logihalo] = ptFindHalos(newImg, erodeDiskSize,haloLevel, 2);

            % If none of the methods is selected (which should be impossible) show error
         else
            fprintf(1, 'One of the methods clustering or segmentation has to be selected.\n');
            return
         end
                      
         if ~isempty(handles.jobs(projNum).coordinatespicone);
            newCoord = handles.jobs(projNum).coordinatespicone;
         else 
            newCoord = checkMinimalCellCell(coordNuc,haloCoord,minDistCellToCell);           
         end
                      
         cellProps = [];
                     
         if clustering
            [cellProps,binaryImage,labeled] = body(seg_img,newCoord,regmax,logihalo,approxDistance,1);
         elseif segmentation
            [cellProps,binaryImage,labeled] = body(newImg,newCoord,regmax,logihalo,approxDistance,2);
         end
                     
         clear regmax;
         clear haloCoord;
         clear coordNuc;
         clear namesnumbers;
         clear uff;
                    
         emptyprop = zeros(1,size(cellProps,2));
                  
         %now we mark these coordinates as good coordinates, meaning their
         %existence is credible. Futhermore cells marked as good cells will
         %be treated with priority, when it comes to allocating cells found in
         %a new picture
         oldCoord = newCoord+goodCellMarker;
                  
         emptyM = zeros(1,4);
                      
       %the following is what we do with pictures that aren't number one
      else
         %adjust level via increment
         levDiffInterpol = levDiffFirst+(countLoops-1)*increment*incrDiff;
         levBaInterpol = levBackFirst+(countLoops-1)*increment*incrBack;
         halointerpol = levHaloFirst+(countLoops-1)*increment*incrHalo;
         levNucInterpol = levNucFirst+(countLoops-1)*increment*incrNuc;
         haloLevel = (halointerpol-levBaInterpol)*2/3+levBaInterpol;
                         
         if clustering
%            background = imopen(newImg,strel('disk',40));
% 	     imgMinusBack = imsubtract(newImg,background); 
%            smallest = min(min(imgMinusBack));
% 	     if sign(smallest)== -1
%               imgMinusBack = imgMinusBack + abs(smallest)
%            end
            [seg_img, obj_val,nothing,mu0] = imClusterSeg(newImg, 0, 'method','kmeans','k_cluster',3,'mu0', mu0);
                          
            coordNuc = [];
            [coordNuc,regmax] = findNucloiTrack(seg_img,levDiffInterpol,minSizeNuc,maxSizeNuc,1);
            [haloCoord,logihalo] = ptFindHalos(seg_img,erodeDiskSize,haloLevel,1);
	                   
         elseif segmentation
            haloLevel = (levHaloFirst-levBackFirst)*2/3+levBackFirst;
            [coordNuc,regmax] = findNucloiTrack(newImg,levDiffInterpol,minSizeNuc,maxSizeNuc,2);
            [haloCoord,logihalo] = ptFindHalos(newImg,erodeDiskSize,haloLevel,2);
         else
            disp('you have to choose one of the methods, clustering or segmentation')
            return
         end

         newCoord = checkMinimalCellCell(coordNuc,haloCoord,minDistCellToCell)  ;  

         %mark coordinates as high quality coordinates. These cells
         %are found in the usual way, so they are superior to cells
         %found for instance by templates
         newCoord = round(newCoord)+goodCellMarker;

         %we all know what the following is for - evading computer idiocity
         oldCoordflo = [];
         coordextrainf = [];
         rowcoord = [];
         firstcoord = [];
         secrowcoord = [];
         secondcoord = [];
         coord2 = [];

         %here we seperate the cells found already long ago from the ones found 
         %only lately (we believe in the existence of the first more than in the 
         %existence of the later)
         oldCoordflo = floor(oldCoord+0.1);
         coordextrainf = round(10*(oldCoord-oldCoordflo));
         rowcoord =  find(coordextrainf(:,1)==templateCell | coordextrainf(:,1)==goodCell);

         %old cells
         firstcoord = oldCoord(rowcoord,:);

         clear rowcoord;
         clear oldCoordflo;

         %new cells
         secrowcoord =  find(coordextrainf(:,1)==newCell | coordextrainf(:,1)==newCellTemplate);
         secondcoord = oldCoord(secrowcoord,:);
                       
         clear coordextrainf;
         clear secrowcoord;
                       
         %find out which cell is which over the gulf of different pictures
         %(track between always TWO pics)
                      
         %first we try to allocate the old cells to any set of coordinates
         %found in the new picture
         tmp1 = [];
         tmp1 = fsmTrackTrackerBMTNN(firstcoord,newCoord,maxSearch);
                      
         clear firstcoord;
         % clear newCoord;
                  
         % if there are new cells, now is the time to try to allocate them
         % to any set of coordinates not yet allocated to an old cell
         if size(find(tmp1(:,1)==0 & tmp1(:,2)==0),1) >0 & size(secondcoord,1) >0
            identnew = find(tmp1(:,1)==0 & tmp1(:,2)==0);
            %in any case the sets coordinates not allocated to an old cell
            %are marked as a new cells, for they will either be allocated
            %to a new cell (found a few pictures ago), or actually become a 
            %new cell
            coord2(:,1) = round(tmp1(identnew,3))+newCellMarker;
            coord2(:,2) = round(tmp1(identnew,4))+newCellMarker;

            %erase the unallocated set of coordinates from tmp1,
            %because they will turn up again in tmp2 no matter what
            %happens
            tmp1(identnew,:) = [];
            clear identnew

            tmp2 = [];
            tmp2 = fsmTrackTrackerBMTNN(secondcoord,coord2,maxSearch);
                                  
            %combine the two allocation matrices to one
            tmp = cat(1,tmp1,tmp2);
            clear tmp2;
         else 
            identnew = find(tmp1(:,1)==0);
            tmp1(identnew,3:4) = round(tmp1(identnew,3:4))+newCellMarker;
                                  
            clear identnew;
                                  
            if size(secondcoord,1) >0
               temzeros = zeros(size(secondcoord,1),4);
               temzeros(:,1:2) = secondcoord(:,1:2);
               tmp = cat(1,tmp1,temzeros);
               clear temzeros;
            else
               %since their are no new cells tmp2 is empty, so tmp1 is
               %complete
               tmp = tmp1;
            end
         end
                      
         if length( find( (tmp(:,3)~=0)&(tmp(:,4)~=0))) ~= size(newCoord,1)
            uuupsyouvelostone=1;
         end
                      
         if length( find( (tmp(:,1)~=0)&(tmp(:,2)~=0))) ~= size(oldCoord,1)
            BIGTROUBLE = 1
         end
                       
         clear tmp1;
         clear secondcoord;
         clear coord2;
                         
         %so far so good, BUT... we still have cells that weren't
         %found in the new picture. Instead of just giving up, we now
         %try to track them with templates. Of course these cells
         %will be marked (within the tracking via templates routine - called 
         %templfindertrack). If they are old cells they will be
         %marked as "old cells tracked via template", if they are new
         %cells as "new cells tracked via template". The
         %actual marker is a digit (fraction) behind the dot
         trouble = [];
         trouble = find(tmp(:,3)==0 & tmp(:,4)==0);
                           
         if ~isempty (trouble)
            for t =  1:size(trouble,1)
               tempcoord = [0,0];
               tempcoord = [tmp(trouble(t),1),tmp(trouble(t),2)];
				
               %first we check if the lost cell was located
               %near the edge (distance to edge smaller than
               %minEdge). If so, we asume the cell has wandered
               %out of the picture and we stop following it
               if abs(tempcoord(1,1))<minEdge | abs(tempcoord(1,1)-img_w)<minEdge | abs(tempcoord(1,2))<minEdge | abs(tempcoord(1,2)-img_h)<minEdge
                  %one never knows what kind of data
                  %somebody may wish
                  %to possess
                  lostOnEdge = lostOnEdge+1;
               else
                  absmaxcorr = [];
                  newcoord = [];
                    
                  %this is the finding via templates
                  %routine. It will return the found
                  %coordinates and the value of the
                  %corralation
                  [newcoord,absmaxcorr] = templfindertrack(tempcoord,oldImg,newImg,levBaInterpol,templateCellMarker,newCell,newCellTemplate,newCellTemplateMarker,percentBackground,sizeTemplate,boxSizeImage);
                                                      
                  if sqrt((tempcoord(1,1)-newcoord(1,1)).^2 + (tempcoord(1,2)-newcoord(1,2)).^2) > 2*maxSearch
                     %the programm made a mistake
                     break
                  end
                                                      
                  %of course we only accept the correlation, if it of a minimal
                  %quality. NOTE WELL: that within the routine the background (everything 
                  %that doesn't belong to a cell) of the template and of the searcharea
                  %is replaced by random noise. So a background - background correlation
                  %will average zero
                            
                  if absmaxcorr > minimalQualityCorr
                     %now we calculate the minimal distance between the cell found by
                     %correlation and any other cell
                     dista = [];
                     rowofit = [];
                     
                     [dista,rowofit] =  min(sqrt((tmp(:,3)-newcoord(1,1)).^2 + (tmp(:,4)-newcoord(1,2)).^2));
                               
                     if dista < round(minDistCellToCell/1.5)
                                             
                        %the following is just extracting the information of
                        %the coordinates found via template and the coordinates
                        %the are very close to it
                        tmpflo = [];
                        extr = [];
                        newflo = [];
                                    
                        tmpflo = floor(tmp(rowofit,3)+0.1);
                        extr = round(10*(tmp(rowofit,3)-tmpflo));
                        
                        clear tmpflo;
                       
                        newflo = floor(newcoord(1,1)+0.1);
                        extrnew = round(10*(newcoord(1,1)-newflo));
                        
                        clear newflo;
                         
                        %if the coordinates found via
                        %template belong to an old cell,
                        %we trust them, otherwise the
                        %cell won't be propagated
                        if extrnew == templateCell
                                                   
                           %if the cell close by is a new cell, we allocate it's coordinates to the 
                           %cell we try to track by template (which is 
                           %an old cell - we've checked that above) and
                           %erase the new cell
                           if extr == newCell
                              tmp(trouble(t),3) = round(tmp(rowofit,3))+goodCellMarker;
                              tmp(trouble(t),4) = round(tmp(rowofit,4))+goodCellMarker;
                              tmp(rowofit,3:4) = 0;
                                        
                               %if it is a new cell, itself propagated via
                               %templates, we erase it and allocate the
                               %coordinates we have (template) found
                               %to the old cell
                           elseif extr == newCellTemplate
                              tmp(rowofit,3:4) = 0;
                              tmp(trouble(t),3) = newcoord(1,1);
                              tmp(trouble(t),4) = newcoord(1,2);
                           end
                        end

                        clear extrnew;
                        clear extr;
                         
                        %if the minimal distance to the
                        %nearest cell is big enough,we allocate the
                        %coordinates we have (template) found
                        %to the old cell
                        tmp(trouble(t),3) = newcoord(1,1);
                        tmp(trouble(t),4) = newcoord(1,2);
                     end
                     clear dista;
                     clear rowofit;
                  end
                  clear newcoord;
                  clear absmaxcorr;
               end
               clear tempcoord;
            end
         end
                       
         %I know all these 'end' are confusing, so for orientation: we 
         %have now got the matrice (tmp), which defines the allocations between the
         %the cells of the last picture and the cells of the picture
         %we are looking at in the moment.
                       
         clear trouble;
                          
         %cosmetics
         stay = find(ismember(tmp(:,1:4),[0 0 0 0],'rows')-1);
         tmp = tmp(stay,:);
                      
         clear stay;
                     
         M(1:size(emptyM,1),1:size(emptyM,2),countLoops-1,1) = emptyM;
         %stuff this information into a stack
         M(1:size(tmp,1),1:size(tmp,2),countLoops-1,1) = tmp;
                          
         if length( find ((tmp(:,1)~=0)&(tmp(:,2)~=0))) ~= size(oldCoord,1)
            BIGTROUBLE=1;
         end
                       
         clear tmp;
         tmp = emptyM;
                       
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %this is a quite elaborate and (if I may say so) sophisticated
         %part of the programm, which basically changes the tracks of
         %certain cells, based on an analysis of the various tracks of
         %cells over the last ? frames (defined in PolyTrack (GUI)). 
                   
         if countLoops > howManyTimeStepSlide+1
                           
            templateRow = [];
            MPMslide = [];
            MPMflo = [];
            xtempl = [];
            ytempl = [];
                               
            MPMslide = trackLinker(M(:,:,(countLoops-howManyTimeStepSlide+1):(countLoops-1)));
            %MPMslide is a matrix which gives the tracks of all cells over
            %the last ? frames. It get's updated after every new frame.
                  
            MPMflo = floor(MPMslide+0.1);
            MPMIdentCell = round(10*(MPMslide-MPMflo));
            %extrinf has the same size as MPMslide and every value reflects the extracted
            %information of the coordinate located at the same spot in
            %MPMslide. The information is already included in MPMslide. It
            %is the digit after the point, which identifies the origin of a
            %coordinate. Look above, marker
                              
            clear MPMflo;
                               
            %After a few frames, new cells will be accepted as old cells, but only if, 
            %after their initial finding, they weren't propagated by templates 
            identNew = find(MPMIdentCell(:,3)==newCell);
            
            if ~isempty(identNew)
               for jnewCell = 1:length(identNew)
                  if isempty(find(MPMIdentCell(identNew(jnewCell),:)==newCellTemplate))  &   MPMIdentCell(identNew(jnewCell),end)~=0; %make sure a new cell wasn't propagated via templates
                     markOld = [];
                     markOld = find(M(:,3,countLoops-1)==MPMslide(identNew(jnewCell),end-1) &  M(:,4,countLoops-1)==MPMslide(identNew(jnewCell),end));
                     M(markOld,3:4,countLoops-1) = round(M(markOld,3:4,countLoops-1))+goodCellMarker;
                     MPMslide(identNew(jnewCell),end-1) =  M(markOld,3,countLoops-1);
                     MPMslide(identNew(jnewCell),end) =  M(markOld,4,countLoops-1);
                                          
                     clear markOld;
                                     
                  elseif  MPMIdentCell(identNew(jnewCell),end)==newCellTemplate & MPMIdentCell(identNew(jnewCell),end-2)==newCellTemplate
                     delCell = [];
	             delCell = find(M(:,3,countLoops-1)==MPMslide(identNew(jnewCell),end-1) &  M(:,4,countLoops-1)==MPMslide(identNew(jnewCell),end));
	             M(delCell,3:4,countLoops-1) = 0;
	             MPMslide(identNew(jnewCell),end-1) = 0;
                     MPMslide(identNew(jnewCell),end) = 0;

                     clear delCell
                  end
               end
            end
                  
            clear identNew;

            %here we search for old cells, that have been propagated via
            %template search over the last five frames. Afterwards we will
            %compare their tracks with tracks of newly found cells. If they
            %corralate strongly, we replace the templatestuff (old cell) by the newly
            %found cell
                
            templateRow = find(ismember(MPMIdentCell(:,3:end),templateCell*ones(1,length(MPMIdentCell(1,3:end))),'rows'));
                  
            if ~isempty(templateRow)
               newCellRow =  find(MPMIdentCell(:,3)==newCell & MPMIdentCell(:,end)==newCell);
               %newly found cells
                      
               if ~isempty(newCellRow)
                           
                  if ~isempty(find(MPMIdentCell(newCellRow,3:end)))
                     BIGISHTROUBLE = 1;
                  end
                                
                  xtempl = [];
                  ytempl = [];
                             
                  %extract the x and y
                  %coordinates of the cells         
                  xtempl(:,:) = (MPMslide(templateRow(:),3:2:end-1))'; 
                  ytempl(:,:) = (MPMslide(templateRow(:),4:2:end))';
                         
                  %we want to correlate the displacementvectors
                  %(displacement of a cell from one frame to the 
                  %next), so we need to calculate it
                  moveXXTemplate = diff(xtempl);
                  moveYYTemplate = diff(ytempl);
                                
                  for indexnewCell = 1:length(newCellRow)
                                                      
                     xnewone = [];
                     ynewone = [];
                           
                     moveXXnewCell = [];
                     moveYYnewCell = [];
                                                
                     xnewone = (MPMslide(newCellRow(indexnewCell),3:2:end-1))';
                     ynewone = (MPMslide(newCellRow(indexnewCell),4:2:end))';
                     
                     moveXXnewCell = diff(xnewone);
                     moveYYnewCell = diff(ynewone);
             
                     %make sure we correlate two things of the same length
                     if size(xtempl,1) == length(xnewone)
                        corrTempRow = [];
                        nomore = [];
                        for indexTempRow = 1:length(templateRow)
                                      
                           lostfound = [];
                           tempfound = [];
                           writeHereRow = [];
                                          
                           distance = min(sqrt((xtempl(:,indexTempRow)-xnewone(:)).^2+(ytempl(:,indexTempRow)-ynewone(:)).^2 ));
                                           
                           %now, if these points at some point are
                           %really close - exchange them
                           if distance < round(minDistCellToCell/2)
                                                         
                              %this is the substitution of the cell
                              %tracked by templates by the newly
                              %found cell. It's kind of tricky,
                              %because we have to change the values
                              %not in MPMslide, but in M. Since M
                              %is not sorted, we have to always
                              %find the right spot
                              %IRGENDWIE_USTUUUUUSCHE = IRGENDWIE_USTUUUUUSCHE +1;
                              templateCoord = [];
                              newCeCoord = [];
                                                       
                              %first we copy out the coordinates of
                              %the old cell, tracked by templates
                              templateCoord = MPMslide(templateRow(indexTempRow),:);
                              %Now we erase the values from MPMslide, so that we won't by accident use them again
                              MPMslide(templateRow(indexTempRow),:) = 0;
                                                     
                              %same procedure for the new cell we
                              %want to link to the old cell
                              newCeCoord = MPMslide(newCellRow(indexnewCell),:);
                              MPMslide(newCellRow(indexnewCell),:) = 0;
                                                  
                              for jFrameCount = 1:((size(MPMslide,2)-2)/2)
                                                          
                                 writeHereRow = [];
                                 eraseHereRow = [];
                                         
                                 %ok,now: in templateCoord is the track of the old 
                                 %cell.newCeCoord is the track of the cell we wish to put at
                                 %it's place. Both are of the form:
                                 %[x1,y1,x2,y2,x3,y3,x4,y4,x5,y5]
                                 %So always two numbers correspond to the coordinates
                                 %in one frame(1,2,3...) mean different frames.
                                 %SOOOO: we find the location of the coordinates of the old
                                 %cell, by searching elements of templateCoord in M. With every
                                 %increase in jFrameCount we increase the index of templateCoord by two
                                 %(because two values (x and y)per picture). M has three
                                 %dimensions. Every element of the third dimension is a
                                 %matrice(other two dimensions)of the form 
                                 %[x1,y1,x2,y2],thus giving the connection
                                 %between two frames, BUUUUTTTT: the next element 
                                 %will be [x2,y2,x3,y3]!!!!!! This is very important!!!!!!. 
                                 %So with every frame (jFrameCount) we make one step into the third dimesion in M.
                                 %with each set of x and y of templateCoord we locate the rigth
                                 %row in M. 
                                 writeHereRow = find(M(:,1,countLoops-5+jFrameCount)== templateCoord(2*jFrameCount-1)  &  M(:,2,countLoops-5+jFrameCount)== templateCoord(2*jFrameCount));
                                                           
                                 %since we overwrite the old
                                 %cell with the new, we have to
                                 %eliminate the new, in order
                                 %not to have it twice
                                 eraseHereRow = find(M(:,1,countLoops-5+jFrameCount)== newCeCoord(2*jFrameCount-1)  &  M(:,2,countLoops-5+jFrameCount)== newCeCoord(2*jFrameCount));
                                  
                                 %now we juggle the information into the right place.
                                 %Since M is
                                 %[x1,y1,x2,y2] , next layer [x2,y2,x3,y3]
                                 %and newCeCoord is
                                 %[x1,y1,x2,y2,x3,y3,x4,y4,x5,y5]
                                 %therefor:  
                                 %M(writeHereRow,:,countLoops-4) replaced by newCeCoord(1:4)
                                 %and
                                 %M(writeHereRow,:,countLoops-3) replaced by newCeCoord(3:6) 
                                 %and so on.
                                 M(writeHereRow,:,countLoops-5+jFrameCount) = newCeCoord(2*jFrameCount-1:2*jFrameCount+2);
                                                       
                                 %erase redundant info
                                 M(eraseHereRow,:,countLoops-5+jFrameCount) = 0;
                              end
                              clear writeHereRow;
                              clear eraseHereRow;
                              clear templateCoord;
                              clear newCeCoord;
                                              
                                              
                              nomore = 1;
                              break
                                                    
                           %here too, we have to stay within a maximal distance
                           elseif distance < distanceCorrel 
                              %now we correlate the displacement vectors of the old cell with
                              %the ones of the new cell
                              correlx = xcorr(moveXXTemplate(:,indexTempRow), moveXXnewCell(:),'coeff');
                              correly = xcorr(moveYYTemplate(:,indexTempRow),moveYYnewCell(:),'coeff');
                              
                              indCorr = round(length(correlx)/2);
                                                            
                              %one correlation value for each row of templates
                              corrTempRow(indexTempRow) = correlx(indCorr)+correly(indCorr);
                           end
                        end
                              
                        %nomore says if this new cell was already asigned to an old cell by distance criteria
                        %if this isn't the case, it's time to see, if there is an old cell (nearby as checked above)
                        %that correlates nicely
                        if isempty(nomore)
                           if ~isempty(corrTempRow);
                              indBC = [];
                              [bestCorr,indBC] = max(corrTempRow) ;
                                                   
                              if bestCorr > minTrackCorr
                                 %this is exactly the same procedure as above. I know it would be clever to make a function out of it
                                 %but there would be so many BIG matrices involved I decided to have it in here twice. Ah well. Anyway look above
                                 %(look for IRGENDWIE_USTUUUSCHE)
                                 %IRGENDWIE_USTUUUSCHE = IRGENDWIE_USTUUUSCHE +1
                                 templateCoord = [];
                                 newCeCoord = [];
                                                                      
                                 templateCoord = MPMslide(templateRow(indBC),:);
                                 MPMslide(templateRow(indBC),:) = 0;
                                                                        
                                 newCeCoord = MPMslide(newCellRow(indexnewCell),:);
                                 MPMslide(newCellRow(indexnewCell),:) = 0;
                                                                      
                                 for jFrameCount = 1:(size(MPMslide,2)-2)/2
                                    writeHereRow = [];
                                    eraseHereRow = [];
                                  
                                    writeHereRow = find(M(:,1,countLoops-5+jFrameCount)== templateCoord(2*jFrameCount-1)  &  M(:,2,countLoops-5+jFrameCount)== templateCoord(2*jFrameCount));
                                    eraseHereRow = find(M(:,1,countLoops-5+jFrameCount)== newCeCoord(2*jFrameCount-1)  &  M(:,2,countLoops-5+jFrameCount)== newCeCoord(2*jFrameCount));
                                          
                                    M(writeHereRow,:,countLoops-5+jFrameCount) = (newCeCoord(2*jFrameCount-1:2*jFrameCount+2));
                                                                 
                                    M(eraseHereRow,:,countLoops-5+jFrameCount) = 0;
                                 end
                                 clear writeHereRow;
                                 clear eraseHereRow;
                                 clear templateCoord;
                                 clear newCeCoord;
                              end
                              clear indBC;
                           end    
                        end
                        clear corrTempRow;
                        clear nomore;
                     else
                        ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRROOOOOOOOOOOOOOOOOOR = 1;
                        %if you've landed here by executing the programm your in TROUBLE. But it's basically impossible
                        %so you haven't - and you needn't worry. UNLESS you've been meddling with the programm. Tut tut.
                        %But actually it isn't THAT bad. The displacement vector of old cell and new cell (moveXXTemplate,moveXXnewCell...)
                        %aren't of the same length, so check what's wrong.
                     end
                     clear xnewone;
                     clear ynewone;
                   
                     clear moveXXnewCell;
                     clear moveYYnewCell;
                  end
               end
            end
            clear templateRow;
            clear MPMslide;
            clear xtempl;
            clear ytempl;
         end
       
         clear newCoord;
          
         oldCoord = [];
         realcoord = [];
         realcoord = find(M(:,3,countLoops-1) | M(:,4,countLoops-1));
         
         %oldCoord = newCoord;
       
         oldCoord(:,1) = M(realcoord,3,countLoops-1);%             tmp(realcoord,3);
         oldCoord(:,2) = M(realcoord,4,countLoops-1);
         clear realcoord;
          
         %run body now, for conveniance (we need binary images of 
         %the cells to run body)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         prope = [];
         ero = [];
                  
         if clustering
            [prope,binaryImage,labeled] = body(seg_img,oldCoord,regmax,logihalo,approxDistance,1);
         elseif segmentation
            [prope,binaryImage,labeled] = body(newImg,oldCoord,regmax,logihalo,approxDistance,2);
         end
                 
         cellProps(1:size(emptyprop,1),1:size(emptyprop,2),countLoops,1) = emptyprop;
         %store this information in a stack
         cellProps(1:size(prope,1),1:size(prope,2),countLoops,1) = prope;
                 
         clear prope;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
      oldImg = newImg;
          
      %we save all the data we found to some directory. It looks awefully
      %puny after so many lines of code.
      cd (saveDirectory)
		
      save('M', 'M')
      save ('cellProps', 'cellProps')
      
      cd ('body')
      % Format
      s=3;
      strg = sprintf('%%.%dd',s);

      % Create numerical index
      indxStr = sprintf(strg, jImageNum);
      nameClust = ['clusters' indxStr]; 
      save (nameClust, 'binaryImage')
        
      %save (['lostOnEdge',date], 'lostOnEdge')
   end
      
   %it's all over now. All we have to do is change the format of the
   %information we have painstakingly gathered
   %cd (saveDirectory)
   %MPM = trackLinker(M);
end
         
% AK: the comments below are kept because Colin has worked so hard on them :-)

%PLEASE NOTE WELL:
%by changing these values, 
%you have nothing to gain.
%no profit, no revenues,
%only TROUBLE that will rain down on your brain.
%By now, my friend, it should be quite plain:
%If you change these values ... you are insane.


