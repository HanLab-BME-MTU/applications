function [MPM,PROPERTIES,M,lostonedge]=trackmater(hObject,projNum,bool)


% % % % % % % % % % % % % % set(handles.GUI_st_bp_mmpixel_pm,'String',num2str(activeJob.mmpixel));




handles = guidata(hObject);

%projNum = get(handles.GUI_st_job_lb,'Value');




path=handles.jobs(projNum).imagedirectory;
main=handles.jobs(projNum).imagename;
ext=main(end-3:end);
bodyfilename=main(1:end-7);


First=handles.jobs(projNum).firstimage;
Last=handles.jobs(projNum).lastimage;


levnuc_fi=handles.jobs(projNum).fi_nucleus;
levback_fi=handles.jobs(projNum).fi_background;
levnuc_la=handles.jobs(projNum).la_nucleus;
levback_la=handles.jobs(projNum).la_background;
               


radius=handles.jobs(projNum).maxsearch;
%range in which direct assignement ,of found coordinates, from frame to
%frame is accepted

Increment=handles.jobs(projNum).increment;


pathbring=handles.jobs(projNum).savedirectory;


percentbackground=handles.jobs(projNum).noiseparameter;
sizetemple=handles.jobs(projNum).sizetemplate;
box_size_img=handles.jobs(projNum).boxsize;


ErodeDiskSize=10;
%how much the blobs found in halosfind shall be eroded. This is an indirect
%size criteria for halosfind. Increase - minimal size of halos will be
%increased, decrease - ... decreased

HaloLevel=handles.jobs(projNum).halolevel;

%minimal/maximal size of the black spot in the cells
minsizenuc=handles.jobs(projNum).minsize;
maxsizenuc=handles.jobs(projNum).maxsize;




MinDistCellCell=handles.jobs(projNum).minsdist;
%an educated guess of the minimal distance between neighbouring cells
%(better to big than to small, for this value is used for the static
%search. Searches with templates are more tolerant.)



LevelChanger=handles.jobs(projNum).leveladjust;
%this value influences the level between the nucloi and the background, as
%calculated from input (clicking in the pictures, which pop up shortly
%after the programm starts rolling

%at least 4!!!
howmanytimestepsslide=handles.jobs(projNum).timestepslide;


MinEdge=handles.jobs(projNum).minedge;
%minimal distance to edge for tracking with template

MinimalQualityCorr=handles.jobs(projNum).mincorrqualtempl;

MinTrackCorr=handles.jobs(projNum).mintrackcorrqual;
  






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that aren't open to discussion

DistanceCorrel=radius*1.5; 
%range within witch to look for corralations of tracks over several pics. 

PlusMinus=round(MinDistCellCell/2);
%Within which distance shall the programm body look for an area, to which
%it can allocate a cell (coordinates of a cell)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DO NOT CHANGE THESE VALUES !!!!!!!!!!!
%They have NOTHING to do with how the programm works!

%PLEASE NOTE WELL:
%by changing these values, 
%you have nothing to gain.
%no profit, no revenues,
%only TROUBLE that will rain down on your brain.
%By now, my friend, it should be quite plain:
%If you change these values ... you are insane.


%these markers identify the origins of coordinates. They get added to the
%values of coordinates. Depending on how the coordinate got come by, a
%certain marker will be added. In this way, the coordinates always carry
%that information

%If you wish, you can add others, to distinguish between more different
%kinds of cells. (A weird phrase meaning: you can add a new CLASS of cells)
%DO NOT USE 0.9 AS A MARKER!!! (up to 0.8 is ok)
%But you will have to add on the new marker to the right cells in the right 
%spot and ,after retrieving them, do something worthwhile with them.

%If you wish to mark cells with more then one property, I suggest you use
%markers placed at the second digit (0.01).

GoodCellMarker=0.1;
GoodCell=round(GoodCellMarker*10);

TempelCellMarker=0.2;
TempelCell=round(TempelCellMarker*10);

%NOOOO, DON'T DO IT!

NewCellMarker=0.3;
NewCell=round(NewCellMarker*10);

NewCelTempelMarker=0.4;
NewCelTempl=round(NewCelTempelMarker*10);

%DO NOT CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=First
while h+Increment<=Last
      h=h+Increment;
end
Last=h

howmanypics=Last-First+1;


levdiff_fi=abs(levnuc_fi-levback_fi)*LevelChanger;
levdiff_la=abs(levnuc_la-levback_la)*LevelChanger;
incr=(levdiff_la-levdiff_fi)/(howmanypics-1);
incrback=(levback_la-levback_fi)/(howmanypics-1);


i=0;
k=0;
a=0;
lostonedge=0;

picnew=[];
MPM=[];
PROPERTIES=[];
M=[];
coord=[];;

emptyM=zeros(1,4);

cd(path);
% mkdir(path,'bin_clust') ;
% 
% bin_clust_path=path,'bin_clust';



if ~Last > First
    disp('Fool, with these values (job#',num2str(projNum),') for the increment, the first and the last picture, not even one step is possible... Look again and choose more wisely');
    return
    
else    
    
    countingalong=0
	for looper=First:Increment:Last
        looper
        countingalong=countingalong+1
        % Format
		s=3; %s=length(num2str(no));
		strg=sprintf('%%.%dd',s);
		
		% Create numerical index
		indxStr=sprintf(strg,looper);
        
        name=[bodyfilename indxStr ext]; 
        
        %get the current picture
        picnew=imread(name);
           
        
     
		
        
	%     
	%     picnew=double(wholeheapimg(i).info);
        
        [img_h,img_w]=size(picnew);
        
        %[labeled,infos,avarea,ero]=body(picnew,round(cat(1,imagedata_struc.Centroid)),binaernuc);
	
        
        
        if looper == First
                     
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%% First picture is treated specially%%%%%%%% 
                      % Why? Because we have to define some filmspecific values and see if
                      % the first results are credible. Ok? here goes...
                    
                
                      %determin the the level of imextendmin by clicking (in the picture,
                      %which will pop up when running the programm)on a lot of nucloi
                      %(followed by ENTER) and several times on the background (followed by ENTER)
                      %This whole procedure has to be gone through twice, once for the
                      %first pic, once for the last. For every picture of the serie a
                      %value for lev will be interpolated
                      
%                  
%                       % Create numerical index
% 				      indxStr=sprintf(strg,Last);
%                         
%                       name=[bodyfilename indxStr ext]; 
%                         
%                       %get the current picture
%                       oldestpic=imread(name);
%            
%                       
%                      
%                      %level increment from pic to pic
% 			
% 	
%                      %incrback is the increment to interpolate the value of the
%                      %background, which we need for function halosfind
%                    
%                      clear levback_la;
%                      
                       
                     %now the results will be presented, so that, should they be bad,
                     %one can stop the programm immediatly. Last pic
%                      coordlast=[];
%                      regmax=[];
%                      [coordlast,regmax]=findnucloitrack(double(oldestpic),levdiff_la,minsizenuc,maxsizenuc);
%       %               figure, imshow(regmax);
%                      figure, imshow(oldestpic);
%                      hold on
%                      plot(coordlast(:,1),coordlast(:,2),'.')
%                      hold off
%                      
%                      clear coordlast;
%                      clear regmax;
% 			         clear levdiff_la;
			
             
%                      %first pic
%                      coordcol=[];
%                      regmax=[];
                     [coordcol,regmax]=findnucloitrack(picnew,levdiff_fi,minsizenuc,maxsizenuc);
%         %             figure, imshow(regmax);
%                      figure, imshow(picnew);
%                      hold on
%                      plot(coordcol(:,1),coordcol(:,2),'.')
%                      hold off
%                      
%                      
%                      altercoor=[];
                      [altercoor,logihalo]=halosfind(picnew,ErodeDiskSize,HaloLevel);
%                       
%                      
%                      %now follows a little something that will ensure a minimal
%                      %distance between two found cells
%                      namesnumbers=[0,0];
%                      if isempty(altercoor) == 0
%                                   for h=1:size(altercoor,1)
%                                                 uff=[];
%                                                 uff= min(sqrt((coordcol(:,1)-altercoor(h,1)).^2+(coordcol(:,2)-altercoor(h,2)).^2));
%                                                 
%                                                 if uff > (1.5*MinDistCellCell)
%                                                           namesnumbers(end+1,1)=altercoor(h,1);
%                                                           namesnumbers(end,2)=altercoor(h,2);
%                                                 end
%                                                 
%                                    end
%                                    
%                                    namesnumbers(1,:)=[];
%                                    
%                                    if isempty(namesnumbers) == 0
%                                                 coordcoll=cat(1,coordcol,namesnumbers);
%                                                 
%                                    else
%                                                 coordcoll=coordcol;
%                                                 
%                                    end
%                                  
%                                    
%                       else
%                                    coordcoll=coordcol;
%                       end
                      
                      coordcoll=handles.jobs(projNum).coordinatespicone;

                      PROPERTIES=[];
                      [PROPERTIES,ero,labeled]= body(picnew,coordcoll,regmax,logihalo,PlusMinus);
                      
                    
                      
                      emptyprop=zeros(1,size(PROPERTIES,2));
                      
                      clear regmax;
                      clear altercoor;
                      clear coordcol;
                      clear namesnumbers;
                      clear uff
                    
                   
                   
                      %manually fill in the missing and erase the wrong cells
                    
                    
			
                      %now we mark these coordinates as good coordinates, meaning their
                      %existence is credible. Futhermore cells marked as good cells will
                      %be treated prioritely, when it comes to allocating cells found in
                      %a new picture
                      coordcoll=coordcoll+GoodCellMarker;
                     
                
                      %%%%numberofnuc=size(coordcoll,2);
                      %%%%coordzuordn=round(cat(1,filmdata_struc(1).Centroid)
                    
                    
                      emptyM=zeros(1,4);
                     
                    
                    
                    
                    
                    
            
        %the following is what we do with pictures that aren't number one
        else
            
                     
        
                    
                      %adjust level via increment
                      levnucinterpol=levdiff_fi+(countingalong-1)*Increment*incr;
                      lebainterpol=levback_fi+(countingalong-1)*Increment*incrback;
                      
                      %find cells that look really dark and nasty
                      coordnuc=[];
                      [coordnuc,regmax]=findnucloitrack(picnew,levnucinterpol,minsizenuc,maxsizenuc);
                     
                      
                      %ensure a minimal distance between them
                      count=1;
                      while count < length(coordnuc)
                               paff=[];
                               paff= min(sqrt((coordnuc(count+1:end,1)-coordnuc(count,1)).^2 + (coordnuc(count+1:end,2)-coordnuc(count,2)).^2));
                               if paff < MinDistCellCell
                                       coordnuc(count,:)=[];
                                       count=count-1;
                               end
                               count=count+1;
                       end
                       clear count;    
                       clear paff;
                       
                       
                       %find cells that look like the third eye (round, big spots of
                       %pure light). We do this because the pictures are of poor
                       %quality and display huge halos around certain cells
                       [altercoor,logihalo]=halosfind(picnew,ErodeDiskSize,HaloLevel);
                       
                       
                       %ensure minimal distance between them
                       count=1;
                       while count < length(altercoor)
                               paffalt=[];
                               paffalt= min(sqrt((altercoor(count+1:end,1)-altercoor(count,1)).^2 + (altercoor(count+1:end,2)-altercoor(count,2)).^2));
                        
                               if paffalt < MinDistCellCell
                                       altercoor(count,:)=[];
                                       count=count-1;
                               end
                               count=count+1;
                       end
                       clear paffalt;
                       clear count;
                       
                       
                       %ensure minimal distance between dark, nasty-looking cells and
                       %third eyes (look above(findnucloitrack,halosfind))
                       namesnumbers=[0,0];
                       if isempty(altercoor) == 0
                            for h=1:size(altercoor,1)
                                    uff=[];
                                     
                                    uff= min(sqrt((coordnuc(:,1)-altercoor(h,1)).^2+(coordnuc(:,2)-altercoor(h,2)).^2));
                                    %note that here the minimal distance is larger than
                                    %between two cells found by the same routine,
                                    %because... ahhhmmm... just to make sure
                                    if uff  >  (1.5*MinDistCellCell)
                                          namesnumbers(end+1,1)=altercoor(h,1);
                                          namesnumbers(end,2)=altercoor(h,2);
                                    end
                             end
                             namesnumbers(1,:)=[];
                             
                             if isempty(namesnumbers) == 0
                                    coord1=cat(1,coordnuc,namesnumbers);
                             else
                                    coord1=coordnuc;
                             end
                                 
                       else
                              coord1=coordnuc;
                       end
                       clear altercoor;    
                       clear namesnumbers;   
                       clear coordnuc;   
                       clear uff;
                       
                       
                                
                               %calculate this stuff now, for conveniance (we need binary images of 
                               %the cells to run body)
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               prope=[];
                               ero=[];
                              
                               [prope,ero,labeled]= body(picnew,coordcoll,regmax,logihalo,PlusMinus);
                             
                               
                               
                               PROPERTIES(1:size(emptyprop,1),1:size(emptyprop,2),countingalong,1)=emptyprop;
                               %stuff this information into a stack
                               PROPERTIES(1:size(prope,1),1:size(prope,2),countingalong,1)=prope;
                               
                               clear prope;
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
                       
                       
                       
                       %mark coordinates as high quality coordinates. These cells
                       %are found in the usual way, so they are superior to cells
                       %found for instance by templates
                       coord1=round(coord1)+GoodCellMarker;
                          
                       
                          
                      
                       %we all know what the following is for - evading computer idiocity
                       coordcollflo=[];
                       coordextrainf=[];
                       rowcoord=[];
                       firstcoord=[];
                       secrowcoord=[];
                       secondcoord=[];
                       coord2=[];
                  
                  
                       %here we seperate the cells found already long ago from the ones found 
                       %only lately (we believe in the existence of the first more than in the 
                       %existence of the later)
                       coordcollflo=floor(coordcoll+0.1);
                       coordextrainf=round(10*(coordcoll-coordcollflo));
                       rowcoord= find(coordextrainf(:,1)==TempelCell | coordextrainf(:,1)==GoodCell);
                       
                       %old cells
                       firstcoord=coordcoll(rowcoord,:);
                       
                       clear rowcoord;
                       clear coordcollflo;
                       
                       
                       %new cells
                       secrowcoord= find(coordextrainf(:,1)==NewCell | coordextrainf(:,1)==NewCelTempl);
                       secondcoord=coordcoll(secrowcoord,:);
                       
                       clear coordextrainf;
                       clear secrowcoord;
                       
                       
                       %find out which cell is which over the gulf of different pictures
                       %(track between always TWO pics)
                      
                       %first we try to allocate the old cells to any set of coordinates
                       %found in the new picture
                       tmp1=[];
                       tmp1=fsmTrackTrackerBMTNN(firstcoord,coord1,radius);
                      
                       clear firstcoord;
                      % clear coord1;
                  
                       % if there are new cells, now is the time to try to allocate them
                       % to any set of coordinates not yet allocated to an old cell
                       if  size(find(tmp1(:,1)==0),1) >0 & size(secondcoord,1) >0
                                  identnew=find(tmp1(:,1)==0 & tmp1(:,2)==0);
                                  %in any case the sets coordinates not allocated to an old cell
                                  %are marked as a new cells, for they will either be allocated
                                  %to a new cell (found a few pictures ago), or actually become a 
                                  %new cell
                                  coord2(:,1)=round(tmp1(identnew,3))+NewCellMarker;
                                  coord2(:,2)=round(tmp1(identnew,4))+NewCellMarker;
                                  
                           
                                  %erase the unallocated set of coordinates from tmp1,
                                  %because they will turn up again in tmp2 no matter what
                                  %happens
                                  tmp1(identnew,:)=[];
                          
                                  clear identnew
                                  
                                  
                                  tmp2=[];
                                  tmp2=fsmTrackTrackerBMTNN(secondcoord,coord2,radius);
                                  
                                  %combine the two allocation matrices to one
                                  tmp=cat(1,tmp1,tmp2)
                                  
                                  clear tmp2;
                       else 
                                  identnew=find(tmp1(:,1)==0);
                                  tmp1(identnew,3:4)=round(tmp1(identnew,3:4))+NewCellMarker;
                                  
                                  clear identnew;
                                  
                                  if size(secondcoord,1) >0
                                          temzeros=zeros(size(secondcoord,1),4);
                                          temzeros(:,1:2)=secondcoord(:,1:2)
                                          tmp=cat(1,tmp1,temzeros);
                                          clear temzeros;
                                      
                                  else
                                          %since their are no new cells tmp2 is empty, so tmp1 is
                                          %complete
                                          tmp=tmp1    
                                  end
                       end
                       
                                   
                                   if length(find(tmp(:,3)~=0 & tmp(:,4)~=0)) ~= size(coord1,1)
                                       uuupsyouvelostone=1
                                   end
                                       
                                   
                                   if length(find(tmp(:,1)~=0 & tmp(:,2)~=0)) ~= size(coordcoll,1)
                                       BIGTROUBLE=1
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
                       %cells as "new cells tracked via template"
                       trouble=[];
                       trouble=find(tmp(:,3)==0);
                           
                           
                       if isempty (trouble) == 0
			
                                  for t= 1:size(trouble,1)
                                            tempcoord=[0,0];
                                                                
                                            tempcoord=[tmp(trouble(t),1),tmp(trouble(t),2)];
                                                   
				
                                             %first we check if the lost cell was located
                                             %near the edge (distance to edge smaller than
                                             %MinEdge). If so, we asume the cell has wandered
                                             %out of the picture and we stop following it
                                             if abs(tempcoord(1,1))<MinEdge | abs(tempcoord(1,1)-img_w)<MinEdge | abs(tempcoord(1,2))<MinEdge | abs(tempcoord(1,2)-img_h)<MinEdge
                                                         
                                                       %one never knows what kind of data
                                                       %somebody may wish to posses
                                                       lostonedge=lostonedge+1;
                                                     
                                             else
                                                      absmaxcorr=[];
                                                      newcoord=[];
                                                      
                                                      %this is the finding via templates
                                                      %routine. It will return the found
                                                      %coordinates and the value of the
                                                      %corralation
                                                      [newcoord,absmaxcorr]=templfindertrack(tempcoord,picold,picnew,lebainterpol,TempelCellMarker,NewCell,NewCelTempl,NewCelTempelMarker,percentbackground,sizetemple,box_size_img);
                                                      
                                                      %of course we only accept the correlation, if it of a minimal
                                                      %quality. NOTE WELL: that within the routine the background (everything 
                                                      %that doesn't belong to a cell) of the template and of the searcharea
                                                      %is replaced by random noise. So a background - background correlation
                                                      %will average zero
                            
                                                      if absmaxcorr > MinimalQualityCorr
                                                      
                                                                %now we calculate the minimal distance between the cell found by
                                                                %correlation and any other cell
                                                                dista=[];
                                                                rowofit=[];
                                                                
                                                                [dista,rowofit]= min(sqrt((tmp(:,3)-newcoord(1,1)).^2 + (tmp(:,4)-newcoord(1,2)).^2));
                                                                 
                                                             
                                                                if dista < round(MinDistCellCell/1.5)
                                                                    
                                                                        %the following is just extracting the information of
                                                                        %the coordinates found via template and the coordinates
                                                                        %the are very close to it
                                                                        tmpflo=[];
                                                                        extr=[];
                                                                        newflo=[];
                                                                        
                                                                        tmpflo=floor(tmp(rowofit,3)+0.1);
                                                                        extr=round(10*(tmp(rowofit,3)-tmpflo));
                                                                        
                                                                        clear tmpflo
                                                                        
                                                                        newflo=floor(newcoord(1,1)+0.1);
                                                                        extrnew=round(10*(newcoord(1,1)-newflo));
                                                                        
                                                                        clear newflo;
                                                                        
                                                                         
                                                                        %if the coordinates found via
                                                                        %template belong to an old cell,
                                                                        %we trust them, otherwise the
                                                                        %cell won't be propagated
                                                                        if extrnew == TempelCell
                                                                          
                                                                                %if the cell close by is a new cell, we allocate it's coordinates to the cell we try to track by template (which is 
                                                                                %an old cell - we've checked that above) and
                                                                                %erase the new cell
                                                                                 if  extr == NewCell
                                                                                        tmp(trouble(t),3)=round(tmp(rowofit,3))+GoodCellMarker;
                                                                                        tmp(trouble(t),4)=round(tmp(rowofit,4))+GoodCellMarker;
                                                                                        tmp(rowofit,3:4)=0;
                                                                                        
                                                                                         %if it is a new cell, itself propagated via
                                                                                         %templates, we erase it and allocate the
                                                                                         %coordinates we have (template) found
                                                                                         %to the old cell
                                                                                 elseif extr == NewCelTempl
                                                                                         tmp(rowofit,3:4)=0;
                                                                                         tmp(trouble(t),3)=newcoord(1,1);
                                                                                         tmp(trouble(t),4)=newcoord(1,2);
                                                                                 end
                                                                                 
                                                                          end
                                                                          clear extrnew;
                                                                          clear extr;
                                                                          
                                                                else
                                                                   
                                                                        %if the minimal distance to the
                                                                        %nearest cell is big enough,we allocate the
                                                                        %coordinates we have (template) found
                                                                        %to the old cell
                                                                        tmp(trouble(t),3)=newcoord(1,1);
                                                                        tmp(trouble(t),4)=newcoord(1,2);
                                                                       
                                                   
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
                      raus=find(ismember(tmp(:,1:4),[0 0 0 0],'rows')-1);
                      tmp=tmp(raus,:);
                      
                      clear raus;
                          
                          
                      
                     
                       M(1:size(emptyM,1),1:size(emptyM,2),countingalong-1,1)=emptyM;
                       %stuff this information into a stack
                       M(1:size(tmp,1),1:size(tmp,2),countingalong-1,1)=tmp;
                       
                          
                                   if length(find(tmp(:,1)~=0 & tmp(:,2)~=0)) ~= size(coordcoll,1)
                                       BIGTROUBLE=1
                                   end
                       
                       
                       clear tmp;
                       tmp=emptyM;
                       
                       
                       %a place for a break point, nothing more
                       if countingalong == 8
                           kooool=1
                       end
                       
                       
                       
                       
                       
                       
                       
                       
                       
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %this is a quite elaborate and (if I may say so) sophisticated
                       %part of the programm, which basically changes the tracks of
                       %certain cells, based on an analysis of the various tracks of
                       %cells over the last ? frames (defined in GUI_start). 
                   
                       if countingalong > howmanytimestepsslide+1
                           
                               templrow=[];
                               MPMslide=[];
                               MPMflo=[];
                               xtempl=[];
                               ytempl=[];
                               sizz=[];
                               
                               MPMslide=alteredfsmTrackLinker(M(:,:,countingalong-howmanytimestepsslide+1:countingalong-1));
                               %MPMslide is a matrix which gives the tracks of all cells over
                               %the last ? frames. It get's updated after every new frame.
                               
                               MPMflo=floor(MPMslide+0.1);
                               extrainf=round(10*(MPMslide-MPMflo));
                               %extrinf has the same size as MPMslide and every value reflects the extracted
                               %information of the coordinate located at the same spot in
                               %MPMslide. The information is already included in MPMslide. It
                               %is the digit after the point, which identifies the origin of a
                               %coordinate. Look above, marker
                              
                               clear MPMflo;
                               
                               
                               %After a few frames, new cells will be accepted as old cells, but only if, 
                               %after their initial finding, they weren't propagated by templates 
                               accept=find(extrainf(:,3)==NewCell);
                               
                               if isempty(accept) == 0
                                       for g=1:length(accept)
                                               if extrainf(accept(g),end)==NewCell;% & extrainf(accept(g),7)==NewCell & extrainf(accept(g),9)==NewCell %make sure a new cell wasn't propagated via templates
                                                       restore=[];
                                                       restore=find(M(:,3,countingalong-1)==MPMslide(accept(g),end-1) &  M(:,4,countingalong-1)==MPMslide(accept(g),end))  ;
                                                       M(restore,3:4,countingalong-1)=round(M(restore,3:4,countingalong-1))+GoodCellMarker;
                                                       MPMslide(accept(g),end-1)= M(restore,3,countingalong-1);
                                                       MPMslide(accept(g),end)= M(restore,4,countingalong-1);
                                                       
                                                       clear restore;
                                                        
                                               elseif  extrainf(accept(g),end)==NewCelTempl & extrainf(accept(g),end-2)==NewCelTempl
								                       annilhil=find(M(:,3,countingalong-1)==MPMslide(accept(g),end-1) &  M(:,4,countingalong-1)==MPMslide(accept(g),end) )  ;
								                       M(annilhil,3:4,countingalong-1)=0;
								                       MPMslide(accept(g),end-1)=0;
                                                       MPMslide(accept(g),end)=0;
                                                       
                                                       clear annilhil
                                                       
                                               end
                                       
                                       end
                               end
                               
                               clear accept;
                               
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            [shrow,shcol]= find(extrainf==0)
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            shadow=unique(shrow);
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            extrainf(shadow,:)=0;
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            clear shadow;
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            clear shrow;
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            clear shcol;
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            %a row containing zeros is no good, so we make sure they won't
	% % % % % % % % % % % % % % % % % % % % % % % % % %                            %be used for any of the following
                               
	
                               %here we search for old cells, that have been propagated via
                               %template search over the last five frames. Afterwards we will
                               %compare their tracks with tracks of newly found cells. If they
                               %corralate strongly, we replace the templatestuff (old cell) by the newly
                               %found cell
                             
                                templrow=find(ismember(extrainf(:,3:end),TempelCell*ones((howmanytimestepsslide-1)*2),'rows'));
                               
                               
                               if ~isempty(templrow)
                                      
                                  
                                     newrows= find(extrainf(:,3)==NewCell & extrainf(:,end)==NewCell);
                                     %newly found cells
                                   
                                     
                                     if ~isempty(newrows)
                                         
                                             if ~isempty(find(extrainf(newrows,3:end)))
                                                 BIGISHTROUBLE=1
                                             end
                                             
                                             xtempl=[];
                                             ytempl=[];
                                            
                                            
                                          
                                             %extract the x and y
                                             %coordinates of the cells         
                                          
                                             xtempl(:,:)=(MPMslide(templrow(:),3:2:end-1))'; 
                                             ytempl(:,:)=(MPMslide(templrow(:),4:2:end))';
						 
                                             sizz=size(xtempl,1);
                                                        
                                      
                                             %we want to correlate the displacementvectors
                                             %(displacement of a cell from one frame to the 
                                             %next), so we need to calculate it
                                             diffxx=diff(xtempl);
                                             diffyy=diff(ytempl);
                                            
                                             
                                             
                                             
                                             
                                             for g=1:length(newrows)
                                                                   
                                                            xnewone=[];
                                                            ynewone=[];
                                                            sizzle=[];
                                                            diffxxnew=[];
                                                            diffyynew=[];
                                                             
                                                            xnewone=(MPMslide(newrows(g),1:2:end-1))';
                                                            ynewone=(MPMslide(newrows(g),2:2:end))';
                                                            sizzle=length(xnewone);
                                                            
                                                            diffxxnew=diff(xnewone);
                                                            diffyynew=diff(ynewone);
                                                            
	% % %                                                                 for q=1:sizzle-1
	% % %                                                                             %vectors of the old location to the new
	% % %                                                                             %(frame to frame)
	% % %                                                                               diffxxnew(sizzle-q,1)=xnewone(sizzle+1-q)-xnewone(sizzle-q);
	% % %                                                                                diffyynew(sizzle-q,1)=ynewone(sizzle+1-q)-ynewone(sizzle-q);
	% % %                                                                 end
                                                             
                                                            
                          
                                                            %make sure we correlate two things of the same
                                                            %length
                                                            if sizz == sizzle
                                                                       correlas=[];
                                                                       nomore=[];
                                                                       for r=1:length(templrow)
                                                                         
                                                                             
                                                                                  lostfound=[];
                                                                                  tempfound=[];
                                                                                  writetemp=[];
                                                                                 
                                                                                  distance=min(sqrt((xtempl(:,r)-xnewone(:,1)).^2+(ytempl(:,r)-ynewone(:,1)).^2 ));
                                                                                 
                                                                                  
                                                                                  %now, if these points at some point are
                                                                                  %really close - exchange them
                                                                                  if distance < round(MinDistCellCell/2)
                                                                                      
                                                                                              %this is the substitution of the cell
                                                                                              %tracked by templates by the newly
                                                                                              %found cell. It's kind of tricky,
                                                                                              %because we have to change the values
                                                                                              %not in MPMslide, but in M. Since M
                                                                                              %is not sorted, we have to always
                                                                                              %find the right spot
                                                                                              IRGENDWIE_USTUUUUUSCHE=1
                                                                                              lostfound=[];
                                                                                              tempfound=[];
                                                                                                             
                                                                                              
                                                                                              %first we copy out the coordinates of
                                                                                              %the old cell, tracked by templates
                                                                                              lostfound=MPMslide(templrow(r),:);
                                                                                              %Now we erase the values from MPMslide, so that we won't by accident use them again
                                                                                              MPMslide(templrow(r),:)=0;
                                                                                           
                                                                                              %same procedure for the new cell we
                                                                                              %want to link to the old cell
                                                                                              tempfound=MPMslide(newrows(g),:)
                                                                                              MPMslide(newrows(g),:)=0;
                                                                                         
                                                                                              for h=1:((size(MPMslide,2)-2)/2)
                                                                                                 
                                                                                                         writetemp=[];
                                                                                                         gotemp=[];
                                                                                                         
                                                                                                         %ok,now: in lostfound is the track of the old 
                                                                                                         %cell.Tempfound is the track of the cell we wish to put at
                                                                                                         %it's place. Both are of the form:
                                                                                                         %[x1,y1,x2,y2,x3,y3,x4,y4,x5,y5]
                                                                                                         %So always two numbers correspond to the coordinates
                                                                                                         %in one frame(1,2,3...) mean different frames.
                                                                                                         %SOOOO: we find the location of the coordinates of the old
                                                                                                         %cell, by searching elements of lostfound in M. With every
                                                                                                         %increase in h we increase the index of lostfound by two
                                                                                                         %(because two values (x and y)per picture). M has three
                                                                                                         %dimensions. Every element of the third dimension is a
                                                                                                         %matrice(other two dimensions)of the form 
                                                                                                         %[x1,y1,x2,y2],thus giving the connection
                                                                                                         %between two frames, BUUUUTTTT: the next element 
                                                                                                         %will be [x2,y2,x3,y3]!!!!!! This is very important!!!!!!. 
                                                                                                         %So with every frame (h) we make one step into the third dimesion in M.
                                                                                                         %with each set of x and y of lostfound we locate the rigth
                                                                                                         %row in M. 
                                                                                                         writetemp=find(M(:,1,countingalong-5+h)== lostfound(2*h-1)  &  M(:,2,countingalong-5+h)== lostfound(2*h));
                                                                                                         
                                                                                                         %since we overwrite the old
                                                                                                         %cell with the new, we have to
                                                                                                         %eliminate the new, in order
                                                                                                         %not to have it twice
                                                                                                          gotemp=find(M(:,1,countingalong-5+h)== tempfound(2*h-1)  &  M(:,2,countingalong-5+h)== tempfound(2*h));
                                                                                                          
                                                                                                          %now we juggle the information into the right place.
                                                                                                          %Since M is
                                                                                                          %[x1,y1,x2,y2] , next layer [x2,y2,x3,y3]
                                                                                                          %and tempfound is
                                                                                                          %[x1,y1,x2,y2,x3,y3,x4,y4,x5,y5]
                                                                                                          %therefor:  
                                                                                                          %M(writetemp,:,countingalong-4) replaced by tempfound(1:4)
                                                                                                          %and
                                                                                                          %M(writetemp,:,countingalong-3) replaced by tempfound(3:6) 
                                                                                                          %and so on.
                                                                                                          M(writetemp,:,countingalong-5+h)=tempfound(2*h-1:2*h+2);
                                                                                              
                                                                                                          %erase redundant info
                                                                                                          M(gotemp,:,countingalong-5+h)=0;
                                                                                                         
                                                                                                         
                                                                                               end
                                                                                               clear writetemp;
                                                                                               clear gotemp;
                                                                                               clear lostfound;
                                                                                               clear tempfound;
                                                                                     
                                                                                     
                                                                                               nomore=1;
                                                                                               break
                                                                                           
                                                                                    %here too, we have to stay within a maximal distance
                                                                                    elseif distance < DistanceCorrel 
                                                                                              
                                                                                                   %now we correlate the displacement vectors of the old cell with
                                                                                                   %the ones of the new cell
                                                                                                
                                                                                                   correlx=xcorr(diffxx(:,r), diffxxnew,'coeff');
                                                                                                   correly=xcorr(diffyy(:,r),diffyynew,'coeff');
                                                                                             
                                                                                                   oops=round(length(correlx)/2);
                                                                                             
                                                                                                   correlas(r)=correlx(oops)+correly(oops);
                                                                                                   
                                                                                    end
                                                                         end
                                                                     
                                                                         %nomore says if this new cell was already asigned to an old cell by distance criteria
                                                                         %if this isn't the case, it's time to see, if there is an old cell (nearby as checked above)
                                                                         %that correlates nicely
                                                                         if isempty(nomore)
                                                                              
                                                                                     if ~isempty(correlas);
                                                                                                 
                                                                                                 important=[];
                                                                                                 [crap,important]=max(correlas) ;
                                                                                                 
                                                                                                 if crap > MinTrackCorr
                                                                                                     
                                                                                                             %this is exactly the same procedure as above. I know it would be clever to make a function out of it
                                                                                                             %but there would be so many BIG matrices involved I decided to have it in here twice. Ah well. Anyway look above
                                                                                                             %(look for IRGENDWIE_USTUUUSCHE)
                                                                                                             IRGENDWIE_USTUUUSCHE=1
                                                                                                             lostfound=[];
                                                                                                             tempfound=[];
                                                                                                             
                                                                                                             lostfound=MPMslide(templrow(important),:);
                                                                                                             MPMslide(templrow(important),:)=0;
                                                                                                               
                                                                                                             tempfound=MPMslide(newrows(g),:);
                                                                                                             MPMslide(newrows(g),:)=0;
                                                                                                             
                                                                                                             for h=1:(size(MPMslide,2)-2)/2
                                                                                                                      writetemp=[];
                                                                                                                      gotemp=[];
                                                                                                                 
                                                                                                                      writetemp=find(M(:,1,countingalong-5+h)== lostfound(2*h-1)  &  M(:,2,countingalong-5+h)== lostfound(2*h));
                                                                                                                      gotemp=find(M(:,1,countingalong-5+h)== tempfound(2*h-1)  &  M(:,2,countingalong-5+h)== tempfound(2*h));
                                                                                                                      
                                                                                                                      
                                                                                                                      M(writetemp,:,countingalong-5+h)=(tempfound(2*h-1:2*h+2));
                                                                                                                      
                                                                                                                  
                                                                                                                      M(gotemp,:,countingalong-5+h)=0;
                                                                                                                      
                                                                                                             end
                                                                                                             clear writetemp;
                                                                                                             clear gotemp;
                                                                                                             clear lostfound;
                                                                                                             clear tempfound;
                                                                                                 end
                                                                                                 clear important;
                                                                                      end    
                                                                         end
                                                                         clear correlas;
                                                                         clear nomore;
                                                            else
                                                                          ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRROOOOOOOOOOOOOOOOOOR=1;
                                                                          %if you've landed here by executing the programm your in TROUBLE. But it's basically impossible
                                                                          %so you haven't - and you needn't worry. UNLESS you've been meddling with the programm. Tut tut.
                                                                          %But actually it isn't THAT bad. The displacement vector of old cell and new cell (diffxx,diffxxnew...)
                                                                          %aren't of the same length, so check what's wrong.
                                                                          
                                                            end
                                                            clear xnewone;
                                                            clear ynewone;
                                                            clear sizzle;
                                                            clear diffxxnew;
                                                            clear diffyynew;
                                                              
                                              end
                                        end
                               end
                               clear templrow;
                               clear MPMslide;
                               clear xtempl;
                               clear ytempl;
                               clear sizz;
                       end
                        
                      
                    
                       clear coord1;
                       
                       
                       
                       coordcoll=[];
                       realcoord=[];
                       realcoord=find(M(:,3,countingalong-1) | M(:,4,countingalong-1));
                      
                       %coordcoll=coord1;
                    
                    
                       coordcoll(:,1)=M(realcoord,3,countingalong-1)   ;%             tmp(realcoord,3);
                       coordcoll(:,2)=M(realcoord,4,countingalong-1)  ;
                       clear realcoord;
             
                       
                       
          end
          picold=picnew;
          
          
          
          
	% Format
	s=3; %s=length(num2str(no));
	strg=sprintf('%%.%dd',s);
	
	% Create numerical index
	indxStr=sprintf(strg,looper);
	
	figure, imshow(ero)
	eval(['print -djpeg100 ',pathbring,filesep,'clusters_binary',indxStr,'.jpg']);
          
	close
          
	end
      
	%it's all over now. All we have to do is change the format of the
	%information we have painstakingly gathered over the last ??? lines
	%MPM=fsmTrackLinker(M);
             
    %we save all the data we found to some directory. It looks awefully
    %puny after so many lines of code.
    cd(pathbring)
    jobvalues=handles.jobs(projNum);
    
    save('M', 'M')
    save ('PROPERTIES', 'PROPERTIES')
    save (['lostonedge',date], 'lostonedge')
    save ('jobvalues','jobvalues')
 
    
end
         