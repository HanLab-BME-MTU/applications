function testbutton(hObject)
% testbutton finds coordinates in an image
%
% SYNOPSIS       trackCells(hObject,projNum)
%
% INPUT          hObject : a handle to the Gui which called the function
%                
% OUTPUT         list of coordinates manually checked and edited by user.
%                Written into handles
%                             
%
% DEPENDENCIES   testbutton uses {imClusterSeg
%								 trackLinker
%								 checkMinimalCellCell
%								 findNucloiTrack
%								 takenkick
%								 ptFindHalos}
%                                  
%                testbutton is used by { PolyTrack }
%
% REMARK         testbutton fetches directly in GUI PolyTrack:
% 					projNum : currently treated job
%                   handles : structure with information used within GUI
% 				  from handles.jobs(projNum):
% 					imagedirectory : where are the images 
% 					imagename : what are the images called
% 					imagenameslist : list of images within imagedirectory with imagename
% 					firstimage : which images shall we start with (refers to imagenameslist)
% 					lastimage : which images shall be the last one (refers to imagenameslist)
% 					intensityMax : highest value image can have (calc from bitdepth)
% 					fi_nucleus : nucloi intensity first image
% 					fi_background : background intensity first image
% 					fi_halolevel : halo intensity first image
% 					leveladjust : factor to adjust intensity difference nucloi/ background
% 					minsize : minimal size of nucloi 
% 					maxsize : maximal size of nucloi
% 					minsdist : minimal distance between two cells
% 					minmaxthresh : onoff - should minima and segmentation be used 
% 					clustering : onoff - should clustering be used
%
% Colin Glass, Feb 04




%determin the the level of imextendmin by clicking (in the picture,
%which will pop up when running the test)several times on the background
%(followed by ENTER), on a lot of nuclei (followed by ENTER) and on a lot of
%halos (followed by ENTER).
%This whole procedure has to be gone through twice, once for the
%first pic, once for the last. For every picture of the serie a
%value for lev will be interpolated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fetch information

handles = guidata(hObject);
projNum = get(handles.GUI_st_job_lb,'Value');

% First tell the user we're busy initializing
fprintf (1, 'Initializing phase has started for job number %d...\n', projNum);

ImageDirectory = handles.jobs(projNum).imagedirectory;
First = handles.jobs(projNum).firstimage;
Last = handles.jobs(projNum).lastimage;
levnuc_fi = handles.jobs(projNum).fi_nucleus;
levback_fi = handles.jobs(projNum).fi_background;
levhalo_fi = handles.jobs(projNum).fi_halolevel;
ImageNamesList = handles.jobs(projNum).imagenameslist;
leveladjust = handles.jobs(projNum).leveladjust;

%minimal/maximal size of the black spot in the cells
MinSizeNuc = handles.jobs(projNum).minsize;
MaxSizeNuc = handles.jobs(projNum).maxsize;

MinDistCellCell = handles.jobs(projNum).minsdist;
%an educated guess of the minimal distance between neighbouring cells
%(better to big than to small, for this value is used for the static
%search. Searches with templates are more tolerant.)

segmentation = handles.jobs(projNum).minmaxthresh;
clustering = handles.jobs(projNum).clustering;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


levdiff_fi = abs(levnuc_fi-levback_fi)* leveladjust;


ErodeDiskSize = round((sqrt(MinSizeNuc))/2) ;
%how much the blobs found in ptFindHalos shall be eroded. This is an indirect
%size criteria for ptFindHalos. Increase - minimal size of halos will be
%increased, decrease - ... decreased


cd(ImageDirectory)


if ~Last > First
    fprintf(1, 'Error: the last image # is before the first image # in job number %d\n', projNum);
    return
else    

        name = char(ImageNamesList(First));
        
        %get the current picture
        firstFrame = imreadnd2(name,0,handles.jobs(projNum).intensityMax);
           
        

        [img_h,img_w] = size(firstFrame);
       
         HaloLevel = (levhalo_fi-levback_fi)*2/3+levback_fi;
        
		
		if clustering
         
%               background = imopen(firstFrame,strel('disk',80));
%               imgMinusBack = imsubtract(firstFrame,background); 
%               smallest = min(min(imgMinusBack));
%               if sign(smallest)== -1
%                   imgMinusBack = imgMinusBack + abs(smallest)
%               end
              
              [seg_img, dummy,mu0] = imClusterSeg(firstFrame, 1, 'method','kmeans','k_cluster',3,'mu0', [levnuc_fi;levback_fi;levhalo_fi]);
			
              %find cells that look really dark and nasty
              [coordNuc,regmax] = findNucloiTrack(seg_img,levdiff_fi,MinSizeNuc,MaxSizeNuc,1);
             
              %find cells that look like the third eye (round, big spots of
              %pure light). We do this because the pictures are of poor
              %quality and display huge halos around certain cells
              [coordHalo,logihalo] = ptFindHalos(seg_img,ErodeDiskSize,HaloLevel,1);
		
		elseif segmentation
           
              [coordNuc,regmax] = findNucloiTrack(firstFrame,levdiff_fi,MinSizeNuc,MaxSizeNuc,2);
              [coordHalo,logihalo] = ptFindHalos(firstFrame,ErodeDiskSize,HaloLevel,2);
			
		
		else
              disp('You have to choose one of the methods, clustering or segmentation')
              return
		end
		
         
        %now follows a little something that will ensure a minimal
        %distance between two found cells
             
        newCoord =  checkMinimalCellCell(coordNuc,coordHalo,MinDistCellCell);           
 

        %manually fill in the missing and erase the wrong cells
        newCoord = takenkick(firstFrame,newCoord);
       
        newCoord = round(newCoord);
         
		
		handles.jobs(projNum).coordinatespicone =  newCoord;
		
		
		% Update handles structure
		guidata(hObject, handles);
		
		%%%%%%%%save altered values to disk%%%%%%%%%%%%
		cd(handles.jobs(projNum).savedirectory)
		jobvalues = handles.jobs(projNum);
		save ('jobvalues','jobvalues')
		clear jobvalues
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
