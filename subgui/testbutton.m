function testbutton(hObject)

%determin the the level of imextendmin by clicking (in the picture,
%which will pop up when running the test)several times on the background
%(followed by ENTER), on a lot of nuclei (followed by ENTER) and on a lot of
%halos (followed by ENTER).
%This whole procedure has to be gone through twice, once for the
%first pic, once for the last. For every picture of the serie a
%value for lev will be interpolated

     

handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');




ImageDirectory = handles.jobs(projNum).imagedirectory;



First = handles.jobs(projNum).firstimage;
Last = handles.jobs(projNum).lastimage;


levnuc_fi = handles.jobs(projNum).fi_nucleus;
levback_fi = handles.jobs(projNum).fi_background;
levhalo_fi = handles.jobs(projNum).fi_halolevel;

ImageNamesList = handles.jobs(projNum).imagenameslist

% levnuc_la = handles.jobs(projNum).la_nucleus;
% levback_la = handles.jobs(projNum).la_background;
levdiff_fi = abs(levnuc_fi-levback_fi)*handles.jobs(projNum).leveladjust;


ErodeDiskSize = 10;
%how much the blobs found in halosfind shall be eroded. This is an indirect
%size criteria for halosfind. Increase - minimal size of halos will be
%increased, decrease - ... decreased


%minimal/maximal size of the black spot in the cells
minsizenuc = handles.jobs(projNum).minsize;
maxsizenuc = handles.jobs(projNum).maxsize;


MinDistCellCell = handles.jobs(projNum).minsdist;
%an educated guess of the minimal distance between neighbouring cells
%(better to big than to small, for this value is used for the static
%search. Searches with templates are more tolerant.)


cd(ImageDirectory)



if ~Last > First
    disp('Fool, with these values (job#',num2str(projNum),') for the increment, the first and the last picture, not even one step is possible... Look again and choose more wisely');
    return
    
else    

%     
%         % Format
% 		s=3; %s=length(num2str(no));
% 		strg = sprintf('%%.%dd',s);
% 		
% 		% Create numerical index
% 		indxStr = sprintf(strg,First);
%         
%         name = [bodyfilename indxStr ext]; 
        
        name = char(ImageNamesList(First))
        
        %get the current picture
        firstFrame = imreadnd2(name,0,handles.jobs(projNum).intensityMax);
           
        

        [img_h,img_w] = size(firstFrame);
       
        [seg_img, obj_val] = imClusterSeg(firstFrame, 1, 'method','kmeans','k_cluster',3,'mu0', [levnuc_fi;levback_fi;levhalo_fi]);

        coordNuc = [];
        regmax = [];
        [coordNuc,regmax] = findnucloitrack(seg_img,levdiff_fi,minsizenuc,maxsizenuc);
    
       
        HaloLevel = (levhalo_fi-levback_fi)*2/3+levback_fi;
        coordHalo = [];
        [coordHalo,logihalo] = halosfind(seg_img,ErodeDiskSize,HaloLevel);
          
         
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