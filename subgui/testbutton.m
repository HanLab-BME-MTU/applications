function testbutton(hObject)



handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');




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





h=First
while h+Increment<=Last
      h=h+Increment;
end
Last=h




levdiff_fi=abs(levnuc_fi-levback_fi)*LevelChanger;



picnew=[];

coord=[];;




if ~Last > First
    disp('Fool, with these values (job#',num2str(projNum),') for the increment, the first and the last picture, not even one step is possible... Look again and choose more wisely');
    return
    
else    

    
        % Format
		s=3; %s=length(num2str(no));
		strg=sprintf('%%.%dd',s);
		
		% Create numerical index
		indxStr=sprintf(strg,First);
        
        name=[bodyfilename indxStr ext]; 
        
        %get the current picture
        picnew=imread(name);
           
        
     
		
        
	%     
	%     picnew=double(wholeheapimg(i).info);
        
        [img_h,img_w]=size(picnew);
       
                 
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%% First picture is treated specially%%%%%%%% 
                      % Why? Because we have to define some filmspecific values and see if
                      % the first results are credible. Ok? here goes...
                    
                
                      %determin the the level of imextendmin by clicking (in the picture,
                      %which will pop up when running the programm)on a lot of nucloi
                      %(followed by ENTER) and several times on the background (followed by ENTER)
                      %This whole procedure has to be gone through twice, once for the
                      %first pic, once for the last. For every picture of the serie a
                      %value for lev will be interpolated
                      
                 
             
                      
	
                       
              
          
                     coordcol=[];
                     regmax=[];
                     [coordcol,regmax]=findnucloitrack(picnew,levdiff_fi,minsizenuc,maxsizenuc);
                
                   
                     
                     altercoor=[];
                     [altercoor,logihalo]=halosfind(picnew,ErodeDiskSize,HaloLevel);
                      
                     
                     %now follows a little something that will ensure a minimal
                     %distance between two found cells
                     namesnumbers=[0,0];
                     if isempty(altercoor) == 0
                                  for h=1:size(altercoor,1)
                                                uff=[];
                                                uff= min(sqrt((coordcol(:,1)-altercoor(h,1)).^2+(coordcol(:,2)-altercoor(h,2)).^2));
                                                
                                                if uff > (1.5*MinDistCellCell)
                                                          namesnumbers(end+1,1)=altercoor(h,1);
                                                          namesnumbers(end,2)=altercoor(h,2);
                                                end
                                                
                                   end
                                   
                                   namesnumbers(1,:)=[];
                                   
                                   if isempty(namesnumbers) == 0
                                                coordcoll=cat(1,coordcol,namesnumbers);
                                                
                                   else
                                                coordcoll=coordcol;
                                                
                                   end
                                 
                                   
                      else
                                   coordcoll=coordcol;
                      end
                      
                     
                    
                      
        
                   
                   
                      %manually fill in the missing and erase the wrong cells
                      coordcoll=takenkick(picnew,coordcoll);
                    
			
              coordcoll=round(coordcoll);
                     
               
                    
                      
                     
                    
                    








handles.jobs(projNum).coordinatespicone= coordcoll;


% Update handles structure
guidata(hObject, handles);

end