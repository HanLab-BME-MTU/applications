function ptSetCellValues (hObject,objectChoice)
% ptSetCellValues let's the user specify values of images
%
% SYNOPSIS       ptSetCellValues (hObject,objectChoice)
%
% INPUT          hObject : handle to an object of GUI (PolyTrack)
%                objectChoice : what kind of value interestes the user 
%                                1 - minimal size nucloi
%                                2 - maximal size nucloi
%                                3 - minimal distance between two nuclei
%
% OUTPUT         none (results get written directly into handles)
%
% DEPENDENCIES   ptSetCellValues uses { nothing }
%                                  
%                ptSetCellValues is used by { PolyTrack }
%
% REMARK         ptSetCellValues fetches directly in GUI PolyTrack:
% 		 jobNumber : current job
%                handles : structure with information used within GUI
% 		   from  handles.jobs (jobNumber):
% 		      imagedirectory : where are the images 
% 		      imagename : what are the images called
% 		      imagenameslist : list of images within imagedirectory with imagename
% 		      firstimage : which images shall we start with (refers to imagenameslist)
% 		      lastimage : which images shall be the last one (refers to imagenameslist)
% 		      intensityMax : highest value image can have (calc from bitdepth)
% 					
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source
% Andre Kerstens        May 04          Renamed to ptSetCellValues.m
% Andre Kerstens        Jul 04          Added parameters to ptGetProcessedImage
% Andre Kerstens        Dec 04          Bugfix in imshow (added []) and
%                                       added try/catch loops

% Get the handles structure from the GUI and determine current job number
handles = guidata (hObject);
jobNumber = get (handles.GUI_st_job_lb,'Value');

change=0;

imageDirectory = handles.jobs(jobNumber).imagedirectory;
firstImage    = handles.jobs(jobNumber).firstimage;
imageNameList = handles.jobs(jobNumber).imagenameslist;
intensityMax   = handles.jobs(jobNumber).intensityMax;

% Fetch the first image from disk
cd (imageDirectory);
fileName = char (imageNameList (firstImage));
tempFirstImg = imreadnd2 (fileName, 0, intensityMax);
[firstImg, backgroundLevel, edgeImage] = ptGetProcessedImage (tempFirstImg, intensityMax, 15, 5);

% Get the image size
[img_h,img_w]=size(firstImg);

% Depending on what we have to do go through one of the code parts below

if objectChoice == 1               % Minimal size of nuclei (1)
   % Show figure with first image
   figure, imshow (firstImg,[]);
   title('Draw a polygon around the smallest nucleus and press ENTER when finished.');

   % Get the polygon coordinates and store in polyImage (this is a binary image 
   % containing the selected region of firstImg)
   try
      polyImage = roipoly;
   catch
      % Close figure
      close;      
      return; 
   end

   % Close figure
   close;

   if length (find (polyImage)) > 0
      % Store the length of the selected region (which is the nucleus diameter) in the handles struct
      handles.jobs(jobNumber).minsize = length (find (polyImage));
      change = 1;
   end

elseif objectChoice == 2            % Maximal size of nuclei (2)
   % Show figure with first image
   figure, imshow (firstImg,[]);
   title('Draw a polygon around the biggest nucleus and press ENTER when finished.');

   % Get the polygon coordinates and store in polyImage (this is a binary image 
   % containing the selected region of firstImg)
   try
      polyImage = roipoly;
   catch
      % Close figure
      close;      
      return;       
   end

    % Close figure
   close;

   if length (find (polyImage)) > 0
      % Store the length of the selected region (which is the nucleus diameter) in the handles struct
      handles.jobs(jobNumber).maxsize = length (find (polyImage));
      change = 1;
   end

elseif objectChoice == 3            % Minimal distance between two nuclei (3)
   % Show figure with first image
   figure, imshow (firstImg,[]);
   title('Click on the centerpoints of the two nuclei closest to each other and press ENTER.');

   % Get the coordinates of the pixels clicked by the user
   try
      [x,y] = getpts;
   catch
      % Close figure
      close;      
      return;
   end

   % Close figure
   close;
     
   if length (x) > 1
      % If more than 2 points clicked, take the first two
      handles.jobs(jobNumber).minsdist = round (sqrt ((x(1) - x(2))^2 + (y(1) - y(2))^2));
      change = 1;
   else
      fprintf (1, 'Warning: Only one nuclei selected (two are needed). Returning...\n');
      return;
   end
else
   fprintf (1, 'Warning: Invalid choice for object. Returning...\n');
   return;
end


%change equals one, if anything has been specified by the user
if change==1
	% Update handles structure
	guidata(hObject, handles);
	
	%%%%%%%%save altered values to disk%%%%%%%%%%%%
	cd(handles.jobs(jobNumber).savedirectory)
	jobvalues=handles.jobs(jobNumber);
	save ('jobvalues','jobvalues')
	clear jobvalues
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%now we update the value in PolyTrack
	ptFillFields (handles,handles.jobs(jobNumber));
end
