function movieMaker (hObject) 
% movieMaker makes movies from the information gathered in the analysis
%
% SYNOPSIS       movieMaker (hObject)
%
% INPUT          hObject : handle to an object within PolyTrack_PP
%
% OUTPUT         no function output: movies are saved to disk          
%
% DEPENDENCIES   movieMaker uses { nothing }
%                                  
%                movieMaker is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04
handles = guidata (hObject);

% Determine the size of the MPM matrix
[numRows, numCols] = size (handles.MPM);

% And get the image name list
imageNameList = handles.jobvalues.imagenameslist;

% How many frames does the user wants to see the tracks of?
dragTailLength = handles.postpro.dragtail;

% Starting frame relativ to first frame analysed with polytrack
startFrame = round ((handles.postpro.moviefirstimg - handles.jobvalues.firstimage) / ...
             handles.jobvalues.increment) + 1;

% We need prior images for the dragTail movie
if startFrame < dragTailLength + 2
    startFrame = dragTailLength + 2;
end


% Last frame relativ to first frame analysed with polytrack
lastFrame = floor ((handles.postpro.movielastimg - handles.jobvalues.firstimage) / ...
            handles.jobvalues.increment + 0.00001) + 1;

% Make sure it doesn't go out of range
if lastFrame > handles.jobvalues.lastimage
    lastFrame = handles.jobvalues.lastimage;
end

% Get image and save directory and go to save dir
imagedirectory = handles.jobvalues.imagedirectory;
saveallpath = handles.postpro.saveallpath;
cd (saveallpath);

% Initialize the movie
makeQTMovie ('start','trackmov.mov');

% Start doing the actual work to create the movie
for movieStep = startFrame : lastFrame
    
   % Use only the cells chosen by the user
   if ~isempty (handles.selectedcells)
      selectedCells = zeros (size (handles.selectedcells, 1), 2);
      selectedCells(:,:) = handles.MPM (handles.selectedcells, (2 * movieStep - 1):(2 * movieStep));
   else
      selectedCells = zeros (size (handles.MPM, 1), 2);
      selectedCells(:,:) = handles.MPM (:, (2 * movieStep - 1):(2 * movieStep));
   end

   % selectedCells defines which cells will be taken into account. This is
   % either defined by the user (PolyTrack_PP) or it will just be all cells
   % within the current picture
     
   selectedCells = find (selectedCells(:,1) & selectedCells(:,2)); 
   
   cd (imagedirectory);

   name = char (imageNameList (movieStep)); 
   nowImgH = imreadnd2 (name, 0, handles.jobvalues.intensityMax);

   [rows,cols] = find (handles.MPM);
		
   % if the user specified a size for the movie, that's the size we are going to use
   if ~isempty (handles.postpro.figureSize)
      figure ('Position', handles.postpro.figureSize);
      imshow (nowImgH, []);
   else
      figure, imshow (nowImgH, []);
   end

   hold on;
   % One colour per time step (dragTailLength tells you how many
   % timesteps there are
   cmap = jet (dragTailLength + 1);

   counter=0;
        
   % Loop through the previous pictures (for the tails)
   for i = (2 * (movieStep - dragTailLength)) : 2 : (2 * movieStep)
      counter = counter + 1;
      vec = handles.MPM (selectedCells, i-3:i);
      [rows,cols] = find (vec==0);
      rows = unique (rows);
      vec(rows,:) = 0;

      for h = 1 : size(vec,1)
         if vec(h,1) ~= 0
	    ph = [];
	    ph = plot (vec (h, 1:2:3), vec (h, 2:2:4));
	    set (ph, 'Color', cmap (counter,:));
	    clear ph;
         end
      end
   end
        
   % Plot the points that actually belong to the current picture
   plot (vec (:,3), vec (:,4), 'r.');
   hold off;

   cd (saveallpath);
  
   % Add the current figure to the movie
   makeQTMovie ('addaxes', gca);
   close;
end

% finalize the movie and write it to disk
makeQTMovie ('finish');
