function ptMovieMaker (hObject) 
% ptMovieMaker makes movies from the information gathered in the analysis
%
% SYNOPSIS       ptMovieMaker (hObject)
%
% INPUT          hObject : handle to an object within PolyTrack_PP
%
% OUTPUT         no function output: movies are saved to disk          
%
% DEPENDENCIES   ptMovieMaker uses { nothing }
%                                  
%                ptMovieMaker is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04
handles = guidata (hObject);

% Determine the size of the MPM matrix
[numRows, numCols] = size (handles.MPM);

% And get the image name list
imageNameList = handles.jobvalues.imagenameslist;

% How many frames does the user wants to see the tracks of?
dragTailLength = handles.postpro.dragtail;

% What is the name of the movie
dragTailFileName = handles.postpro.dragtailfile;

% Starting frame for the movie
startFrame = handles.postpro.moviefirstimg;

% We need prior images for the dragTail movie
if startFrame < dragTailLength + 2
    startFrame = dragTailLength + 2;
end

% Last frame for the movie
lastFrame = handles.postpro.movielastimg;

% Make sure it doesn't go out of range
if lastFrame > handles.jobvalues.lastimage
   lastFrame = handles.jobvalues.lastimage;
elseif lastFrame < startFrame
   h = errordlg('First frame, last frame and dragtail length are not compatible! Please change these.');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end;

% Get image and save directory and go to save dir
imageDirectory = handles.jobvalues.imagedirectory;
savePath = handles.postpro.saveallpath;
cd (savePath);

% Initialize the movie
makeQTMovie ('start', dragTailFileName);

% Initialize counter
frameCounter = startFrame - 1;

% Let the user know we are starting
fprintf (1, 'ptMovieMaker: starting to generate movie frames %d to %d...\n', startFrame, lastFrame);

% Start doing the actual work to create the movie
for movieStep = startFrame : lastFrame
    
   % counter
   frameCounter = frameCounter + 1;
   fprintf (1, 'ptMovieMaker: Creating movie frame # %d ...\n', frameCounter); 
    
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
   
   cd (imageDirectory);

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
   colorMap = jet (dragTailLength + 1);

   % Initialize the dragtail counter
   colorCount = 0;
        
   % Loop through the previous pictures to generate the dragtails
   for iCount = (2 * (movieStep - dragTailLength)) : 2 : (2 * movieStep)
      colorCount = colorCount + 1;
      vec = handles.MPM (selectedCells, iCount-3 : iCount);
      [rows, cols] = find (vec == 0);
      rows = unique (rows);
      vec (rows,:) = 0;

      for hCount = 1 : size (vec,1)
         if vec (hCount,1) ~= 0
	        ph = [];
	        ph = plot (vec (hCount, 1:2:3), vec (hCount, 2:2:4));
	        set (ph, 'Color', colorMap (colorCount,:));
	        clear ph;
         end   % if vec
      end   % for hCount
   end   % for iCount
        
   % Plot the points that actually belong to the current picture
   plot (vec (:,3), vec (:,4), 'r.');
   hold off;

   cd (savePath);
  
   % Add the current figure to the movie
   makeQTMovie ('addaxes', gca);
   
   % Close the figure
   close;
end   % for movieStep

% finalize the movie and write it to disk
makeQTMovie ('finish');

% Let the user know we have finished
fprintf (1, 'ptMovieMaker: finished generating movie.\n');
