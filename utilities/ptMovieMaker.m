function ptMovieMaker (ptPostpro, MPM) 
% ptMovieMaker makes movies from the information gathered in the analysis
%
% SYNOPSIS       ptMovieMaker (ptPostpro, MPM)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI
%                MPM       : matrix containing the cell tracks
%
% OUTPUT         no function output: movies are saved to file specified by the user         
%
% DEPENDENCIES   ptMovieMaker uses { nothing }
%                                  
%                ptMovieMaker is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Jun 04          Cleaned up source and renamed file.
%                                       Also made it independent of gui handles

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
movieStartFrame = ptPostpro.moviefirstimg;
movieEndFrame = ptPostpro.movielastimg;
numberOfFrames = ceil((movieEndFrame - movieStartFrame) / increment) + 1;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;
imageDirectory = ptPostpro.imagepath;
imageNameList = ptPostpro.imagenameslist;
intensityMax = ptPostpro.intensitymax;
figureSize = ptPostpro.figureSize;

% How many frames does the user wants to see the tracks of?
dragTailLength = ptPostpro.dragtail;

% What is the name of the movie
dragTailFileName = ptPostpro.dragtailfile;

% We need prior images for the dragTail movie
if ceil ((movieStartFrame - startFrame + 1) / increment) < dragTailLength + 1
    movieStartFrame = ((dragTailLength + 1) * increment) + 1;
end
 
% Go to the directory where the image will be saved
cd (savePath);

% Initialize the movie
makeQTMovie ('start', dragTailFileName);
%mov = avifile ('dragtail.avi')

% Initialize MPM counter
MPMCount = ceil ((movieStartFrame - startFrame) / increment);

% Let the user know we are starting
fprintf (1, '\nptMovieMaker: Generating movie frames: ');

% Start doing the actual work to create the movie
for movieStep = movieStartFrame : increment : movieEndFrame
    
   % Let the user know where we are
   fprintf (1, '%d ', movieStep); 
   
   % Increase MPM counter
   MPMCount = MPMCount + 1;
    
   selectedCells = zeros (size (MPM, 1), 2);
   selectedCells(:,:) = MPM (:, (2 * MPMCount - 1):(2 * MPMCount));
      
   % Remove the zero rows
   selectedCells = find (selectedCells(:,1) & selectedCells(:,2)); 
   
   % Go the image directory and fetch the current frame
   cd (imageDirectory);
   name = char (imageNameList (movieStep)); 
   nowImgH = imreadnd2 (name, 0, intensityMax);
	
   % If the user specified a size for the movie, that's the size we are going to use
   if ~isempty (figureSize)
      figure ('Position', figureSize), imshow (nowImgH, []);
   else
      figure, imshow (nowImgH, []);
   end

   % Hold on the figure for the dragtail drawing
   hold on;
   
   % One colour per time step (dragTailLength tells you how many
   % timesteps there are)
   colorMap = jet (dragTailLength + 1);

   % Initialize the dragtail counter
   colorCount = 0;
        
   % Loop through the previous frames to generate the dragtails
   for iCount = (2 * (MPMCount - dragTailLength)) : 2 : (2 * MPMCount)
       
      % Increase dragtail counter
      colorCount = colorCount + 1;
      
      % Get the needed frames from MPM
      selectedFrames = MPM (selectedCells, iCount-3 : iCount);
      
      % Find the ones that contain zeros and make the whole row 0
      [rows, cols] = find (selectedFrames == 0);
      rows = unique (rows);
      selectedFrames (rows,:) = 0;

      % Overlay the dragtails on the figure
      for hCount = 1 : size (selectedFrames,1)
         if selectedFrames (hCount,1) ~= 0
	        ph = [];
	        ph = plot (selectedFrames (hCount, 1:2:3), selectedFrames (hCount, 2:2:4));
	        set (ph, 'Color', colorMap (colorCount,:));
	        clear ph;
         end   % if selectedFrames
      end   % for hCount
   end   % for iCount
        
   % Plot the generated dragtail points on the figure
   plot (selectedFrames (:,3), selectedFrames (:,4), 'r.');
   hold off;

   % Go to the directory where the movies are stored
   cd (savePath);
  
   % Add the current figure to the movie
   makeQTMovie ('addaxes', gca);
   %F = getframe (gca);
   %mov = addframe (mov, F);
   
   % Close the figure
   close;
   
end   % for movieStep

% finalize the movie and write it to disk
makeQTMovie ('finish');
%mov = close(mov);

% Let the user know we have finished
fprintf (1, '\nptMovieMaker: finished generating movie.\n');
