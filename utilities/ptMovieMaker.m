function result = ptMovieMaker (radioButtons, handles) 
% ptMovieMaker makes movies from the information gathered in the analysis
%
% SYNOPSIS       ptMovieMaker (ptPostpro, MPM)
%
% INPUT          radioButtons : a structure which contains the status of the radio
%                               buttons on the GUI
%                handles      : struct containing the GUI data
%
% OUTPUT         result : result of the operation (0 = ok, 1 = error)         
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
% Andre Kerstens        Aug 04          Bug fix: forgot to close fig after addframe QT movie
% Andre Kerstens        Sep 04          Rewrite for the new data structures

% Initialize result
result = 0;

% Make sure only 1 MPM is present
if length(handles.allMPM) > 1
    result = 1;
    return;
end

% Get the latest data from the handles
MPM = handles.allMPM{1};    % We should have received only one MPM
cellProps = handles.allCellProps;
clusterProps = handles.allClusterProps;
frameProps = handles.allFrameProps;
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
movieStartFrame = guiData.moviefirstimg;
movieEndFrame = guiData.movielastimg;
if movieEndFrame > jobData(1).lastimg
    movieEndFrame = jobData(1).lastimg;
end

% Get start and end frames and increment value
startFrame = jobData(1).firstimg;
endFrame = jobData(1).lastimg;
increment = jobData(1).increment;
numberOfFrames = ceil((movieEndFrame - movieStartFrame) / increment) + 1;

% First assign all the jobdata fields to a meaningful variable
imageName = jobData(1).imagename;
imageDirectory = jobData(1).imagefilepath;
imageNameList = jobData(1).imagenameslist;
intensityMax = jobData(1).intensitymax;
figureSize = jobData(1).figuresize;

% Determine the movie type
if radioButtons.movietypeavi
   movieType = 1; % AVI
elseif radioButtons.movietypeqt
   movieType = 2; % QT
else
   movieType = 1; % AVI
end

% How many frames does the user wants to see the tracks of?
dragTailLength = guiData.dragtail;

% We need enough frames to show the length of the dragtail
if (movieEndFrame - startFrame + 1) < dragTailLength
   h = errordlg (['Not enough frames for dragtail length of ' num2str(dragTailLength) '. Exiting...']);
   uiwait (h);
   result = 1;
   return;
end

% What is the name of the movie
dragTailPath = guiData.dragtailfile;
[dir,dragTailFile,dummy,ext] = fileparts(dragTailPath);

% We need prior images for the dragTail movie
if ceil ((movieStartFrame - startFrame + 1) / increment) < dragTailLength + 1
    movieStartFrame = ((dragTailLength + 1) * increment) + movieStartFrame;
end

% Get save path
SaveDir = guiData.savedatapath;

% If it doesn't exist yet, create it in the results directory
if ~exist (SaveDir, 'dir')
   mkdir (SaveDir);
end

% Go to the directory where the movie will be saved
cd (SaveDir);

% Initialize the movie
if movieType == 1   % AVI
   mov = avifile ([dragTailPath '.avi']);
elseif movieType == 2   % QT
   MakeQTMovie ('start', [dragTailPath '.mov']);
else
   h = errordlg (['Unknown movie type. Please choose QT or AVI. Exiting...']);
   uiwait (h);
   result = 1;
   return;
end

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
    
   selectedCells = zeros(size (MPM, 1), 2);
   selectedCells(:,:) = MPM(:, (2 * MPMCount - 1):(2 * MPMCount));
   
   % Store the coordinates for later
   selectedDots = selectedCells; 
   selectedDots(find(selectedDots(:,1) == 0 & selectedDots(:,2) == 0),:) = [];
   
   % Remove the zero rows
   selectedCells = find (selectedCells(:,1) & selectedCells(:,2)); 
   
   % Go the image directory and fetch the current frame
   cd (imageDirectory);
   name = char (imageNameList (movieStep)); 
   nowImgH = imreadnd2 (name, 0, intensityMax);
	
   if movieStep == movieStartFrame
       if ~isempty (figureSize)
          h_fig = figure ('Position', figureSize);
       else
          h_fig = figure;
       end
   end
   
   % If the user specified a size for the movie, that's the size we are going to use
   if ~isempty (figureSize)
      imshow (nowImgH, []);
   else
      if movieStep == movieStartFrame
          imshow (nowImgH, []);
      else
          imgPtr=findall(gca,'Type','Image');
          set(imgPtr,'CData',nowImgH);
          set(imgPtr,'CDataMapping','scaled')
          set(gca,'CLimMode','auto');
          refresh;
      end
   end  % ~isempty (figureSize)

   % Check whether there are already dots plotted and if yes delete them
   currentH=findall(gca,'Tag','dots');
   if ~isempty(currentH)
      delete(currentH);
   end
   
   % Check whether there are already tails plotted and if yes delete them
   currentH=findall(gca,'Tag','tails');
   if ~isempty(currentH)
      delete(currentH);
   end
   
   % Hold on the figure for the dragtail drawing
   hold on;
   
   if get(handles.GUI_fm_incltracks_rb,'Value') == 1
   
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

          % Overlay the dragtails on the figure; each with its own color
          for hCount = 1 : size (selectedFrames,1)
             if selectedFrames (hCount,1) ~= 0
                %ph = [];
                ph = plot (selectedFrames (hCount, 1:2:3), selectedFrames (hCount, 2:2:4));

                % Mark the tails as tails
                set (ph, 'Color', colorMap(colorCount,:),'Tag','tails','LineWidth',2);
                %clear ph;
             end   % if selectedFrames
          end   % for hCount
       end   % for iCount
   end  % get(handles.GUI_fm_incltracks_rb,'Value') == 1
   
   if get(handles.GUI_fm_inclcentromers_rb,'Value') == 1       
       % Plot the nuclei coordinates on the figure (red dots) and mark them
       dots = plot (selectedDots (:,1), selectedDots (:,2), 'r.');
       set (dots, 'Tag', 'dots');
   end  % get(handles.GUI_fm_inclcentromers_rb,'Value') == 1
   
   hold off;

   % Go to the directory where the movies are stored
   cd (SaveDir);
  
   % Save movie as tiff if needed
   if get(handles.GUI_fm_saveastiff_rb,'Value') == 1 
       print (h_fig, [SaveDir filesep [dragTailFile '_' num2str(movieStep) '.tif']],'-dtiff');
   end
   
   % Add the current figure to the movie
   if movieType == 1   % QT
      F = getframe (gcf);
      mov = addframe (mov, F);
   elseif movieType == 2   % AVI
      MakeQTMovie ('addaxes', gcf);
   else
      h = errordlg (['Unknown movie type. Please choose QT or AVI. Exiting...']);
      uiwait (h);
      result = 1;
      return;
   end    
end   % for movieStep

% finalize the movie and write it to disk
if movieType == 1   % QT
   mov = close(mov);
elseif movieType == 2   % AVI
   MakeQTMovie ('finish');
else
   h = errordlg (['Unknown movie type. Please choose QT or AVI. Exiting...']);
   uiwait (h);
   result = 1;
   return;
end

% Close figure
close;

% Let the user know we have finished
fprintf (1, '\nptMovieMaker: finished generating movie.\n');
