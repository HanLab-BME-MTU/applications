function [allMPM, allCellProps, allClusterProps, allFrameProps, allValidFrames, jobData, result] = ptRetrieveJobData (fileList, filesSelected)
% ptRetrieveJobData loads data from the files indicated by
% filePath and stores these in several structs and cells
%
% SYNOPSIS       [allMPM, allCellProps, allClusterProps, allFrameProps, jobData, result] = 
%                                ptRetrieveJobData (fileList, filesSelected)
%
% INPUT          fileList : list of MPM files 
%                filesSelected : which files in the list should be processed
%
% OUTPUT         allMPM : cell containing all MPM read
%                allCellProps : cell containing all cellProps matrices
%                allClusterProps : cell containing all clusterProps matrices
%                allFrameProps : cell containing all frameProps matrices
%                allValidFrames : cell containing all validFrames arrays
%                jobData : struct containing all other data related to jobs
%                result : result of the operation (0 = good, 1 = error)
%
% DEPENDENCIES   ptRetrieveJobData uses {nothing}
%                                  
%                ptRetrieveJobData is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Aug 04          Initial Release
% Andre Kerstens        Nov 04          build in validFrames support

% Initialize result of the whole operation (0 = okay)
result = 0;
filePath = cell(1);

% filesSelected could be 'all' or a number. Find out first
if ischar(filesSelected)
   if num2str(filesSelected) == 'all'
      fileNumbers = length(fileList);
      method = 1;
   else
      result = 1;  % Error
      allMPM = [];
      allCellProps = [];
      allClusterProps = [];
      allFrameProps = [];
      allValidFrames = [];
      jobData = [];
      return;
   end
elseif isnumeric(filesSelected)
   fileNumbers = length(filesSelected);
   method = 2;
end

% Get the data for all jobs
for iCount = 1 : fileNumbers
    
    % Get the filepath
    if method == 1
       filePath(iCount) = cellstr(fileList{iCount});
    elseif method == 2
       filePath(iCount) = cellstr(fileList{filesSelected(iCount)});
    end

    % Load MPM matrix and assign to handles struct
    load (char(filePath(iCount)));
    allMPM{iCount} = MPM;

    % Start with the default post processing structure
    %handlesIn.postpro = handlesIn.defaultPostPro;

    % Get the directory part of the string
    [pathString, filename, ext, version] = fileparts (char(filePath(iCount)));
    
    % Load the jobvalues file
    [allJobvalues{iCount}, result] = ptReadValues (pathString, 'jobvalues.mat', 'jobvalues');
    if result == 1  % error
       return
    end

    % Load the cell properties file
    [allCellProps{iCount}, result] = ptReadValues (pathString, 'cellProps.mat', 'cellProps');
    if result == 1  % error
       return
    end
    
    % Load the cluster properties file
    [allClusterProps{iCount}, result] = ptReadValues (pathString, 'clusterProps.mat', 'clusterProps');
    if result == 1  % error
       return
    end
    
    % Load the frame properties file
    [allFrameProps{iCount}, result] = ptReadValues (pathString, 'frameProps.mat', 'frameProps');
    if result == 1  % error
       return
    end

    % Load the validFrames file
    [allValidFrames{iCount}, result] = ptReadValues (pathString, 'validFrames.mat', 'validFrames');
    if result == 1  % error
       return
    end
    
    % Set image filename
    jobData(iCount).imagename = allJobvalues{iCount}.imagename;    
    
    % Set image sizes
    if isfield(allJobvalues{iCount},'rowsize')
        jobData(iCount).rowsize = allJobvalues{iCount}.rowsize;
    end
    if isfield(allJobvalues{iCount},'colsize')
        jobData(iCount).colsize = allJobvalues{iCount}.colsize;
    end
    
    % Determine image file path from MPM path
    currentPath = pwd;
    cd (pathString); cd ('..');
    
    if isempty(dir(jobData(iCount).imagename))
       % Image not found: let's try one more directory down
       cd ('..');
    end
    
    if isempty(dir(jobData(iCount).imagename))
        % Still no images found
        % We continue loading, but certain functions like manual postpro
        % will not be possible
        jobData(iCount).imagefilepath = '';
        jobData(iCount).imagesavailable = 0;
    else   
        jobData(iCount).imagefilepath = pwd;
        jobData(iCount).imagesavailable = 1;

        % Store the size of the image
        %if ~isfield(jobData(iCount),'rowsize') | ~isfield(jobData(iCount),'colsize')
            info = imfinfo ([jobData(iCount).imagefilepath filesep jobData(iCount).imagename]);
            jobData(iCount).rowsize = info.Height;
            jobData(iCount).colsize = info.Width;
        %end

        % Figure size for movie generation
        jobData(iCount).figuresize = [];

        % Get the imagename without .tif
        jobData(iCount).imagenamenotiff = regexprep(jobData(iCount).imagename, '.tif', '', 'ignorecase');
    end
    
    % Reset the original path
    cd (currentPath);
    
    % Now we have to fill up the rest of the jobData structure with
    % our previously found data and parameters
    jobData(iCount).selectedcells = [];
    jobData(iCount).increment = allJobvalues{iCount}.increment;
    jobData(iCount).firstimg = allValidFrames{iCount}(1,1);
    jobData(iCount).lastimg = allValidFrames{iCount}(1,end);
    jobData(iCount).imagenameslist = allJobvalues{iCount}.imagenameslist;
    jobData(iCount).intensitymax = allJobvalues{iCount}.intensityMax;
    jobData(iCount).jobpath = pathString;
    jobData(iCount).timeperframe = allJobvalues{iCount}.timeperframe;
    jobData(iCount).mmpixel = allJobvalues{iCount}.mmpixel;
    jobData(iCount).nrbadframes = allValidFrames{iCount}(1,end)-length(allValidFrames{iCount}(1,:));
end

%-----------------------------------------------------------------------

function [outputValues, result] = ptReadValues (pathString, fileName, varName)
% ptReadValues reads a mat file and returns the values in outputValues
% result shows whether the operation was successful

% Initialize result
result = 0;

% Load the mat file
if exist ([pathString filesep fileName], 'file')
   load ([pathString filesep fileName]);
   if exist(varName,'var')
      outputValues = eval(varName);
   else
      result = 1;
      outputValues = [];
      h = errordlg ([varName 'could not be read...']);
      uiwait(h);          % Wait until the user presses the OK button
      return;
   end
else
   % It might be one directory back if it is a processed MPM file
   if exist ([pathString filesep '..' filesep fileName], 'file')
      load ([pathString filesep '..' filesep fileName]);
      if exist(varName,'var')
         outputValues = eval(varName);
      else
         result = 1;
         outputValues = [];
         h = errordlg ([varName 'could not be read...']);
         uiwait(h);          % Wait until the user presses the OK button
         return;
         end
   else
      result = 1;
      outputValues = [];
      h = errordlg (['The file ' fileName ' does not exist...']);
      uiwait(h);          % Wait until the user presses the OK button
      return;
   end
end
