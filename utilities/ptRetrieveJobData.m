function [allMPM, allCellProps, allClusterProps, allFrameProps, jobData, result] = ptRetrieveJobData (fileList, filesSelected)
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

% Initialize result of the whole operation (0 = okay)
result = 0;
filePath = cell(1);

% Get the data for all jobs
for iCount = 1 : length (filesSelected)
    
    % Get the filepath
    filePath(iCount) = cellstr(fileList{filesSelected(iCount)});

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
    
    % Set imagefile path and name
    jobData(iCount).imagefilepath = allJobvalues{iCount}.imagedirectory;
    jobData(iCount).imagename = allJobvalues{iCount}.imagename;
    
    % Store the size of the image
    info = imfinfo ([jobData(iCount).imagefilepath filesep jobData(iCount).imagename]);
    jobData(iCount).rowsize = info.Height;
    jobData(iCount).colsize = info.Width;
    
    % Figure size for movie generation
    jobData(iCount).figuresize = [];

    % Get the imagename without .tif
    jobData(iCount).imagenamenotiff = regexprep(jobData(iCount).imagename, '.tif', '', 'ignorecase');

    % Now we have to fill up the rest of the jobData structure with
    % our previously found data and parameters
    jobData(iCount).selectedcells = [];
    jobData(iCount).increment = allJobvalues{iCount}.increment;
    jobData(iCount).firstimg = allJobvalues{iCount}.firstimage;
    jobData(iCount).lastimg = allJobvalues{iCount}.lastimage;
    jobData(iCount).imagenameslist = allJobvalues{iCount}.imagenameslist;
    jobData(iCount).intensitymax = allJobvalues{iCount}.intensityMax;
    jobData(iCount).jobpath = pathString;
    jobData(iCount).timeperframe = allJobvalues{iCount}.timeperframe;
    jobData(iCount).mmpixel = allJobvalues{iCount}.mmpixel;
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
