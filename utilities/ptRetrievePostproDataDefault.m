function [handlesOut, resultOut] = ptRetrievePostproDataDefault (filePath, handlesIn)
% ptRetrievePostproDataDefault loads ptPostpro data from the files indicated by
% filePath and stores these in the handles struct
%
% SYNOPSIS       ptRetrievePostproDataDefault (filePath, handlesIn)
%
% INPUT          filePath : path to and name of the file to process
%                handlesIn : GUI handle struct containing the postpro info
%
% OUTPUT         handlesOut : GUI handle struct containing the updated postpro info
%                result : result of the operation (0 = good, 1 = error)
%
% DEPENDENCIES   ptRetrievePostproDataDefault uses {nothing}
%                                  
%                ptRetrievePostproDataDefault is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Aug 04          Initial Release

% Get the directory part of the string
[pathString, filename, ext, version] = fileparts (filePath);

% Initialize result
resultOut = 0;
handlesOut = [];

% Do nothing in case the user doesn't select anything
if filename == 0
    return
else
   
    % Determine image file path from jobvalues path
    cd (pathString); cd ('..');
    imageFilePath = pwd;

    % Load MPM matrix and assign to handles struct
    load (filePath);
    handlesIn.MPM = MPM;

    % Start with the default post processing structure
    handlesIn.postpro = handlesIn.defaultPostPro;
    
    % Load the jobvalues file
    if exist ([pathString filesep 'jobvalues.mat'], 'file')
       load ([pathString filesep 'jobvalues.mat']);
       handlesIn.jobvalues = jobvalues;
    else
       % It might be one directory back if it is a processed MPM file
       cd ..;
       if exist ([pathString filesep 'jobvalues.mat'], 'file')
          load ([pathString filesep 'jobvalues.mat']);
          handlesIn.jobvalues = jobvalues;
       else
          resultOut = 1;
          h = errordlg ('The file jobvalues.mat does not exist...');
          uiwait(h);          % Wait until the user presses the OK button
          return;
       end
    end

    % Load the cell properties file if it exists
    if exist ([pathString filesep 'cellProps.mat'],'file')
       load([pathString filesep 'cellProps.mat']);
       if exist('cellProps','var')
          handlesIn.postpro.cellProps = cellProps;
       end
    else
       resultOut = 1; 
       h = errordlg('The file cellProps.mat does not exist. Please make sure it is present...');
       uiwait(h);          % Wait until the user presses the OK button
       return;
    end

    % Load the cluster properties file if it exists
    if exist ([pathString filesep 'clusterProps.mat'],'file')
       load([pathString filesep 'clusterProps.mat']);
       if exist('clusterProps','var')
          handlesIn.postpro.clusterProps = clusterProps;
       end
    else
       resultOut = 1; 
       h = errordlg('The file clusterProps.mat does not exist. Please make sure it is present...');
       uiwait(h);          % Wait until the user presses the OK button
       return;
    end

    % Load the frame properties file if it exists
    if exist ([pathString filesep 'frameProps.mat'],'file')
       load([pathString filesep 'frameProps.mat']);
       if exist('frameProps','var')
          handlesIn.postpro.frameProps = frameProps;
       end
    else
       resultOut = 1; 
       h = errordlg('The file frameProps.mat does not exist. Please make sure it is present...');
       uiwait(h);          % Wait until the user presses the OK button
       return;
    end

    % Set the HOME env variable if not set already
    home = getenv('HOME');
    if isempty (home)
       if ispc
          home = 'H:';
       else
          home = '/tmp';
       end
    end
    cd (home);
    
    % Here is where a new data subdirectory has to be created. 
    dataDirectory = 'ptData';

    % If it doesn't exist yet, create it in the results directory
    if ~exist (dataDirectory, 'dir')
       mkdir (home, dataDirectory);
    end

    % Save it in the handlesIn struct and tell the loop we're done
    handlesIn.saveallpath = [home filesep dataDirectory];

    % Store the size of the image
    cd (imageFilePath);
    info = imfinfo (handlesIn.jobvalues.imagename);
    handlesIn.postpro.rowsize = info.Height;
    handlesIn.postpro.colsize = info.Width;

    % Get the imagename without .tif
    imageNameNoTiff = regexprep(handlesIn.jobvalues.imagename, '.tif', '', 'ignorecase');

    % Now we have to fill up the rest of the postpro structure with
    % our previously found data and parameters
    handlesIn.selectedcells = [];
    handlesIn.postpro.imagepath = imageFilePath;
    handlesIn.postpro.increment = handlesIn.jobvalues.increment;
    handlesIn.postpro.firstimg = handlesIn.jobvalues.firstimage;
    handlesIn.postpro.lastimg = handlesIn.jobvalues.lastimage;
    handlesIn.postpro.maxdistpostpro = handlesIn.jobvalues.maxsearch;
    handlesIn.postpro.plotfirstimg = handlesIn.jobvalues.firstimage;
    handlesIn.postpro.plotlastimg = handlesIn.jobvalues.lastimage;
    handlesIn.postpro.selectedcells = [];
    handlesIn.postpro.moviefirstimg = handlesIn.jobvalues.firstimage;
    handlesIn.postpro.movielastimg = handlesIn.jobvalues.lastimage;
    handlesIn.postpro.jobpath = pathString;
    handlesIn.postpro.imagename = handlesIn.jobvalues.imagename;
    handlesIn.postpro.imagenamenotiff = imageNameNoTiff;
    handlesIn.postpro.imagenameslist = handlesIn.jobvalues.imagenameslist;
    handlesIn.postpro.intensitymax = handlesIn.jobvalues.intensityMax;    

    % These are new additions which won't be in the older jobs
    if ~isempty (handlesIn.jobvalues.timeperframe)
       handlesIn.postpro.timeperframe = handlesIn.jobvalues.timeperframe;
    end

    if ~isempty (handlesIn.jobvalues.mmpixel)
       handlesIn.postpro.mmpixel = handlesIn.jobvalues.mmpixel;
    end
end

% Return the modified handlesIn struct
handlesOut = handlesIn;