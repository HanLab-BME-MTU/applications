function cdCopyMovies
%CDCOPYMOVIES copies deltavision movies (and logfiles) to a new directory
%
% SYNOPSIS: cdCopyMovies
%
% INPUT 
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.6.0.324 (R2008a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 03-Apr-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ask user for old movie dir
[currentDir, oldDir] = cdBiodata(1);
oldMovieDir = uigetdir(pwd,'Please select a movie top-level directory');
if ~ischar(oldMovieDir)
    return
end
newMovieDir = uigetdir(oldMovieDir,'Please select a new directory');
if ~ischar(newMovieDir)
   return
end

% lookfor .r3d or .dv files
fileList = searchFiles('r3d$|dv$','^DIC',oldMovieDir);

selectIdx = listSelectGUI(fileList(:,1));
% perform selection
fileList = fileList(selectIdx,:);

% copyfile is a function that is close to the OS, thus it doesn't accept
% arrays as input. Loop.

nFiles = size(fileList,1);
progressText(0,'Copying movies') % Create text
for iFile = 1:nFiles
    % update progress text
    progressText((iFile-1)/nFiles+eps,sprintf('Copying (current movie: %s)',fileList{iFile,1}))
    % copy movie
    copyfile(fullfile(fileList{iFile,2},fileList{iFile,1}),newMovieDir);
    % copy associated files (such as log files)
    logFile = searchFiles([fileList{iFile,1}(1:end-4),'.+log'],'-data',fileList{iFile,2},0);
    if ~isempty(logFile)
        copyfile(fullfile(logFile{1,2},logFile{1,1}),newMovieDir);
    end
    % copy DIC movie
    dicFile = searchFiles('^DIC','',fileList{iFile,2},0);
    if ~isempty(dicFile)
        copyfile(fullfile(dicFile{1,2},dicFile{1,1}),newMovieDir);
    end
end
% close progressText
progressText(1,sprintf('Copying %i movies',nFiles));