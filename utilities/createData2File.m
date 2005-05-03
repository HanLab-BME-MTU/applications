function [data2File, data2FileName] = createData2File(movieName,directory)
%CREATEDATA2FILE sets up the necessary fields for a data2-file
%
% SYNOPSIS [data2File,data2FileName] = createData2File(movieName,directory)
%
% INPUT movieName   movie for which the data2 file should be created
%       directory  (opt) directory name. Default: pwd
%
% OUTPUT data2File: structure with fields into which to fill the location
%                   of all the files of the directory
%        data2FileName : name of the data2 file for saving. data2FileName
%                        will be [movieName, '-data2-', nowString, '.mat']
%                        at the same time, also a logfile will be created,
%                        named [data2FileName, '.log']
%
% c: jonas 5/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================
% TEST INPUT
%======================

if nargin < 1 || isempty(movieName)
    error('you have to specify a movie name!')
end

if nargin < 2 || isempty(directory)
    directory = pwd;
end

%=======================

%=============================
% CREATE STRUCTURE AND NAMES
%=============================

% check movieName for fileExtension
extension = findstr('movieName','.');
if ~isempty(extension)
    % check if it's really an extension
    if extension == length(movieName)-3
        movieName = movieName(1:end-4);
    else
        error('Potentially improper movieName!')
    end
end

data2FileName = [movieName,'-data2-',nowString,'.mat'];
data2LogFileName = [movieName,'-data2-',nowString,'.mat.log'];

% these are all the fields I can think of right now. 
data2File = struct('dataProperties','',...
    'r3dMovieHeader','','correctionData','',...
    'movieName','','filteredMovieName','',...
    'slist','','idlist','','idlist_L','',...
    'idlisttrack','','idlisttrack_L','',...
    'lastResult','','synthIdlist','','synthSlist','');

% make files
oldDir = cd(directory);
save(data2FileName,'data2File');
logID = fopen([data2LogFileName],'w');
fprintf(logID, '%s: logFile created/n/n',nowString);
fclose(logID);

% go back to old Dir
cd(oldDir);