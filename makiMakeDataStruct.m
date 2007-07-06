function dataStruct = makiMakeDataStruct(rawMovieName, rawMoviePath, projectName, dataFilePath)
%MAKIMAKEDATASTRUCT creates the overall data structure for mammalian kinetochore analysis
%
% SYNOPSIS: dataStruct = makiMakeDataStruct
%
% INPUT   : rawMovieName : filename of raw movie.
%           rawMoviePath : pathname of raw movie
%           projectName  (opt): identifier of the project. If empty,
%                        projectName will be rawMovieName minus the file
%                        extension
%           dataFilePath (opt): pathname for dataFile. Specify if the
%                        dataFile stored in a different path from the
%                        rawMovie 
%
% OUTPUT dataStruct: Data structure with fields 
%           projectName : identifier for the project. DataFileName will be
%                       projectName-data-currentDate.mat
%           rawMovieName/path : file and path name for raw movie (can be
%                       different from dataFilePath
%			dataFilePath: path where the data file is stored. This will
%                       always be a full path. However, makiLoadDataFile
%                       will keep overwriting this field on each load
%           dataFileName: filename of the data file
%           xxx/xxxName; xxx=e.g. dataProperties: data/fileName pairs. When
%                       makiSaveDataFile.m will save xxx into
%                       xxxName, makiLoadDataFile will write the content of
%                       the file xxxName into xxx. The path is assumed to
%                       be the same as the dataFilePath (can be changed in
%                       the future)
%			status: status of the image analysis. 1: done, 0: not done, -1:
%                       to be done
%           statusHelp: cell array with row for every entry in status
%                       {name of analysis, date done, comments}
%     
%
% REMARKS 
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 27-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================
%% CHECK INPUT
%===================

% check input arguments
if nargin < 2 || isempty(rawMovieName) || isempty(rawMoviePath)
    error('not enough input arguments - need nonempty moviename/path')
end
if nargin < 3 || isempty(projectName)
    [dummy,fileBody] = fileparts(rawMovieName);
    % insert switches here to make names more complicated 
    % Also, talk about naming conventions 
    projectName = fileBody; 
end
    
if nargin < 4 || isempty(dataFilePath)
    dataFilePath = rawMoviePath;
end


%====================
%% MAKE DATASTRUCT
%====================

% comments
% -makiData- will distinguish this data from the yeast data files
% logfile-name is [projectName,'_analysis.log'] to avoid confusion
% NEED: name/structure of file that will store the tracks. idlist would be
%       nice for visualization with labelgui, in case that is useful
% expand status along with expansion of dataFile fields. Use date instead
%       of nowString to allow checking for earlier/later at some point

dataStruct = struct('projectName',projectName,...
    'rawMovieName',rawMovieName,...
    'rawMoviePath',rawMoviePath,...
    'dataFileName',[projectName,'-makiData-',nowString,'.mat'],...
    'dataFilePath',dataFilePath,...
    'dataPropertiesName',['dataProperties_',projectName,'.mat'],...
    'dataProperties',[],...
    'movieHeaderName',['movieHeader_',projectName,'.mat'],...
    'movieHeader',[],...
    'initCoordName',['initCoord_',projectName,'.mat'],...
    'initCoord',[],...
    'tracksName',['tracks_',projectName,'.mat'],...
    'tracks',[],...
    'planeFitName',['planeFit_',projectName,'.mat'],...
    'planeFit',[],...
    'sisterListName',['sisterList_',projectName,'.mat'],...
    'sisterList',[],...
    'slistName',['slist_',projectName,'.mat'],...
    'slist',[],...
    'status',[0,0,0,0,0,0,0]',...
    'statusHelp',...
    {{'cropMovie','','';...
    'dataFileSaved','','';...
    'initCoord','','';...
    'tracks','','';...
    'planeFit','','';...
    'sisterList','','';...
    'slist','',''}});

% discontinued
%  'filteredMovieName',['filtered_',projectName,'.mat'],...