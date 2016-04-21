function importCellMaskOO(movieData,varargin)
%importCellMask imports a previously-defined cell mask to movieData
%
%Khuloud Jaqaman, March 2015

%% Input

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn = ip.Results.paramsIn;

%Get the indices of any previous cell mask import process                                                                        
%If the process doesn't exist, create it
iProc = movieData.getProcessIndex('ImportCellMaskProcess',1,0);
if isempty(iProc)
    iProc=numel(movieData.processes_)+1;
    movieData.addProcess(ImportCellMaskProcess(movieData));
end
maskProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(maskProc,paramsIn);

%% --------------- Initialization ---------------%%

%Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
mkClrDir(p.OutputDirectory)
for i = p.ChannelIndex;
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
maskProc.setOutFilePaths(outFilePaths);

%% --------------- cell mask import ---------------%%%

for i = p.ChannelIndex
    
    %Ask user where cell mask is currently stored
    [fileName,filePath,filterIndx] = uigetfile('*.tif',['Cell mask file for Channel ' num2str(i) ' of Movie ' movieData.movieDataFileName_(1:end-4)],movieData.roiMaskPath_);
    
    %save cell mask file and its original location
    if filterIndx == 0
        error('No file selected');
    else
        copyfile(fullfile(filePath,fileName),fullfile(p.OutputDirectory,['cellMask_channel_' num2str(i) '.tif']));
        save(outFilePaths{1,i} ,'fileName','filePath');
    end
    
end

