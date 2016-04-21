function importCellMaskOO(processOrMovieData,varargin)
%importCellMask imports a previously-defined cell mask to movieData
%
%Khuloud Jaqaman, March 2015

%% Input

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(processOrMovieData,varargin{:});
paramsIn = ip.Results.paramsIn;

% Use new API to allow the Process to be passed directly
[movieData, maskProc] = getOwnerAndProcess(processOrMovieData,'ImportCellMaskProcess',true);


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

