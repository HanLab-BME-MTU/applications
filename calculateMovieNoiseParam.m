function calculateMovieNoiseParam(movieData,varargin)
% calculateMovieNoiseParam calculates the noise model of a movie
%
% calculateMovieNoiseParam returns the mean backgroung intensity (I0), the
% mean standard deviation (sDN) and GaussRatio, which indicates how well
% the dark noise of the camera approximates a normal distribution. More 
% precisely, GaussRatio corresponds to the ratio:
%
%                                std(background image)
%                           -------------------------------,
%                            std(low-pass filtered image)
%
% where 'low-pass filtered image' is the result of the convolution of the 
% image with a  Gaussian kernel with standard deviation = sigma.
%
% SYNOPSIS [I0,sDN,GaussRatio]=calculateMovieNoiseParam(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT  
%
% REMARK   calculateMovieNoiseParam loads all the selected images into 
%          memory. For this reason, adapt image number and size to your 
%          machine's amount of memory.
%
% Sebastien Besson 5/2011
% Adapted from fsmCalcNoiseParam

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous noise estimation processes                                                                              
iProc = movieData.getProcessIndex('NoiseEstimationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(NoiseEstimationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%% --------------- Initialization ---------------%%
wtBar = waitbar(0,'Initializing...');

% Get channel paths and initialize process paths and output dirs
imDirs = movieData.getChannelPaths;
inFilePaths = cell(1,numel(movieData.channels_));
outFilePaths = cell(1,numel(movieData.channels_));
cropImagesDir = cell(1,numel(movieData.channels_));

% Define string extension
pfName = 'noise_for_channel_'; 
dName = 'cropped_images_for_channel_'; 

for j = p.ChannelIndex;    
    %Create string for current directory
    inFilePaths{j} = imDirs{j};
    cropImagesDir{j} = [p.OutputDirectory filesep dName num2str(j)];    
    outFilePaths{j} = [p.OutputDirectory filesep pfName num2str(j) '.mat'];
   
    %Check/create directory
    mkClrDir(cropImagesDir{j})               
end
movieData.processes_{iProc}.setOutFilePaths(outFilePaths)

% Load external file (to be checked!!)
if ~isempty(p.loadExternalFile)
    for i = 1:numel(p.ChannelIndex)
        iChan = p.ChannelIndex(i);
        waitbar(0,wtBar,['Please wait, loading noise data for channel ' ...
            num2str(iChan) ' ...']);
        disp('Loading external file:')
        load(p.loadExternalFile)
        save(outFilePaths{iChan},'I0','sDN','GaussRatio')
    end
end

%% --------------- Estimating noise ---------------%%% 

% Retrieve various useful constants
imageFileNames = movieData.getImageFileNames;
nFrames = p.lastImage - p.firstImage+1;

% Anonymous functions for reading input/output
inImage=@(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];
outImage=@(chan,frame) [cropImagesDir{chan} filesep imageFileNames{chan}{frame}];
logMsg = @(chan) ['Please wait, calculating noise for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];

disp('Starting estimating noise...')
tic;
nTot = sum(nFrames(p.ChannelIndex));
for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Display log
    disp(logMsg(iChan))
    disp(inFilePaths{iChan});
    disp('Cropped images will be saved as :')
    disp(cropImagesDir{iChan});
    
    % Read the first image, get the crop dimensions and initialize stack
    frameRange = p.firstImage(iChan):p.lastImage(iChan);
    currImage = imread(inImage(iChan,1));
    dummy=imcrop(currImage,p.cropROI(iChan,:));
    cropSize = size(dummy);
    stack = zeros([cropSize frameRange(end)-frameRange(1)],class(currImage));

    fprintf(1,'Loading stack...');
    for j = 1:numel(frameRange)
        if mod(j,5)==0
            tj=toc;
            nj = sum(nFrames(1:i-1))+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) ...
                timeMsg(tj*nTot/nj-tj)]));
        end
                
        %Load the current image, crop and save it
        iFrame = frameRange(j);   
        currImage = imread(inImage(iChan,iFrame));
        stack(:,:,j) = imcrop(currImage,p.cropROI(iChan,:));
        imwrite(stack(:,:,j),outImage(iChan,iFrame),'tif');
    end

    % Calculate noise model
    fprintf(1,'Calculating noise parameters...');
    [I0,sDN,GaussRatio] = calculateStackNoiseParam(stack,...
        movieData.camBitdepth_,p.filterSigma(iChan));
    
    % Save results
    save(outFilePaths{iChan},'I0','sDN','GaussRatio');
end
% Close waitbar
close(wtBar);

%% ------ Finish - Save parameters and movieData ----- %%

%Set process date and time and save the movieData
movieData.processes_{iProc}.setDateTime;
movieData.save; 

disp('Finished calculating noise!')
