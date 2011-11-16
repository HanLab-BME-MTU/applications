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

noiseProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(noiseProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',noiseProc.getName());
else
    wtBar=-1;
end
% Get channel paths and initialize process paths and output dirs
nChan = numel(movieData.channels_);

% Setup the  input directories
inFilePaths = cell(1,nChan);
for j = p.ChannelIndex;    
    %Create string for current directory
    inFilePaths{1,j} = movieData.getChannelPaths{j};
end
noiseProc.setInFilePaths(inFilePaths);

% Setup the output directories
outFilePaths = cell(2,nChan);
pfName = 'noise_for_channel_'; 
dName = 'cropped_images_for_channel_'; 
for  j = p.ChannelIndex;   
    outFilePaths{1,j} = [p.OutputDirectory filesep pfName num2str(j) '.mat'];
    outFilePaths{2,j} = [p.OutputDirectory filesep dName num2str(j)];
    %Check/create directory
    mkClrDir( outFilePaths{2,j})               
end
noiseProc.setOutFilePaths(outFilePaths);

%% --------------- Estimating noise ---------------%%% 

% Retrieve various useful constants
imageFileNames = movieData.getImageFileNames;
nFrames = p.lastImage - p.firstImage+1;

% Anonymous functions for reading input/output
outImage=@(chan,frame) [outFilePaths{2,chan} filesep imageFileNames{chan}{frame}];
logMsg = @(chan) ['Please wait, calculating noise for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];

disp('Starting estimating noise...')
nTot = sum(nFrames(p.ChannelIndex));
noiseLog=cell(numel(p.ChannelIndex),1);
for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Display log
    noiseLog{i} = sprintf('Channel %g: %s\n',iChan,inFilePaths{1,iChan});
    disp(logMsg(iChan))
    disp(inFilePaths{1,iChan});
    disp('Result will be saved as :')
    disp(outFilePaths{1,iChan});
    disp('Cropped images will be saved under :')
    disp(outFilePaths{2,iChan});
    
    % Read the first image, get the crop dimensions and initialize stack
    frameRange = p.firstImage(iChan):p.lastImage(iChan);
    currImage = movieData.channels_(iChan).loadImage(1);
    dummy=imcrop(currImage,p.cropROI);
    cropSize = size(dummy);
    stack = zeros([cropSize frameRange(end)-frameRange(1)],class(currImage));

    disp('Loading stack...');
    tic
    for j = 1:numel(frameRange)                
        %Load the current image, crop and save it
        iFrame = frameRange(j);   
        stack(:,:,j) = imcrop(movieData.channels_(iChan).loadImage(iFrame),p.cropROI);
        imwrite(stack(:,:,j),outImage(iChan,iFrame),'tif');
        
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = sum(nFrames(1:i-1))+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) ...
                timeMsg(tj*nTot/nj-tj)]));
        end  
    end

    % Calculate noise model
    disp('Calculating noise parameters...');
    [I0,sDN,GaussRatio] = calculateStackNoiseParam(stack,...
        movieData.camBitdepth_,p.filterSigma(iChan)); %#ok<ASGLU,NASGU>
    
    % Save results
    save(outFilePaths{1,iChan},'I0','sDN','GaussRatio');
    
    % Create channel log fot output
    noiseLog{i} = [noiseLog{i} ...
        sprintf(['Noise model parameters\n' ...
        'Average background intensity\t: %2.5f +/- %2.5f\n'...
        'Gauss ratio\t\t\t: %2.4f\n'],I0,sDN,GaussRatio)];
    
end
% Close waitbar
if ishandle(wtBar), close(wtBar); end

% Create process report
procLog=[sprintf('Noise model calibration detection summary\n\n') noiseLog{:}];
disp(procLog);
fid=fopen([p.OutputDirectory filesep 'NoiseModelCalibrationSummary.txt'],'w+');
fprintf(fid,procLog);
fclose(fid);

disp('Finished calculating noise!')