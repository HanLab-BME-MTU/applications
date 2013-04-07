function calc3DMoviePhotobleaching(movieData,paramsIn)

%Another quick and sloppy function to get this shit finished!!!!

%% ------------------- Parameters ---------------- %%

iProcChan = 1;%Hard-coded channels for non-channel specific mask and post-processing association. Yeah I know this is a stupid way to do it, but who fucking cares.
%iPostProcChan = 2;
fName = 'bleach_calc';%File name for output

%Figure saving parameters
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTif = {'-r100',...% dpi = 100
        '-dtiff'};% use tiff format
    
%Field names and function handles for image int stats
fieldNames = {'meanIntensity','medianIntensity','stdIntensity','nVoxels'};
fHans = {@mean,@median,@std,@numel};
nFields = numel(fieldNames);    
iFit = 1;%Index of the parameter to fit to.

fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));     %Double-exponential function for fitting
bInit = [1 0 1 0]; %Initial guess for fit parameters.

figName = 'photobleach correction fit'; %Name for saving figure to files

%% -------------------- Input -------------- %%

if nargin < 1 || ~isa(movieData,'MovieData3D')
    error('You must input a valid MovieData3D object as the first input!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous processes from this function
iProc = movieData.getProcessIndex('PhotoBleachCalcProcess3D',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(PhotoBleachCalcProcess3D(movieData));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);


%% ------------------- Init --------------- %%


%Make sure the movie has been segmented and has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(iProcChan)%TEMP? how do deal with three-input channels for mask and one output?? ALways associate with channel 1?
    error('The input MovieData does not have valid masks for the selected channel!');
end

%Get mask file names and directory
maskDir = movieData.processes_{iSegProc}.outFilePaths_{1,iProcChan};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(iProcChan);
% %And the post-process added pixel directory and file-names
% ppDir = movieData.processes_{iSegProc}.outFilePaths_{1,iPostProcChan};
% ppNames = movieData.processes_{iSegProc}.getOutMaskFileNames(iPostProcChan);

%Get image locations and info
imDir = movieData.getChannelPaths(p.ChannelIndex);
imNames = movieData.getImageFileNames(p.ChannelIndex);
nChan = numel(p.ChannelIndex);
if ~isempty(movieData.eventTimes_)
    nFrames = eventTimes_(2);
else
    nFrames = movieData.nFrames_;
end

if ~isempty(movieData.timeInterval_)
    timePoints = (0:1:nFrames-1) * movieData.timeInterval_;     %time data
else
    timePoints = (0:1:nFrames-1);
end


%Get the xy and z pixel sizes for scaling the masks
if ~isempty(movieData.pixelSize_) && ~isempty(movieData.zSpacing_)
    pixXY = movieData.pixelSize_;
    pixZ = movieData.zSpacing_;
    %hasSizes = true;    
else
    %warn the user, and assume unit pixel sizes.
    warning('Migration3D:MissingVoxelDimensions',...
        'Pixel XY size and Z spacing not specified in input movieData! Intensity analysis will assume symmetric voxels of size 1nm!');
    pixXY = 1;
    pixZ = 1;
    %hasSizes = false;
end


outDir = p.OutputDirectory;
mkClrDir(outDir);

disp('Starting photobleaching rate calculation...')

for j = 1:nFields
    bleachStats.(fieldNames{j}) = nan(nFrames,nChan);
end


if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, calculating bleach rate....');
end


%% --------------- Intensity Stats ----------------- %%
%Go through each frame and get the masked intensity statistics


for iFrame = 1:nFrames
    
    %Load the mask and post processed pixels and make the voxels symmetric to get the stretched size
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);            
    
    for iChan = 1:nChan            
        currIm = stackRead([imDir{iChan} filesep imNames{iChan}{iFrame}]);                                    
        %Get the pixels so we only perform the indexing operation once,
        %cast to double to avoid rounding error
        currPix = double(currIm(currMask(:)));
        
        for iField = 1:nFields
            bleachStats.(fieldNames{iField})(iFrame,iChan) = fHans{iField}(currPix);            
        end
        
    end  
    
    if ~p.BatchMode
        waitbar(iFrame/nFrames,wtBar)
    end    
    
end

%% ---------- Bleach Rate Calculation ------------ %%

%Fit function to ratio timeseries
fitOptions = statset('Robust','off','MaxIter',500,'Display','off');

for iChan = 1:nChan
    
    %Setup the parameter initial guesses based on the first and last
    %intensities, assuming equal contributioin from both exponentials
    bInitCurr([1 3]) = bleachStats.(fieldNames{iFit})(1,iChan)/2;%Aplitudes based on first frame
    bInitCurr([2 4]) = log(bleachStats.(fieldNames{iFit})(end,iChan) / (2*bInitCurr(1)))/timePoints(end);
    
    [bleachStats.bFit(:,iChan),bleachStats.resFit(:,iChan),bleachStats.jacFit(:,:,iChan),bleachStats.covFit(:,:,iChan),bleachStats.mseFit(iChan)] ...
        = nlinfit(timePoints(:),bleachStats.(fieldNames{iFit})(:,iChan),fitFun,bInitCurr,fitOptions);
    %Get confidence intervals of fit and fit values
    [bleachStats.fitValues(:,iChan),bleachStats.deltaFit(:,iChan)] = nlpredci(...
        fitFun,timePoints(:),bleachStats.bFit(:,iChan),bleachStats.resFit(:,iChan),'covar',bleachStats.covFit(iChan),'mse',bleachStats.mseFit(iChan));

    %Check the fit jacobian
    [dummy,bleachStats.R(iChan)] = qr(bleachStats.jacFit(iChan),0); %#ok<ASGLU>
    if ~p.BatchMode && condest(bleachStats.R(iChan)) > 1/(eps(class(bleachStats.bFit(:,iChan))))^(1/2)        
        warndlg('WARNING: The photobleach correction fit is not very good!!!!! SHIT!!!!!')
    end
end


%% ---------- Figure Creation ------------- %%

for iChan = 1:nChan

    if p.BatchMode
        fitFig = figure('Visible','off');
    else
        fitFig = figure;
    end

    
    subplot(2,1,1)
    hold on
    title('Photobleach Correction Fit')
    if ~isempty(movieData.timeInterval_)
        xlabel('Time, seconds')
    else
        xlabel('Frame Number')
    end
    ylabel(fieldNames{iFit})
    plot(timePoints,bleachStats.(fieldNames{iFit})(:,iChan),'.-')
    plot(timePoints,bleachStats.fitValues(:,iChan),'r')
    yl = ylim;
    plot(timePoints,bleachStats.fitValues(:,iChan)+bleachStats.deltaFit(iChan),'--r')
    legend(fieldNames{iFit},'Fit','Fit 95% C.I.')
    plot(timePoints,bleachStats.fitValues(:,iChan)-bleachStats.deltaFit(iChan),'--r')
    ylim(yl);
    xlim(timePoints([1 end]))
    
    subplot(2,1,2)
    hold on
    plot(timePoints,bleachStats.resFit,'.-')
    if ~isempty(movieData.timeInterval_)
        xlabel('Time, seconds')
    else
        xlabel('Frame Number')
    end
    ylabel('Fit Residuals')
    xlim(timePoints([1 end]))
    plot(xlim,[0 0],'--k')
    
    
    
    currFigFile = [p.OutputDirectory filesep figName ' channel ' num2str(iChan)];                   
    hgsave(fitFig,currFigFile)
    print(fitFig,currFigFile,pOptTif{:})
    print(fitFig,currFigFile,pOpt{:})
    
    %Log this file name in the parameter structure
    p.figName = figName;

    if p.BatchMode && ishandle(fitFig) %make sure user hasn't closed it.
        close(fitFig)
    end
end


%% ------Output/Finalization----- %%

outFile = [outDir filesep fName '.mat'];
save(outFile,'bleachStats');

%Set the out file paths, times etc
movieData.processes_{iProc}.setOutFilePath(iProcChan,outFile);
movieData.processes_{iProc}.setDateTime;
movieData.processes_{iProc}.setPara(p);
movieData.save; %Save the new movieData to disk

disp('Finished photobleach rate calc!')

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end
