function analyze3DMovieMaskedIntensities(movieData,paramsIn)
%ANALYZE3DMOVIEMASKEDINTENSITIES calls analyze3DImageMaskedIntensities on each frame of the input movie 
% 
% analyze3DMovieMaskedIntensities(movieData)
% analyze3DMovieMaskedIntensities(movieData,paramsIn)
% 
% This function is just a wrapper that checks and loads all input required
% and then calls analyze3DMovieMaskedIntensities on each image/frame of the
% input movie.
% 
% Input:
% 
%   movieData - A valid MovieData3D object describing the movie to analyze.
% 
%   paramsIn - A structure containing optional parameters, or parameter
%   name/value pairs. The possible paramter field names and values are:
% 
%       (FieldName->fieldValue)
%
%       ('CurvSampRad' -> Positive scalar) Specifies the local curvature
%       and intensity sampling radius in nanometers.
%       Optional. Default is 500nm.
%
%       ('PhotoBleachMeth'-> 'None' MovieData or 'Self')
%       String or MovieData object specifying the method to use for photobleach correction.
%           'None' - No correction performed.
%           MovieData - The photobleach correction from the input moviedata
%                       (e.g. a different ROI or cell movieData from the same
%                       experiment)
%           'Self' - The existing photobleach correction from the movie
%                    being processed is used. (Must be run previously. This
%                    uses the entire masked intensity)
%           'SelfCortical' - This uses the cortical intensity samples
%                   themselves for the fit. (Fit and correction is performed
%                   here, doesn't require having already been run. Only
%                   uses cortex)
%
%       ('TrendRemoval'->'None','Linear') Determines what trend removal is
%       performed on time-series before cross correlation analysis. Default
%       is Linear.
%
%       ('ChannelIndex' -> positive integer) Integer index of the
%       channel to analyze intensities from.
%       Optional. Default is all channels.
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the analysis to.
%       Optional. Default is a sub-directory of the movie's outputDirectory
%       called "masked_intensity_analysis"
%
%       (BatchMode->true/false) If true, all graphical output is
%       suppressed, such as progress bars, figures etc..
% 
% Output:
% 
% The analysis results will be written to file in the
% location specified by the OutputDirectory parameter, including any
% figures created.
% 
% Hunter Elliott
% 1/2012
% 

%% ------------------- Parameters ---------------- %%

iProcChan = 1;%Hard-coded channels for non-channel specific mask and post-processing association. Yeah I know this is a stupid way to do it, but who fucking cares.
%iPostProcChan = 2;
fName = 'intensity_analysis';%File name for output
%Figure saving parameters
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTif = {'-r100',...% dpi = 100
        '-dtiff'};% use tiff format
    

         
%Intensity fields
[intTypes,intNames] = getIntTypeFields;
isIntTypeMean = find(~cellfun(@isempty,strfind(intTypes,'Mean')));        %To avoid making a zillion figures, get index of means...
iIntForPB = 1;%Intensities to use for p.b. correct
iStatForPB = 3;%Use median since dist often non normal.
iPBChan = 1;%We won't have multi-channel live-cell data anytime soon, so fuck you I'm hard-coding all this.

nIntTypes = numel(intTypes);        
        


nCurvBins = 100;%It's just for display so arbitrarily use 100 
         
         
pfStr = 'PerFrame'; %String for naming per-frame stats

pfStatNames = {'Mean','STD','Median'};
pfStatFuns = {@mean,@std,@median};
nPfStat = numel(pfStatFuns);    


%% -------------------- Input -------------- %%

if nargin < 1 || ~isa(movieData,'MovieData3D')
    error('You must input a valid MovieData3D object as the first input!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous processes from this function
iProc = movieData.getProcessIndex('MaskedIntensity3DProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(MaskedIntensity3DProcess(movieData));                                                                                                 
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
nFrames = movieData.nFrames_;
%imSize = [movieData.imSize_ movieData.nSlices_];

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

%Convert sample radius to pixels
sampRadPix = p.CurvSampRad / pixXY;

%Branch profile fields and functions    
[curvTypes,curvNames,curvUnits,curvConv] = getCurveTypeFields(pixXY,false);
         
nCurvTypes = numel(curvTypes); 


iPruneProc = movieData.getProcessIndex('SkeletonPruningProcess',1,~p.BatchMode);
if isempty(iPruneProc) || ~movieData.processes_{iPruneProc}.checkChannelOutput(iProcChan)        
    error('No valid branch-pruning was found for this movie! Please run branch pruning first!')
end

%Make sure the mask geometry analysis has been run.
iMgProc = movieData.getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
if isempty(iMgProc) || ~movieData.processes_{iMgProc}.checkChannelOutput(iProcChan)
    error('No valid mask geometry analysis found! Please run mask geometry analysis first!')    
elseif movieData.processes_{iMgProc}.funParams_.PhysicalUnits
    %Lazy way to avoid un-scaling all the properties here.
    error('Please run mask geometry analysis with the physical units option set to false!')
end

outDir = p.OutputDirectory;
mkClrDir(outDir);

disp('Starting masked intensity analysis...')


%Check if this movie is an ROI
if ~isempty(movieData.parent_)
    roiInf = load(movieData.roiMaskPath_);                
else
    roiInf = [];
end


if nFrames <= 3
    p.PhotoBleachMeth = 'None';
end

switch p.PhotoBleachMeth
    
    case isa(p.PhotoBleachMeth,'MovieData3D')
    
        %If a photbleach correction ROI was input, load the fit data
        iPBProc = p.PhotoBleachMeth.getProcessIndex('PhotoBleachCalcProcess3D',1,~p.BatchMode);
        pbStat = p.PhotoBleachMeth.processes_{iPBProc}.loadChannelOutput(iProcChan);
        
    case 'Self'
        
        iPBProc = movieData.getProcessIndex('PhotoBleachCalcProcess3D',1,~p.BatchMode);
        pbStat = movieData.processes_{iPBProc}.loadChannelOutput(iProcChan);
        
    case 'SelfCortical'
        
        %This can only be calculated after we have the cortical intensity
        %samples, so we do nothing now.
        pbStat = [];
        
    otherwise        
        pbStat = [];
end

timePoints = 0:movieData.timeInterval_:(movieData.timeInterval_*(movieData.nFrames_-1));
timeUnits = 'Seconds';

%Use rule-of-thumb based on trajectory length to determine delays range for
%crosscorr
maxLag = floor(nFrames/4);
timeLags = -(movieData.timeInterval_*maxLag):movieData.timeInterval_:(movieData.timeInterval_*maxLag);

%TEMP - initialize per-frame mean stuctures, branch profile structures etc
%etc bla bla bla fuck i can't wait to be done with this


if p.BatchMode
    figArgs = {'Visible','off'};
else
    figArgs = {};
end    


%% --------------- Processing ----------------- %%

for iFrame = 1:nFrames
    
    %Load the mask and post processed pixels and make the voxels symmetric to get the stretched size
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
    
    if ~isempty(roiInf)
        %Keep only the ROI, so that the distance transforms will take into
        %account the ROI boundary.
        currMask = currMask(roiInf.cropX(1):roiInf.cropX(2),roiInf.cropY(1):roiInf.cropY(2),roiInf.cropZ(1):roiInf.cropZ(2));    
    end
    
    %currPP = tif3Dread([ppDir filesep ppNames{1}{iFrame}]);
    currMask = make3DImageVoxelsSymmetric(currMask,pixXY,pixZ);
    %currPP = make3DImageVoxelsSymmetric(currPP,pixXY,pixZ);            
    
    for iChan = 1:nChan    
        if iChan == 1
            currIm = zeros(size(currMask),'uint16');
        end
        
        tmp = stackRead([imDir{iChan} filesep imNames{iChan}{iFrame}]);        
        if ~isempty(roiInf)            
            currIm(:,:,:,iChan) = make3DImageVoxelsSymmetric(tmp(roiInf.cropX(1):roiInf.cropX(2),roiInf.cropY(1):roiInf.cropY(2),roiInf.cropZ(1):roiInf.cropZ(2)),pixXY,pixZ);
        else
            currIm(:,:,:,iChan) = make3DImageVoxelsSymmetric(tmp,pixXY,pixZ);
        end
    end         
        
    currMaskProp = movieData.processes_{iMgProc}.loadChannelOutput(iProcChan,iFrame);
    currSkelGraph = movieData.processes_{iPruneProc}.loadChannelOutput(iProcChan,iFrame);
    
    if ~isempty(roiInf)
        %Adjust the skeleton and mask geometry to agree with the cropped
        %image & mask   
        currMaskProp = adjustROIMaskGeometry(currMaskProp,roiInf,pixXY,pixZ);
        currSkelGraph = adjustROISkeletonGraph(currSkelGraph,roiInf,pixXY,pixZ);
    end        
        
    branchProfiles(iFrame) = analyze3DImageMaskedIntensities(currIm,currMask,currSkelGraph,currMaskProp,sampRadPix);
    
    %Calculate per-frame averages, apply photobleach correction if
    %necessary
    if ~isempty(pbStat)
        currPBCorr = pbStat.fitValues(iFrame);
    else
        currPBCorr = 1;
    end        
    %TEMP- DOESN"T WORK FOR MULTICHANNEL YOU IDIOT!!!!
    for k = 1:nPfStat
        for j = 1:nCurvTypes
            curvStats.([pfStr pfStatNames{k} curvTypes{j}])(iFrame,iChan) = pfStatFuns{k}(branchProfiles(iFrame).(curvTypes{j}));
            
        end
        for j = 1:nIntTypes
            intStats.([pfStr pfStatNames{k} intTypes{j}])(iFrame,iChan) = pfStatFuns{k}(branchProfiles(iFrame).(intTypes{j})(:,iChan)) / currPBCorr;
        end
    end
end

outVars = {'branchProfiles','intStats','curvStats'};

%% ---- Apply Sef-cortical correction if requested --- %%


if strcmp(p.PhotoBleachMeth,'SelfCortical')
    
    
    
    % ---- Calc the FIt ----- %
    
    
    
    %TEMP - THIS IS STUPID, AND A DUPLICATION OF CODE BUT I JUST GIVE ZERO
    %FUCKS RIGHT NOW. move to its own function at some point 
    fitOptions = statset('Robust','off','MaxIter',500,'Display','off');
    fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));     %Double-exponential function for fitting
    intTsforPB = intStats.([pfStr pfStatNames{iStatForPB} intTypes{iIntForPB}])';
    bInitCurr([1 3]) = intTsforPB(1)/2;%Aplitudes based on first frame
    bInitCurr([2 4]) = log(intTsforPB(end) / (2*bInitCurr(1)))/timePoints(end);
    
    [pbStat.bFit(:,iChan),pbStat.resFit(:,iChan),pbStat.jacFit(:,:,iChan),pbStat.covFit(:,:,iChan),pbStat.mseFit(iChan)] ...
        = nlinfit(timePoints(:),intTsforPB(:),fitFun,bInitCurr,fitOptions);
    %Get confidence intervals of fit and fit values
    [pbStat.fitValues(:,iChan),pbStat.deltaFit(:,iChan)] = nlpredci(...
        fitFun,timePoints(:),pbStat.bFit(:,iChan),pbStat.resFit(:,iChan),'covar',pbStat.covFit(iChan),'mse',pbStat.mseFit(iChan));
    
    
    
    % ---- Apply the FIt ----- %
    
    for iFrame = 1:nFrames
        for k = 1:nPfStat            
            for j = 1:nIntTypes
                intStats.([pfStr pfStatNames{k} intTypes{j}])(iFrame) =  pfStatFuns{k}(branchProfiles(iFrame).(intTypes{j})) / pbStat.fitValues(iFrame);
            end
        end                
    end
    
    
    
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
    ylabel(intTypes{iIntForPB})
    plot(timePoints,intTsforPB,'.-')
    plot(timePoints,pbStat.fitValues(:,iChan),'r')
    yl = ylim;
    plot(timePoints,pbStat.fitValues(:,iChan)+pbStat.deltaFit(iChan),'--r')
    legend(intTypes{iIntForPB},'Fit','Fit 95% C.I.')
    plot(timePoints,pbStat.fitValues(:,iChan)-pbStat.deltaFit(iChan),'--r')
    ylim(yl);
    xlim(timePoints([1 end]))
    
    subplot(2,1,2)
    hold on
    plot(timePoints,pbStat.resFit,'.-')
    if ~isempty(movieData.timeInterval_)
        xlabel('Time, seconds')
    else
        xlabel('Frame Number')
    end
    ylabel('Fit Residuals')
    xlim(timePoints([1 end]))
    plot(xlim,[0 0],'--k')
    
    
    
    currFigFile = [p.OutputDirectory filesep 'self-cortical photobleach correct'];
    hgsave(fitFig,currFigFile)
    print(fitFig,currFigFile,pOptTif{:})
    print(fitFig,currFigFile,pOpt{:})
    
    %Log this file name in the parameter structure    

    if p.BatchMode && ishandle(fitFig) %make sure user hasn't closed it.
        close(fitFig)
    end
end
outVars = [outVars{:} {'pbStat'}];

%% --------------- Figure Creation ---------------- %%

% --- Init --- %


chanNames = cellfun(@(x)(x(max(strfind(x,filesep))+1:end)),imDir,'Unif',0);

%We assume the channels are in order of wavelength, so reverse the colors
if nChan <= 3
    chanCols = [0 0 1;
                0 1 0;
                1 0 0];
else
    chanCols = jet(nChan);
    chanCols = chanCols(end:-1:1,:);
end
endVal = .95;
chanFrameCols = zeros(nChan,3,nFrames);

if nFrames > 1
    for j = 1:nChan
        for k = 1:3
            chanFrameCols(j,k,:) = linspace(chanCols(j,k),endVal,nFrames);
        end
    end
else
    chanFrameCols = chanCols;
end

%For single-channel figures
frameCols = jet(nFrames);


%% ---- Avg. Intensity vs. distance line overlay plot ----%


% ----- Int Vs. Dist, Normalized to first frame ---- %

dLabel = {'Distance from membrane [nm]','Negative is outside cell, positive inside'};
figure(figArgs{:});
hold on
for iFrame = 1:nFrames               
    for k = 1:nChan
        if iFrame == 1
            %Normalize to first frame
            chanMin(k) = nanmin(branchProfiles(iFrame).wholeMaskMeanVsDist(k,:));                
            chanMax(k) = nanmax(branchProfiles(iFrame).wholeMaskMeanVsDist(k,:));
        end
        tmp = branchProfiles(iFrame).wholeMaskMeanVsDist(k,:) - chanMin(k);
        tmp = tmp ./ (chanMax(k)-chanMin(k));
        plot(branchProfiles(iFrame).wholeMaskDists * movieData.pixelSize_,tmp,'.-','color',chanFrameCols(k,:,iFrame),'LineWidth',2,'MarkerSize',10);    
    end    
    if iFrame == 1        
        legend(chanNames,'Location','NorthWest');
    end
end


% %Transparency crashes orchestra, so I disabled this for now.
% %SOOOO LAZY AND TIRED RIGHT NOW - DO IT THIS WAY SO LEGEND MAKES SENSE
% for k = 1:nChan
%     tmp = branchProfiles(iFrame).wholeMaskMeanVsDist(k,:) - nanmin(branchProfiles(iFrame).wholeMaskMeanVsDist(k,:));            
%     tmp2 = branchProfiles(iFrame).wholeMaskSTDVsDist(k,:) / nanmax(tmp);
%     tmp = tmp ./ nanmax(tmp);        
%     plotTransparent(branchProfiles(iFrame).wholeMaskDists * movieData.pixelSize_,tmp,tmp2,chanCols(k,:),.3,0);        
% end        

xlabel(dLabel)
ylabel('Normalized Average Fluoresence, a.u.')
title({'Distance - Intensity Profile, Whole-Cell','Bands show +/- STD','Intensities normalized to first frame'})
plot([0 0],ylim,'--k')
figName = 'distance vs masked intensity overlay';
set(gca,'Color',[.7 .7 .7])
print(pOpt{:},[outDir filesep figName '.eps']);
print(pOptTif{:},[outDir filesep figName '.tif']);
hgsave(gcf,[outDir filesep figName '.fig']);

% ----- Int Vs. Dist, Normalized Each frame ---- %

dLabel = {'Distance from membrane [nm]','Negative is outside cell, positive inside'};
figure(figArgs{:});
hold on
for iFrame = 1:nFrames               
    for k = 1:nChan        
        %Normalize each frame
        chanMin(k) = nanmin(branchProfiles(iFrame).wholeMaskMeanVsDist(k,:));                
        chanMax(k) = nanmax(branchProfiles(iFrame).wholeMaskMeanVsDist(k,:));        
        tmp = branchProfiles(iFrame).wholeMaskMeanVsDist(k,:) - chanMin(k);
        tmp = tmp ./ (chanMax(k)-chanMin(k));
        plot(branchProfiles(iFrame).wholeMaskDists * movieData.pixelSize_,tmp,'.-','color',chanFrameCols(k,:,iFrame),'LineWidth',2,'MarkerSize',10);    
    end    
    if iFrame == 1        
        legend(chanNames,'Location','NorthWest');
    end
end


% %Transparency crashes orchestra, so I disabled this for now.
% %SOOOO LAZY AND TIRED RIGHT NOW - DO IT THIS WAY SO LEGEND MAKES SENSE
% for k = 1:nChan
%     tmp = branchProfiles(iFrame).wholeMaskMeanVsDist(k,:) - nanmin(branchProfiles(iFrame).wholeMaskMeanVsDist(k,:));            
%     tmp2 = branchProfiles(iFrame).wholeMaskSTDVsDist(k,:) / nanmax(tmp);
%     tmp = tmp ./ nanmax(tmp);        
%     plotTransparent(branchProfiles(iFrame).wholeMaskDists * movieData.pixelSize_,tmp,tmp2,chanCols(k,:),.3,0);        
% end        

xlabel(dLabel)
ylabel('Normalized Average Fluoresence, a.u.')
title({'Distance - Intensity Profile, Whole-Cell','Bands show +/- STD','Intensities normalized Each frame. Darker is first frame, lighter is last'})
plot([0 0],ylim,'--k')
figName = 'distance vs masked intensity overlay';
set(gca,'Color',[.7 .7 .7])
print(pOpt{:},[outDir filesep figName '.eps']);
print(pOptTif{:},[outDir filesep figName '.tif']);
hgsave(gcf,[outDir filesep figName '.fig']);


%% ----2D Intensity histogram vs. Distance figure -----%

iFrame = 1;%Only doing first frame for this one because no one cares about this figure
figure(figArgs{:})
for k = 1:nChan
    subplot(1,nChan,k)
    hold on
    imagesc(branchProfiles(iFrame).wholeMaskDists * movieData.pixelSize_,branchProfiles(iFrame).wholeMaskIntHistBins(k,:),squeeze(branchProfiles(iFrame).wholeMaskNormHistVsDist(k,:,:))')
    ylim(branchProfiles(iFrame).wholeMaskIntHistBins(k,[1 end]))%For some reason the limits sometimes exceed these ranges, so force them.
    xlim(branchProfiles(iFrame).wholeMaskDists([1 end]) * movieData.pixelSize_)
    caxis([0 .1])
    plot([0 0 ],ylim,'w--')
    xlabel(dLabel)
    ylabel('Intensity, A.U.')
    title({'Normalized whole-mask int vs. dist histograms',chanNames{k} ,'FIRST FRAME ONLY'})    
end
print(pOpt{:},[outDir filesep 'distance vs masked intensity histogram.eps']);
print(pOptTif{:},[outDir filesep 'distance vs masked intensity histogram.tif']);
hgsave(gcf,[outDir filesep 'distance vs masked intensity histogram.fig']);

%% ----------- Per-Branch Intensity Figures ------ %%


avgBranchInt = cell(nChan,1);
for k = 1:nChan    
    f1 = figure(figArgs{:});
    hold on
    f2 = figure(figArgs{:});
    hold on
    
    
    for iFrame = 1:nFrames
        
        figure(f1);
                
        %Quick workaround for failure in per-branch intensity analysis... Shitty
        %but hey it works.
        hasLvD = ~cellfun(@isempty,branchProfiles(iFrame).branchTipPixelLenVsDepth);
        branchProfiles(iFrame).branchTipPixelLenVsDepth(~hasLvD) = deal({NaN(1,2)});
        avgBranchRadius = cellfun(@(x)(mean(x(:,2))),branchProfiles(iFrame).branchTipPixelLenVsDepth) .* movieData.pixelSize_;
        branchLen = cellfun(@(x)(max(x(:,1))),branchProfiles(iFrame).branchTipPixelLenVsDepth) .* movieData.pixelSize_;

        avgBranchInt{k} = cellfun(@(x)(mean(x(:,k))),branchProfiles(iFrame).branchTipPixelInt);        
        plot(avgBranchRadius,avgBranchInt{k},'o','MarkerSize',15,'color',chanCols(k,:));    
        
        figure(f2)
        plot(branchLen,avgBranchInt{k},'o','MarkerSize',15,'color',chanCols(k,:));    
        
        
    end
    figure(f1)
    title('Intensity vs. Radius, All Frames')
    ylabel(['Mean ' chanNames{k} ' intensity, a.u.'])
    xlabel('Mean branch radius, nm')   
    
    figName = ['branch radius vs masked intensity in ' chanNames{k}];
    print(pOpt{:},[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);
    
    figure(f2)
    
    title('Intensity vs. Approximate Length, All Frames')
    ylabel(['Mean ' chanNames{k} ' intensity, a.u.'])
    xlabel('Approximate Branch Length, nm')   

    figName = ['branch length vs masked intensity in ' chanNames{k}];
    print(pOpt{:}   ,[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);
               
end


%% ----------- Intensity vs. Curvature Figures ------------- %


%TEMP - convert to physical units!!!


iFrame = 1;%TEEEMMMPPP

for l = 1:nCurvTypes
    
    currCurv = real(branchProfiles(iFrame).(curvTypes{l}));
    curvBins = linspace(prctile(currCurv,.5),prctile(currCurv,99.5),nCurvBins);
    curvBins = sort(curvBins);%For the strictly negative measures, we need to make sure these are increasing
    
    for k = 1:nChan                                       
        
        
        %currInt = branchProfiles(iFrame).intForCurvSampMean(:,k);
        currMaxInt = branchProfiles(iFrame).intForCurvSampMax(:,k);
        currMaxIntDN = branchProfiles(iFrame).intForCurvSampMaxDepthNorm(:,k);                        
            
        % --- Intensity - vs - Curvature PDF -- %
        
        figure(figArgs{:})
        subplot(2,1,1)
        hold on
        N = hist3([currMaxInt currCurv],{branchProfiles(iFrame).wholeMaskIntHistBins(k,:), curvBins});
        for m = 1:size(N,1)
            N(m,:) = N(m,:) ./ sum(N(m,:));
        end       
        imagesc(branchProfiles(iFrame).wholeMaskIntHistBins(k,:),curvBins,log10(N'));
        colorbar
        %caxis([0 prctile(N(:),95)])
        xlim(branchProfiles(iFrame).wholeMaskIntHistBins(k,[1 end]))
        ylim(curvBins([1 end]))
        xlabel('Intensity, a.u.')
        ylabel(curvNames{l})
        title({['2D PDF of ' curvNames{l} ' vs. Max Local intensity in ' chanNames{k}],...
                'Color indicates log10(probability)','FIRST FRAME ONLY'})        
                        
        % --- Depth-Normalized Intensity - vs - Curvature PDF -- %
        
        subplot(2,1,2)
        hold on
        N = hist3([currMaxIntDN currCurv],{branchProfiles(iFrame).wholeMaskIntHistBins(k,:), curvBins});
        for m = 1:size(N,1)
            N(m,:) = N(m,:) ./ sum(N(m,:));
        end       
        imagesc(branchProfiles(iFrame).wholeMaskIntHistBins(k,:),curvBins,log10(N'));
        colorbar
        %caxis([0 prctile(N(:),95)])
        xlabel('Depth-Normalized Intensity, a.u.')
        ylabel(curvNames{l})
        xlim(branchProfiles(iFrame).wholeMaskIntHistBins(k,[1 end]))
        ylim(curvBins([1 end]))
        title({['2D PDF of ' curvNames{l} ' vs. Depth-Normalized Max Local Intensity in Channel ' chanNames{k}],...            
            'Color indicates probability','FIRST FRAME ONLY'})   
        figName = ['max local intensity in ' chanNames{k} ' vs ' curvNames{l} ' histogram and PDF '];
        print(pOpt{:},[outDir filesep figName '.eps']);
        print(pOptTif{:},[outDir filesep figName '.tif']);
        hgsave(gcf,[outDir filesep figName '.fig']);

        nBinsCur = numel(branchProfiles(iFrame).wholeMaskIntHistBins(k,:));
        avgCurvVsInt = nan(1,nBinsCur-1);
        stdCurvVsInt = nan(1,nBinsCur-1);
        semCurvVsInt = nan(1,nBinsCur-1);
        nCurvVsInt = nan(1,nBinsCur-1);
        avgCurvVsIntDN = nan(1,nBinsCur-1);
        stdCurvVsIntDN = nan(1,nBinsCur-1);
        semCurvVsIntDN = nan(1,nBinsCur-1);
        nCurvVsIntDN = nan(1,nBinsCur-1);
        for i = 1:nBinsCur-1
            currCurrVals = currMaxInt >= branchProfiles(iFrame).wholeMaskIntHistBins(k,i) & currMaxInt < branchProfiles(iFrame).wholeMaskIntHistBins(k,i+1);
            avgCurvVsInt(i) = mean(currCurv(currCurrVals));
            stdCurvVsInt(i) = std(currCurv(currCurrVals));
            nCurvVsInt(i) = nnz(currCurrVals);
            semCurvVsInt(i) = stdCurvVsInt(i) / sqrt(nCurvVsInt(i));            
            
            currCurrVals = currMaxIntDN >= branchProfiles(iFrame).wholeMaskIntHistBins(k,i) & currMaxIntDN < branchProfiles(iFrame).wholeMaskIntHistBins(k,i+1);
            avgCurvVsIntDN(i) = mean(currCurv(currCurrVals));
            stdCurvVsIntDN(i) = std(currCurv(currCurrVals));
            nCurvVsIntDN(i) = nnz(currCurrVals);
            semCurvVsIntDN(i) = stdCurvVsIntDN(i) / sqrt(nCurvVsIntDN(i));            
        end

        figure(figArgs{:});
        subplot(2,1,1)
        hold on
        title({['Average ' curvNames{l} ' vs. Max Local Intensity in Channel ' chanNames{k}],'FIRST FRAME ONLY'})                    
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),avgCurvVsInt)
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),avgCurvVsIntDN,'r')
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),avgCurvVsInt + 1.96*semCurvVsInt,'--')
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),avgCurvVsInt - 1.96*semCurvVsInt,'--')
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),avgCurvVsIntDN + 1.96*semCurvVsIntDN,'--r')
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),avgCurvVsIntDN - 1.96*semCurvVsIntDN,'--r')
        xlabel('Intensity, a.u.')
        ylabel(['Average ' curvNames{l}])
        legend('Average Intensity','Average Depth-Normalized Intensiyt','+/- 1.96SEM')

        subplot(2,1,2)
        hold on
        title({['Local STD of ' curvNames{l} ' vs. Max Local Intensity in Channel ' chanNames{k}],'FIRST FRAME ONLY'})                         
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),stdCurvVsInt)
        plot(branchProfiles(iFrame).wholeMaskIntHistBins(k,1:end-1),stdCurvVsIntDN)
        xlabel('Intensity, a.u.')
        ylabel(['STD of ' curvNames{l}])                                  
        legend('Average Intensity','Average Depth-Normalized Intensiyt')

        figName = ['average and STD of ' curvNames{l} ' vs max local intensity in ' chanNames{k}];
        print(pOpt{:},[outDir filesep figName '.eps']);
        print(pOptTif{:},[outDir filesep figName '.tif']);
        hgsave(gcf,[outDir filesep figName '.fig']);

    end
    
end

%% ------ Intensity and Curvature over time figures ------ %%


%% ---- Curv and Int over Time --- %%

if nFrames > 1

    %Only do mean intensity to avoid making a million figures
    ccIntCurv = nan(1,nCurvTypes,maxLag*2+1);
    for iIntType = 1

        for iCurvType = 1:nCurvTypes

            for iChan = 1:nChan


                % --- CUrv and INt dist over time --- %

                figure(figArgs{:})            
                hold on
                for iFrame = 1:nFrames
                    subplot(2,1,1)
                    hold on
                    [intDist,xPos] = ksdensity(branchProfiles(iFrame).(intTypes{iIntType})(:,iChan));
                    plot(xPos,intDist,'color',frameCols(iFrame,:))

                    subplot(2,1,2)
                    hold on
                    [intDist,xPos] = ksdensity(real(branchProfiles(iFrame).(curvTypes{iCurvType})));%Real is for k1 and k2
                    plot(xPos,intDist,'color',frameCols(iFrame,:))


                end
                subplot(2,1,1)    
                title({'Intensity and Curvature Distributions over time','Blue = First frame, Red = Last'})
                xlabel([intNames{iIntType} ', Channel ' chanNames{iChan} ', a.u.'])
                ylabel('Probability')
                subplot(2,1,2)       
                xlabel(curvNames{iCurvType})
                ylabel('Probability')

                figName = ['Dist over Time ' curvNames{iCurvType} ' vs ' intNames{iIntType} ' in channel ' chanNames{k}];
                print(pOpt{:},[outDir filesep figName '.eps']);
                print(pOptTif{:},[outDir filesep figName '.tif']);
                hgsave(gcf,[outDir filesep figName '.fig']);


                % --- Curv and int stats over time --- %


                for iStatType = 1:nPfStat

                    intTs = intStats.([pfStr pfStatNames{iStatType} intTypes{iIntType}])(:);
                    curvTs = real(curvStats.([pfStr pfStatNames{iStatType} curvTypes{iCurvType}])(:));

                    if strcmp(p.TrendRemoval,'Linear') && nFrames > 2
                        intTs = removeLinearTrend(intTs);
                        curvTs = removeLinearTrend(curvTs);                    
                    end

                    figure(figArgs{:})
                    [ax,h1,h2] = plotyy(timePoints,intTs,...
                                        timePoints,curvTs);
                    xlabel(['Time, ' timeUnits])
                    set(get(ax(1),'Ylabel'),'String',[pfStatNames{iStatType} ' ' intNames{iIntType}])
                    set(get(ax(2),'Ylabel'),'String',[pfStatNames{iStatType} ' ' curvNames{iCurvType}])
                    title(['Curvature and Intensity ' pfStatNames{iStatType} ' over Time'])
                    figName = [pfStatNames{iStatType} ' over Time of ' curvNames{iCurvType} ' and ' intNames{iIntType} ' in channel ' chanNames{k}];
                    print(pOpt{:},[outDir filesep figName '.eps']);
                    print(pOptTif{:},[outDir filesep figName '.tif']);
                    hgsave(gcf,[outDir filesep figName '.fig']);


                    tmp = crossCorr(intTs,curvTs,maxLag);
                    ccIntCurv(iIntType,iCurvType,:) = tmp(:,1);

                    figure(figArgs{:})
                    hold on
                    plot(timeLags,tmp(:,1))
                    xlim([min(timeLags) max(timeLags)])
                    ylim([-1 1])
                    plot(xlim,[1 1]*1.96/sqrt(nFrames),'--r')
                    legend('Correlation','95% Signficance bounds')                
                    plot(xlim,[1 1]*-1.96/sqrt(nFrames),'--r')                
                    plot(xlim,[0 0 ],'--k')
                    plot([0 0], ylim,'--k')
                    title({['Cross-Corr ' pfStatNames{iStatType} ' of ' intNames{iIntType} ' and ' curvNames{iCurvType} ' over Time'],...
                           'Positive Delay means curvature follows intensity'})
                    xlabel(['Delay, ' timeUnits])
                    ylabel('Correlation Coefficient')
                    figName = ['CrossCorr ' pfStatNames{iStatType} ' of ' intNames{iIntType} ' and ' curvNames{iCurvType} ' in channel ' chanNames{k}];
                    print(pOpt{:},[outDir filesep figName '.eps']);
                    print(pOptTif{:},[outDir filesep figName '.tif']);
                    hgsave(gcf,[outDir filesep figName '.fig']);


                end    



            end
        end        
    end


    outVars = [outVars{:} {'ccIntCurv','intNames','curvNames','pfStatNames','intTypes','curvTypes','chanNames'}];

else
    outVars = [outVars{:} {'intNames','curvNames','pfStatNames','intTypes','curvTypes','chanNames'}];
    
end


%% ---------------- Finalize ------------------ %%

if p.BatchMode
    close all
end



%Save analysis to disk, store directories and parameters in process object
outFile = [outDir filesep fName '.mat'];
save(outFile,outVars{:});
movieData.processes_{iProc}.setOutFilePath(iProcChan,outFile);
movieData.processes_{iProc}.setDateTime;
movieData.save;

disp('Finished analyzing masked intensities!')









