function analyze3DMovieMaskedIntensities(movieData,varargin)
%ANALYZE3DMOVIEMASKEDINTENSITIES calls analyze3DImageMaskedIntensities on each frame of the input movie 
% 
% analyze3DMovieMaskedIntensities(movieData)
% analyze3DMovieMaskedIntensities(movieData,paramsIn)
% analyze3DMovieMaskedIntensities(movieData,'ParamName1',paramValue1,...)
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
%       ('ChannelIndex' -> positive integer) Integer index of the
%       channel to analyze intensities from.
%       Optional. Default is all channels.
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the analysis to.
%       Optional. Default is a sub-directory of the movie's outputDirectory
%       called "intensity_analysis"
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


%% -------------------- Input -------------- %%

if nargin < 1 || ~isa(movieData,'MovieData3D')
    error('You must input a valid MovieData3D object as the first input!')
end

%FINISH THIS SHIT YOU LAZY BASTARD!

%TEMP OOOOOOOBBBBBBBVVVVVIIIOOOUUUSSSSSLLLLYYYYY
p.ChannelIndex = 1:3;
p.BatchMode = true;
iProcChan = 1;%Hard-coded channels for non-channel specific mask and post-processing association. Yeah I know this is a stupid way to do it, but who fucking cares.
iPostProcChan = 2;

%% ------------------- Init --------------- %%

%Make sure the movie has been segmented and has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(iProcChan)%TEMP? how do deal with three-input channels for mask and one output?? ALways associate with channel 1?
    error('The input MovieData does not have valid masks for the selected channel!');
end

%Get mask file names and directory
maskDir = movieData.processes_{iSegProc}.outFilePaths_{1,iProcChan};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(iProcChan);
%And the post-process added pixel directory and file-names
ppDir = movieData.processes_{iSegProc}.outFilePaths_{1,iPostProcChan};
ppNames = movieData.processes_{iSegProc}.getOutMaskFileNames(iPostProcChan);

%Get image locations and info
imDir = movieData.getChannelPaths(p.ChannelIndex);
imNames = movieData.getImageFileNames(p.ChannelIndex);
nChan = numel(p.ChannelIndex);
nFrames = movieData.nFrames_;
imSize = [movieData.imSize_ movieData.nSlices_];
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

iPruneProc = movieData.getProcessIndex('SkeletonPruningProcess',1,1);
if isempty(iPruneProc) || ~movieData.processes_{iPruneProc}.checkChannelOutput(iProcChan)        
    error('No valid branch-pruning was found for this movie! Please run branch pruning first!')
end

%Make sure the mask geometry analysis has been run.
iMgProc = movieData.getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
if isempty(iMgProc) || ~movieData.processes_{iMgProc}.checkChannelOutput(iProcChan)
    error('No valid mask geometry analysis found! Please run mask geometry analysis first!')    
elseif movieData.processes_{iMgProc}.funParams_.PhysicalUnits
    %Lazy way to avoid re-scaling all the properties here.
    error('Please run mask geometry analysis with the physical units option set to false!')
end


%% --------------- Processing ----------------- %%

for iFrame = 1:nFrames
    
    %Load the mask and post processed pixels and make the voxels symmetric to get the stretched size
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
    currPP = tif3Dread([ppDir filesep ppNames{1}{iFrame}]);
    currMask = make3DImageVoxelsSymmetric(currMask,pixXY,pixZ);
    currPP = make3DImageVoxelsSymmetric(currPP,pixXY,pixZ);
    
    for iChan = 1:nChan    
        if iChan == 1
            currIm = zeros(size(currMask),'uint16');
        end
        currIm(:,:,:,iChan) = make3DImageVoxelsSymmetric(stackRead([imDir{iChan} filesep imNames{iChan}{iFrame}]),pixXY,pixZ);                                    
    end         
    
    currMaskProp = movieData.processes_{iMgProc}.loadChannelOutput(iProcChan,iFrame);
    currSkelGraph = movieData.processes_{iPruneProc}.loadChannelOutput(iProcChan,iFrame);
    
    branchProfiles(iFrame) = analyze3DImageMaskedIntensities(currIm,currMask,currSkelGraph,currMaskProp);
    
    
end


%TEEEEEEMMMPPPP!!!
outDir = [movieData.outputDirectory_ filesep 'intensity analysis'];
mkdir(outDir);
save([outDir filesep 'intensity analysis.mat'],'branchProfiles','currMask','currIm','currMaskProp','currSkelGraph','movieData');

%% --------------- Figure Creation ---------------- %%

%TEMP - stupid figure creation problem las tminute fix
p.ChannelIndex = 1:3;
p.BatchMode = false;
iProcChan = 1;
nChan = size(currIm,4);
outDir = [movieData.outputDirectory_ filesep 'intensity analysis'];

pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTif = {'-r300',...% dpi = 300
        '-dtiff'};% use eps format
    
    
%dLabel = {'Distance from membrane [pixels]','Negative is outside cell, positive inside'};
dLabel = {'Distance from membrane [nm]','Negative is outside cell, positive inside'};
branchProfiles.wholeMaskDists = branchProfiles.wholeMaskDists .* movieData.pixelSize_;
plotChans = [1 3];%TEMP!!!!
chanNames = {'EGFP-Myosin','Membrane Dye','Phalloidin-Actin'};

%----- Overall intensity correlation histogram ----- %%

myoIm = currIm(:,:,:,1);
dyeIm = currIm(:,:,:,2);
actIm = currIm(:,:,:,3);
[N,C] = hist3([myoIm(currMask(:)) actIm(currMask(:))],[256 256]);


if p.BatchMode
    figArgs = {'Visible','off'};
else
    figArgs = {};
end


figure(figArgs{:})
%TEEEEEMMMPpp!!TEMP
imagesc(C{2},C{1},N);axis xy
caxis([0 10])
xlabel(chanNames{plotChans(2)})
ylabel(chanNames{plotChans(1)})

print(pOpt{:},[outDir filesep 'channel masked intensity bivariate histogram.eps']);
print(pOptTif{:},[outDir filesep 'channel masked intensity bivariate histogram.tif']);
hgsave(gcf,[outDir filesep 'channel masked intensity bivariate histogram.fig']);

%---- Avg. Intensity vs. distance line overlay plot ----%


figure(figArgs{:})
hold on
chanCols = [0 0 1; %TEMP!
            0 1 0;
            1 0 0];
        
        
chanStr = cell(nChan,1);
for k = plotChans
    tmp = branchProfiles.wholeMaskMeanVsDist(k,:) - nanmin(branchProfiles.wholeMaskMeanVsDist(k,:));        
    tmp2 = branchProfiles.wholeMaskSTDVsDist(k,:) / nanmax(tmp);
    tmp = tmp ./ nanmax(tmp);
    plot(branchProfiles.wholeMaskDists,tmp,'.-','color',chanCols(k,:),'LineWidth',2,'MarkerSize',10);
    chanStr{k} = ['Channel ' num2str(k)];
end    
legend(chanNames{plotChans});

%SOOOO LAZY AND TIRED RIGHT NOW - DO IT THIS WAY SO LEGEND MAKES SENSE
for k = plotChans
    tmp = branchProfiles.wholeMaskMeanVsDist(k,:) - nanmin(branchProfiles.wholeMaskMeanVsDist(k,:));        
    %tmp2 = branchProfiles.wholeMaskSTDVsDist(k,:) / nanmax(tmp) ./ sqrt(branchProfiles.wholeMaskNumPixVsDist(:) * 1.96);
    tmp2 = branchProfiles.wholeMaskSTDVsDist(k,:) / nanmax(tmp);
    tmp = tmp ./ nanmax(tmp);        
    plotTransparent(branchProfiles.wholeMaskDists,tmp,tmp2,chanCols(k,:),.3,0);        
end        

xlabel(dLabel)
ylabel('Normalized Average Fluoresence')
%title({'Distance - Intensity Profile, Whole-Cell','Bands show +/- 1.96/sqrt(n)'})
title({'Distance - Intensity Profile, Whole-Cell','Bands show +/- STD'})
plot([0 0],ylim,'--k')

print(pOpt{:},[outDir filesep 'distance vs masked intensity overlay.eps']);
print(pOptTif{:},[outDir filesep 'distance vs masked intensity overlay.tif']);
hgsave(gcf,[outDir filesep 'distance vs masked intensity overlay.fig']);
%----2D Intensity histogram vs. Distance figure -----%

figure(figArgs{:})
for k = plotChans
    subplot(1,nChan,k)
    hold on
    imagesc(branchProfiles.wholeMaskDists,branchProfiles.wholeMaskIntHistBins(k,:),squeeze(branchProfiles.wholeMaskNormHistVsDist(k,:,:))')
    ylim(branchProfiles.wholeMaskIntHistBins(k,[1 end]))%For some reason the limits sometimes exceed these ranges, so force them.
    xlim(branchProfiles.wholeMaskDists([1 end]))
    caxis([0 .1])
    plot([0 0 ],ylim,'w--')
    xlabel(dLabel)
    ylabel('Intensity, A.U.')
    title({'Normalized whole-mask int vs. dist histograms',chanNames{k}})    
end
print(pOpt{:},[outDir filesep 'distance vs masked intensity histogram.eps']);
print(pOptTif{:},[outDir filesep 'distance vs masked intensity histogram.tif']);
hgsave(gcf,[outDir filesep 'distance vs masked intensity histogram.fig']);

% ----------- Per-Branch Intensity Figures ------ %%

avgBranchRadius = cellfun(@(x)(mean(x(:,2))),branchProfiles.branchTipPixelLenVsDepth) .* movieData.pixelSize_;
branchLen = cellfun(@(x)(max(x(:,1))),branchProfiles.branchTipPixelLenVsDepth) .* movieData.pixelSize_;

for k = plotChans    
    figure(figArgs{:});
    avgBranchInt{k} = cellfun(@(x)(mean(x(:,k))),branchProfiles.branchTipPixelInt);
    plot(avgBranchRadius,avgBranchInt{k},'o','MarkerSize',15,'color',chanCols(k,:));    
    ylabel(['Mean ' chanNames{k} ' intensity, a.u.'])
    xlabel('Mean branch radius, nm')   
    
    print(pOpt{:},[outDir filesep 'branch radius vs masked intensity overlay.eps']);
    print(pOptTif{:},[outDir filesep 'branch radius vs masked intensity overlay.tif']);
    hgsave(gcf,[outDir filesep 'branch radius vs masked intensity overlay.fig']);
    
    figure(figArgs{:});
    plot(branchLen,avgBranchInt{k},'o','MarkerSize',15,'color',chanCols(k,:));    
    ylabel(['Mean ' chanNames{k} ' intensity, a.u.'])
    xlabel('Approximate Branch Length, nm')   

    print(pOpt{:}   ,[outDir filesep 'branch length vs masked intensity overlay.eps']);
    print(pOptTif{:},[outDir filesep 'branch length vs masked intensity overlay.tif']);
    hgsave(gcf,[outDir filesep 'branch length vs masked intensity overlay.fig']);
end

if p.BatchMode
    close all
end












