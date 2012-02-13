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
iPostProcChan = 2;
fName = 'intensity_analysis';%File name for output
%Figure saving parameters
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTif = {'-r100',...% dpi = 100
        '-dtiff'};% use tiff format
    


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
%And the post-process added pixel directory and file-names
ppDir = movieData.processes_{iSegProc}.outFilePaths_{1,iPostProcChan};
ppNames = movieData.processes_{iSegProc}.getOutMaskFileNames(iPostProcChan);

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
mkdir(outDir);
%mkClrDir(outDir);

disp('Starting masked intensity analysis...')


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
    
    %TEMP - for testing for biases in sampling
    %currIm(:,:,:,end) = bwdist(~currMask);
    
    currMaskProp = movieData.processes_{iMgProc}.loadChannelOutput(iProcChan,iFrame);
    currSkelGraph = movieData.processes_{iPruneProc}.loadChannelOutput(iProcChan,iFrame);
    
    %branchProfiles(iFrame) = analyze3DImageMaskedIntensities(currIm,currMask,currSkelGraph,currMaskProp,currPP);
    
%TEMP TEMP TEEEEMMMPPP!!
branchProfiles = load('L:\nih\4D gfpMIIB fix and stain\temp distNorm branch profiles for debug.mat');
branchProfiles = branchProfiles.bp2Norm;
% 
% branchProfiles = load('L:\nih\4D gfpMIIB fix and stain\temp branch profiles for debug.mat');
% branchProfiles = branchProfiles.bp;

end


%% --------------- Figure Creation ---------------- %%

if p.BatchMode
    figArgs = {'Visible','off'};
else
    figArgs = {};
end    
chanNames = cellfun(@(x)(x(max(strfind(x,filesep))+1:end)),imDir,'Unif',0);

%TEMP - for testing for biases in sampling
%chanNames{end} = 'DistX';
chanNames{end} = 'Depth-Normalized 488';
chanNames{2} = 'Depth-Normalized 560';

%We assume the channels are in order of wavelength, so reverse the colors
if nChan <= 3
    chanCols = [0 0 1;
                0 1 0;
                1 0 0];
else
    chanCols = jet(nChan);
    chanCols = chanCols(end:-1:1,:);
end
%---- Avg. Intensity vs. distance line overlay plot ----%

dLabel = {'Distance from membrane [nm]','Negative is outside cell, positive inside'};
figure(figArgs{:});
hold on
               
for k = 1:nChan
    tmp = branchProfiles.wholeMaskMeanVsDist(k,:) - nanmin(branchProfiles.wholeMaskMeanVsDist(k,:));    
    tmp = tmp ./ nanmax(tmp);
    plot(branchProfiles.wholeMaskDists * movieData.pixelSize_,tmp,'.-','color',chanCols(k,:),'LineWidth',2,'MarkerSize',10);    
end    
legend(chanNames);

%SOOOO LAZY AND TIRED RIGHT NOW - DO IT THIS WAY SO LEGEND MAKES SENSE
for k = 1:nChan
    tmp = branchProfiles.wholeMaskMeanVsDist(k,:) - nanmin(branchProfiles.wholeMaskMeanVsDist(k,:));            
    tmp2 = branchProfiles.wholeMaskSTDVsDist(k,:) / nanmax(tmp);
    tmp = tmp ./ nanmax(tmp);        
    plotTransparent(branchProfiles.wholeMaskDists * movieData.pixelSize_,tmp,tmp2,chanCols(k,:),.3,0);        
end        

xlabel(dLabel)
ylabel('Normalized Average Fluoresence, a.u.')
title({'Distance - Intensity Profile, Whole-Cell','Bands show +/- STD'})
plot([0 0],ylim,'--k')
figName = 'distance vs masked intensity overlay';
print(pOpt{:},[outDir filesep figName '.eps']);
print(pOptTif{:},[outDir filesep figName '.tif']);
hgsave(gcf,[outDir filesep figName '.fig']);

%----2D Intensity histogram vs. Distance figure -----%

figure(figArgs{:})
for k = 1:nChan
    subplot(1,nChan,k)
    hold on
    imagesc(branchProfiles.wholeMaskDists * movieData.pixelSize_,branchProfiles.wholeMaskIntHistBins(k,:),squeeze(branchProfiles.wholeMaskNormHistVsDist(k,:,:))')
    ylim(branchProfiles.wholeMaskIntHistBins(k,[1 end]))%For some reason the limits sometimes exceed these ranges, so force them.
    xlim(branchProfiles.wholeMaskDists([1 end]) * movieData.pixelSize_)
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

avgBranchInt = cell(nChan,1);
for k = 1:nChan    
    figure(figArgs{:});
    avgBranchInt{k} = cellfun(@(x)(mean(x(:,k))),branchProfiles.branchTipPixelInt);
    plot(avgBranchRadius,avgBranchInt{k},'o','MarkerSize',15,'color',chanCols(k,:));    
    ylabel(['Mean ' chanNames{k} ' intensity, a.u.'])
    xlabel('Mean branch radius, nm')   
    
    figName = ['branch radius vs masked intensity in ' chanNames{k}];
    print(pOpt{:},[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);
    
    figure(figArgs{:});
    plot(branchLen,avgBranchInt{k},'o','MarkerSize',15,'color',chanCols(k,:));    
    ylabel(['Mean ' chanNames{k} ' intensity, a.u.'])
    xlabel('Approximate Branch Length, nm')   

    figName = ['branch length vs masked intensity in ' chanNames{k}];
    print(pOpt{:}   ,[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);
               
end


% ----------- Intensity vs. Curvature Figures ------------- %


%TEMP - convert to physical units!!!

%Field names for different curvature measures
curvTypes = {'gaussCurvSampMean',...
             'gaussCurvSampSTD',...
             'meanCurvSampMean',...
             'meanCurvSampSTD',...
             'PC1CurvSampMean',...
             'PC1CurvSampSTD',...
             'PC2CurvSampMean',...
             'PC2CurvSampSTD',...             
             'MaxAbsPCCurvSampMean',...
             'MaxAbsPCCurvSampSTD'};
         
%And intuitive names for them...
curvNames = {'Gaussian Curvature',...
             'STD of Gaussian Curvature',...
             'Mean Curvature',...
             'STD of Mean Curvature',...             
             'k1',...
             'STD of k1',...
             'k2',...
             'STD of k2',...
             'Max Absolute Curvature',...
             'STD of Max Absolute Curvature'};

nCurvTypes = numel(curvTypes);         

nCurvSamp = numel(branchProfiles.gaussCurvSampMean);
nCurvBins = 100;    

for l = 1:nCurvTypes
    
    currCurv = real(branchProfiles.(curvTypes{l}));
    curvBins = linspace(prctile(currCurv,.5),prctile(currCurv,99.5),nCurvBins);
    curvBins = sort(curvBins);%For the strictly negative measures, we need to make sure these are increasing
    
    for k = 1:nChan                                       
        
        
        %currInt = branchProfiles.intForCurvSampMean(:,k);
        currMaxInt = branchProfiles.intForCurvSampMax(:,k);        
            
        figure(figArgs{:})
        subplot(2,1,1)
        hold on        
        N = hist3([currMaxInt currCurv],{branchProfiles.wholeMaskIntHistBins(k,:), curvBins});
        imagesc(branchProfiles.wholeMaskIntHistBins(k,:),curvBins,N');
        colorbar
        caxis([0 prctile(N(:),95)])
        xlim(branchProfiles.wholeMaskIntHistBins(k,[1 end]))
        ylim(curvBins([1 end]))
        xlabel('Intensity, a.u.')
        ylabel(curvNames{l})
        title({['2D Histogram of ' curvNames{l} ' vs. Max Local intensity in Channel ' chanNames{k}],...                
                'Color indicates pixel count'})        
        subplot(2,1,2)
        hold on
        for m = 1:size(N,1)
            N(m,:) = N(m,:) ./ sum(N(m,:));
        end
        imagesc(branchProfiles.wholeMaskIntHistBins(k,:),curvBins,N');
        colorbar
        caxis([0 prctile(N(:),95)])
        xlabel('Intensity, a.u.')
        ylabel(curvNames{l})
        xlim(branchProfiles.wholeMaskIntHistBins(k,[1 end]))
        ylim(curvBins([1 end]))
        title({['2D PDF of ' curvNames{l} ' vs. Max Local Intensity in Channel ' chanNames{k}],...            
            'Color indicates probability'})   
        figName = ['max local intensity in ' chanNames{k} ' vs ' curvNames{l} ' histogram and PDF '];
        print(pOpt{:},[outDir filesep figName '.eps']);
        print(pOptTif{:},[outDir filesep figName '.tif']);
        hgsave(gcf,[outDir filesep figName '.fig']);

        nBinsCur = numel(branchProfiles.wholeMaskIntHistBins(k,:));
        avgCurvVsInt = nan(1,nBinsCur-1);
        stdCurvVsInt = nan(1,nBinsCur-1);
        semCurvVsInt = nan(1,nBinsCur-1);
        nCurvVsInt = nan(1,nBinsCur-1);
        for i = 1:nBinsCur-1
            currCurrVals = currMaxInt >= branchProfiles.wholeMaskIntHistBins(k,i) & currMaxInt < branchProfiles.wholeMaskIntHistBins(k,i+1);
            avgCurvVsInt(i) = mean(currCurv(currCurrVals));
            stdCurvVsInt(i) = std(currCurv(currCurrVals));
            nCurvVsInt(i) = nnz(currCurrVals);
            semCurvVsInt(i) = stdCurvVsInt(i) / sqrt(nCurvVsInt(i));            
        end

        figure(figArgs{:});
        subplot(2,1,1)
        hold on
        title(['Average ' curvNames{l} ' vs. Max Local Intensity in Channel ' chanNames{k}])                    
        plot(branchProfiles.wholeMaskIntHistBins(k,1:end-1),avgCurvVsInt)
        plot(branchProfiles.wholeMaskIntHistBins(k,1:end-1),avgCurvVsInt + 1.96*semCurvVsInt,'--')
        plot(branchProfiles.wholeMaskIntHistBins(k,1:end-1),avgCurvVsInt - 1.96*semCurvVsInt,'--')
        xlabel('Intensity, a.u.')
        ylabel(['Average ' curvNames{l}])
        legend('Average','+/- 1.96SEM')

        subplot(2,1,2)
        hold on
        title(['STD of ' curvNames{l} ' vs. Max Local Intensity in Channel ' chanNames{k}])                         
        plot(branchProfiles.wholeMaskIntHistBins(k,1:end-1),stdCurvVsInt)
        xlabel('Intensity, a.u.')
        ylabel(['STD of ' curvNames{l}])                                               

        figName = ['average and STD of ' curvNames{l} ' vs max local intensity in ' chanNames{k}];
        print(pOpt{:},[outDir filesep figName '.eps']);
        print(pOptTif{:},[outDir filesep figName '.tif']);
        hgsave(gcf,[outDir filesep figName '.fig']);

    end
    
end

if p.BatchMode
    close all
end

%% ---------------- Finalize ------------------ %%

%Save analysis to disk, store directories and parameters in process object
outFile = [outDir filesep fName '.mat'];
save(outFile,'branchProfiles');
movieData.processes_{iProc}.setOutFilePath(iProcChan,outFile);

movieData.processes_{iProc}.setDateTime;
movieData.save;

disp('Finished analyzing masked intensities!')









