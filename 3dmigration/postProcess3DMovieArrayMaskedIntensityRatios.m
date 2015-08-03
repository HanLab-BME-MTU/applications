function postProcess3DMovieArrayMaskedIntensityRatios(MA,p)

%This function is horribly written and shouldn't be used.

nChan = numel(MA(1).channels_);

p.DilationRadius = 10;
p.ChannelIndex= 1:nChan;
p.SampleRadius = 2e3;
p.CorrChan = 1:nChan; %Only the actin channel need


iProcChan = 1;
iFrame = 1;%All fixed so 1 frame...
nMov = numel(MA);



for iMov = 1:nMov
    
    currOutDir = [MA(iMov).outputDirectory_ filesep 'intensity ratio analysis'];
    mkClrDir(currOutDir);
    
    %Get image locations and info
    imDir = MA(iMov).getChannelPaths(p.ChannelIndex);
    imNames = MA(iMov).getImageFileNames(p.ChannelIndex);

    
    %nZ = MA(iMov).nSlices_;
    
    iSegProc = MA(iMov).getProcessIndex('SegmentationProcess3D',1,1);
    %Get mask file names and directory
    maskDir = MA(iMov).processes_{iSegProc}.outFilePaths_{1,iProcChan};
    maskNames = MA(iMov).processes_{iSegProc}.getOutMaskFileNames(iProcChan);
    mask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
    mask = make3DImageVoxelsSymmetric(mask,MA(iMov).pixelSize_,MA(iMov).zSpacing_);
    nZ = size(mask,3);
    %backMask = ~imdilate(mask,binarySphere(p.DilationRadius));                
    cortMask = bwdist(~mask) <= (p.SampleRadius / MA(iMov).pixelSize_) & mask;
    
    for iChan = 1:nChan
        
        tmp = stackRead([imDir{iChan} filesep imNames{iChan}{iFrame}]);               
        tmp = make3DImageVoxelsSymmetric(tmp,MA(iMov).pixelSize_,MA(iMov).zSpacing_);
        
        %backMask = estimateBackgroundArea(tmp);
        %
    
        
        clear backMask backMeanVsZ backSTDVsZ imCorr 
        for iZ = 1:nZ
           
            backMask(:,:,iZ) = estimateBackgroundArea(tmp(:,:,iZ));
            currPlane = tmp(:,:,iZ);
            backMeanVsZ(iZ) = mean(currPlane(backMask(:,:,iZ)));
            backSTDVsZ(iZ) = std(single(currPlane(backMask(:,:,iZ))));
                 
            if any(iChan ==p.CorrChan)
                imCorr(:,:,iZ) = single(currPlane) - backMeanVsZ(iZ);        
            else
                imCorr(:,:,iZ) = single(currPlane);
            end
            
        end               
        
        cf = figure;
        %plot(backMeanVsZ,1:nZ,'.-')
        errorbar(backMeanVsZ,backSTDVsZ)
        xlabel('Z Plane')
        ylabel('Mean Background Intensity, a.u.')
        mfFigureExport(cf,[currOutDir filesep 'Background vs z channel ' num2str(iChan)]);
        
        allSamps{iMov}(:,iChan) = imCorr(mask(:));
        allCortSamps{iMov}(:,iChan) = imCorr(cortMask(:));
        
        jkl=1;
             
    end
    
    
    
end

allInt= vertcat(allCortSamps{iMov});
for j = 1:nChan
    allInt(:,j) = allInt(:,j) - min(allInt(:,j));
end

cf = figure;

nBins = 100;
%xl = prctile(allInt(:,3),[0.01 99.9]);
%yl = prctile(allInt(:,1),[0.01 99.9]);
xl = [0 450];
yl = [0 1e3];

xBins = linspace(xl(1),xl(2),nBins);
yBins = linspace(yl(1),yl(2),nBins);

[N,C] = hist3(allInt(:,[1 3]),{yBins xBins});
nLog = log10(N);
nLog(~isfinite(nLog)) = 0;
imagesc(C{2},C{1},nLog');


% xlim([0 450])
% ylim([0 1e3])


cl = colorbar;
%set(cl,'YLabel','Log(Voxel Count)')
saturateImageColormap(gca,5)
ylabel('MII-GFP Intensity (a.u.)')
xlabel('phalloidin Intensity (a.u.)')
axis xy

%figName = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Figure Panels for Revision\Bleb actin vs MII color is log10 pixel count';
figName = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Figure Panels for Revision\Bleb actin vs MII color is log10 pixel count equal axes';
%figName = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Figure Panels for Revision\WT actin vs MII color is log10 pixel count equal axes';
mfFigureExport(cf,figName);

save([pwd filesep 'WT MII and actin ratio wkspc.mat'])

