function [frmMeanI]=spcklMovWrtIm(params,fluorData,not2keep)
% write the images as 16-bit where dynamic range is 100-500 gray levels

% imaging parameters
rAiry=(1.22*params.lambdaEmit)/(2*params.na); % Airy disk radius (nm)
sigma=rAiry/3; % sigma for gaussian approx of PSF (nm)
frmMeanI=zeros(params.nFrames,1);

% get correct string for number of files
s=length(num2str(params.nFrames));
strg=sprintf('%%.%dd',s);

% create image - each pixel is now pixNM x pixNM
for fm=1:params.nFrames
    smoothSampledIm=zeros(not2keep.nSamplesL,not2keep.nSamplesW);
    startStep=(fm-1)*not2keep.nTmPtsPerFm+1; % first time pt per frame
    col=find(not2keep.pts2Store==startStep); % column in fluorData.state in which first time point is stored
    
    % add up gaussian-smoothed images over integration time for current
    % frame
    for i=col:col+not2keep.nTmPtsStoredPerFm-1
        % convert location from nm to which sub-box the fluorophore is in
        visBoundLoc=xy2index(fluorData.px(fluorData.state(:,i)~=0,i),...
            fluorData.py(fluorData.state(:,i)~=0,i),not2keep.nSamplesL,not2keep.nSamplesW,params.sampleScale);

        % get how many fluor. fell into each sampled area
        nPts=zeros(not2keep.nSampleSquares,1);
        [a b]=countEntries(visBoundLoc);
        nPts(a(~isnan(a)))=b(~isnan(a));
        sampledIm=zeros(not2keep.nSamplesL,not2keep.nSamplesW);
        sampledIm(:)=nPts; % let intensity of each fluorophore be 1
        smoothSampledIm=smoothSampledIm+Gauss2D(sampledIm,sigma/params.sampleScale); % do Gaussian smoothing

    end
    % resample the image to pixels (digitization step)
    ccdIm=imResample(smoothSampledIm, [params.sampleScale params.sampleScale], [params.pixNM params.pixNM]);
    % crop off the borders
    ccdIm=ccdIm(params.border.top+1:params.border.top+params.imL,params.border.left+1:params.border.left+params.imW);
    % save as a .mat file
    indxStr=sprintf(strg,fm);
    save([params.outputDirTifs filesep 'ccd_image' indxStr],'ccdIm');
end

[listOfImages] = searchFiles('.mat',[],params.outputDirTifs);
nImTot=length(listOfImages);

% get movie min/max values over all the images
movMin=10^9; movMax=0;
for fm=1:params.nFrames
    ccdImFile=[char(listOfImages(fm,2)) filesep char(listOfImages(fm,1))];
    ccdIm=load(ccdImFile);
    
    movMin=min(movMin,min(ccdIm.ccdIm(:)));
    movMax=max(movMax,nanmax(ccdIm.ccdIm(:)));
end

% normalize images to 100-500 and save as 16-bit tiffs
for fm=1:params.nFrames
    ccdImFile=[char(listOfImages(fm,2)) filesep char(listOfImages(fm,1))];
    ccdIm=load(ccdImFile);
    ccdIm=ccdIm.ccdIm;
    
    % put in range 100-500 and convert to 16-bit image
    ccdIm=400*(ccdIm-movMin)/(movMax-movMin)+100;
    ccdIm=uint16(ccdIm);

    % write the file
    indxStr=sprintf(strg,fm);
    % to write labeling ratio into file name
    extra=['_',];
    imwrite(ccdIm,[params.outputDirTifs filesep 'image' indxStr '.tif']);

    % get mean intensity of the frame
    frmMeanI(fm)=mean(ccdIm(:));
end

delete([params.outputDirTifs filesep '*.mat'])