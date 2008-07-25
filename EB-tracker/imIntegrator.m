function imIntegrator(runInfo,nFms2avg)

if nargin<2
    error('imIntegrator: need 2 input arguments')
end

% check runInfo format and directory names between file systems
if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('eb1SpotDetector: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir]=formatPath(runInfo.anDir);
    [runInfo.imDir]=formatPath(runInfo.imDir);
end

% get roi edge pixels and make region outside mask NaN
if ~isfield(runInfo,'roiMask')
    roiMask=1;
    edgePix=[];
else
    polyEdge=bwmorph(runInfo.roiMask,'remove');
    edgePix=find(polyEdge);
    roiMask=double(runInfo.roiMask);
    %roiMask=swapMaskValues(roiMask,0,NaN);
end

% make spot directory if it doesn't exist from batch
intImDir=[runInfo.anDir filesep 'spot' filesep 'intIm'];
if ~isdir(intImDir)
    mkdir(intImDir);
end

% count number of images in image directory
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
nImTot=size(listOfImages,1);

s1=length(num2str(nImTot));
strg1=sprintf('%%.%dd',s1);


for i=1:nImTot
    
    % set start/end frame numbers to integrate over
    if i<ceil(nFms2avg/2) 
        %first few: take first nFms2avg frames
        sF=1;
        eF=nFms2avg;
    elseif i>nImTot-floor(nFms2avg/2)
        %last few: take last nFms2avg frames
        sF=nImTot-nFms2avg+1;
        eF=nImTot;
    else
        %middle frames: take frames before and after i
        sF=i-(nFms2avg-1)/2; 
        eF=i+(nFms2avg-1)/2; 
    end

    % get mean image
    img=0;
    for j=sF:eF
        fileNameIm=[char(listOfImages(j,2)) filesep char(listOfImages(j,1))];
        img=img+double(imread(fileNameIm));
    end
    meanImg=img./nFms2avg;
    
    % normalize, mask out, and save as 8-bit image
    meanImg=uint8(roiMask.*round(255.*(meanImg-min(meanImg(:)))./(max(meanImg(:))-min(meanImg(:)))));
    indxStr1=sprintf(strg1,i);
    imwrite(meanImg,[intImDir filesep 'meanImg' indxStr1 '.tif']);
    

end

