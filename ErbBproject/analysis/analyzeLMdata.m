function success=analyzeLMdata(dataDirectory,wavelength,varargin)
%ANALYZELMDATA analyzes localization microscopy data from DeltaVision files
%   Input:
%         required arguments
%           dataDirectory  -> directory of experimental data
%              wavelength  -> wavelength(s) [nm], must be in proper ordering
%
%         optional arguments
%            NA  ->  numerical aperture, default 1.49
%           MAG  ->  magnification, default 100x
%           rep  ->  number of frames in activation cycle, default 10
%    createMask  ->  should a mask be created? default: false
%  
%   Output:
%         success  ->  1/0: analysis was successful/unsuccessful
%
% Ulrich Schmidt, March 13, 2012
%

success=1;

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('dataDirectory',@ischar);
ip.addRequired('wavelength',@isnumeric);

ip.addOptional('NA',1.49,@isscalar);
ip.addOptional('MAG',100,@isscalar);
ip.addOptional('rep',10,@isnumeric);
ip.addOptional('createMask',false,@islogical);

ip.addOptional('doMMF',false,@islogical);

ip.parse(dataDirectory,wavelength,varargin{:});

NA=ip.Results.NA;
MAG=ip.Results.MAG;
rep=ip.Results.rep;
createMask=ip.Results.createMask;

doMMF=ip.Results.doMMF;

% remove trailing '/' from dataDirectory string
if dataDirectory(end) == filesep
    dataDirectory=dataDirectory(1:end-1);
end

% get list of experiment files
fileList=dir([dataDirectory filesep '*.dv']);
nFiles=numel(fileList);

wavelength=ip.Results.wavelength;
nChannels=numel(wavelength);
psfSigmaTheo=ones(size(wavelength));
psfSigmaInt=ones(size(wavelength));


for iFile=5:6
    % load movie data and set output directory
    movieFile=[dataDirectory filesep fileList(iFile).name];
    outputDir=[dataDirectory filesep];
    MD=MovieData.load(movieFile,false,'outputDirectory',dataDirectory);
    % set numerical aperture and magnification manually
    MD.numAperture_=NA;
    MD.magnification_=MAG;
    
    %MD.sanityCheck();
    % a DV file stores the pixelsize in nm
    pixelSize=MD.pixelSize_;
    
    if createMask
        % create ROI and store in MD
        img=double(MD.channels_(1).loadImage(1));
        mask=roipoly(img/max(img(:)));
        close all;
        maskFile=MD.getFilename();
        maskFile(end-2:end)='tif';
        maskPath=MD.getPath();
        fullMaskPath=[maskPath filesep maskFile];
        imwrite(mask,[maskPath filesep maskFile],'tif');
        MD.addROI(fullMaskPath,maskPath);
        
        maskFile=MD.getFilename();
        maskFile(end:end+1)='sk';
        MD.rois_(1).setFilename(maskFile);
        MD.rois_(1).setPath(maskPath);
        MD.rois_(1).save();
    else
        mask=true(MD.imSize_);
    end
    
    MD.save();
    
    for iChannel=1:nChannels
        % calculate theoretical PSF sigma
        lambda=wavelength(iChannel);
        psfSigmaTheo(iChannel)=getGaussianPSFsigma(NA,MAG,pixelSize*1e-9*MAG,lambda*1e-9);
        psfSigmaInt(iChannel)=max(1,ceil(psfSigmaTheo(iChannel)));
    end
    
    % total number of frames in DV file
    nFrames=MD.nFrames_;
    
    features=cell(nFrames,1);

    iFrame=100;
    shift=10;
    fprintf(1,'movie being analyzed: %s\n',fileList(iFile).name);
    while iFrame < nFrames
 
        k=mod(floor(iFrame/rep),nChannels);
        sigma=psfSigmaInt(k+1);
        
        img=double(MD.channels_(1).loadImage(iFrame+1));
        locMax=findLocalMaxima(img,sigma,'alpha',1e-6,'mask',mask);
        
        if ~numel(locMax.x)
            iFrame=iFrame+1;
            continue;
        end
        
        sigma=psfSigmaTheo(k+1);
        f=mixtureModelFitting(img,locMax,sigma,'doMMF',doMMF);
        
        features{iFrame+1}=f;
        
        iFrame=iFrame+1;
        
        
        %if( iFrame == shift )
        %    iFrame=iFrame+10;
        %    shift=shift+20;
        %end
        
        progressText(iFrame/nFrames,'All work and no play makes Jack a dull boy');
    end
    
    fprintf(1,'\n');
    
    % save results
    outputDir=MD.outputDirectory_;
    outputDir=[outputDir filesep 'results'];
    if ~exist(outputDir,'dir');
        mkdir(outputDir);
    end
    outputFile=[MD.movieDataFileName_(1:end-8) '.mat'];
    save([outputDir filesep outputFile],'features');
        
end