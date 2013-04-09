function success=analyzeLMdataTracking(dataDirectory,wavelength,varargin)
%ANALYZELMDATA analyzes localization microscopy data from DeltaVision files
%   using pointSourceDetection followed by trackCloseGapsKalmanSparse and
%   saves the results
%
%   Input:
%         required arguments
%           dataDirectory  -> directory of experimental data
%              wavelength  -> wavelength(s) [nm], must be in proper ordering
%
%         optional arguments
%            NA  ->  numerical aperture, default: 1.49
%           MAG  ->  magnification, default: 100x
%           rep  ->  number of frames in activation cycle, default: 10
%    createMask  ->  should a mask be created? default: false
%         doMMF  ->  use mixture model fitting default: false
%        GapLen  ->  maximum gap size to be closed, default: 4
%        Radius  ->  maxium search radius for linking defalt: 2
%  
%   Output:
%         success  ->  1/0: analysis was successful/unsuccessful
%
% Ulrich Schmidt, March 13, 2012
% Jeffrey Werbin, March 2013

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
ip.addOptional('GapLen',4,@isnumeric);
ip.addOptional('Radius',2,@isnumeric);

ip.parse(dataDirectory,wavelength,varargin{:});

NA=ip.Results.NA;
MAG=ip.Results.MAG;
rep=ip.Results.rep;
createMask=ip.Results.createMask;
doMMF=ip.Results.doMMF;
GapLen=ip.Results.GapLen;
Radius=ip.Results.Radius;

%% Creates parameter structures for tracking
%% general gap closing parameters
gapCloseParam.timeWindow = GapLen; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatStationaryLink';

%parameters
parameters.searchRadius = Radius;

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatStationaryCloseGaps';

%parameters
parameters.searchRadius = Radius; 
parameters.gapPenalty = 2;

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions = [];

%% additional input
saveResults = 0;

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% Make MovieData files for each .dv movie in the directory and run detection and tracking one at a time

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

for iFile=1:nFiles
    % load movie data and set output directory
    movieFile=[dataDirectory filesep fileList(iFile).name];
    outputDir=[dataDirectory filesep fileList(iFile).name(1:end-3) filesep];
    MD=MovieData.load(movieFile,false,'outputDirectory',outputDir);
    % set numerical aperture and magnification manually
    MD.numAperture_=NA;
    MD.magnification_=MAG;
    
    MD.sanityCheck();
    % a DV file stores the pixelsize in nm
    %pixelSize=MD.pixelSize_;
    pixelSize = 62.8100;
    
    if createMask
        % create ROI and store in MD
        img=double(MD.channels_(1).loadImage(1));
        mask=double(roipoly(img/max(img(:))));
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
        mask=ones(MD.imSize_);
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

    iFrame=0;
    %shift=10;
    
    
    fprintf(1,'movie being analyzed: %s\n',fileList(iFile).name);
    while iFrame < nFrames
 
        k=mod(floor(iFrame/rep),nChannels);
        sigma=psfSigmaTheo(k+1);
        
        img=double(MD.channels_(1).loadImage(iFrame+1));
        features{iFrame+1}=pointSourceDetection(img,sigma,'alpha',1e-9,...
            'mask',mask,'FitMixtures',doMMF);
        
       if ~isempty(features{iFrame+1})
        movieInfo(iFrame+1).xCoord = [features{iFrame+1}.x',features{iFrame+1}.x_pstd'];
        movieInfo(iFrame+1).yCoord = [features{iFrame+1}.y',features{iFrame+1}.y_pstd'];
        movieInfo(iFrame+1).amp = [features{iFrame+1}.A',features{iFrame+1}.A_pstd'];
       else
        movieInfo(iFrame+1).xCoord = [[],[]];
        movieInfo(iFrame+1).yCoord = [[],[]];
        movieInfo(iFrame+1).amp = [[],[]];
       end
           
        iFrame=iFrame+1;
        
        progressText(iFrame/nFrames,'All work and no play makes Jack a dull boy');
    end
    
    fprintf(1,'\n');
    
    %Apply tracking and Gap closing to localized data
    [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
    
    
    % save results
    outputDir=MD.outputDirectory_;
    outputDir=[outputDir filesep 'results'];
    if ~exist(outputDir,'dir');
        mkdir(outputDir);
    end
    outputFile=[MD.movieDataFileName_(1:end-7) '_tracking.mat'];
    save([outputDir filesep outputFile],'features','tracksFinal','GapLen','Radius');
        
end