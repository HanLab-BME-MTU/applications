function success=analyzeLMdata(dataDirectory,wavelength,varargin)
%ANALYZELMDATA analyzes localization microscopy data
%   Input:
%         required arguments
%           dataDirectory -> directory of experimental data (required)
%              wavelength -> wavelength(s), must be in proper ordering
%
%         optional arguments
%           NA  -> numerical aperture, default 1.49
%           MAG -> magnification, default 100x
%  
%   Output:
%         success -> 1/0: analysis finished successfully/unsuccessfully
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

ip.parse(dataDirectory,varargin{:});

% remove trailing '/' from dataDirectory string
if dataDirectory(end) == filesep
    dataDirectory=dataDirectory(1:end-1);
end

% get list of experiment files
fileList=dir([dataDirectory filesep '*.dv']);
nFiles=numel(fileList);

nChannels=numel(ip.Results.wavelength);

for iFile=1:nFiles
    % load movie data and set output directory
    movieFile=[dataDirectory filesep fileList(iFile).name];
    MD=MovieData.load(movieFile,false);
    % set numerical aperture and magnification manually
    MD.numAperture_=ip.Results.NA;
    MD.magnification_=ip.Results.MAG;
    
    MD.sanityCheck();
    
    NA=MD.numAperture_;
    MAG=MD.magnification_;
    pixelSize=MD.pixelsize_;
    
    for iChannel=1:nChannels
        lambda=ip.Results.wavelength(iChannel);
        psfSigma=getgaussianPSFsigma(NA,MAG,pixelSize*MAG*1e-6,lambda*1e-9);
        psfSigmaInt=min(1,ceil(psfSigma));
        
    end
    
end

end

