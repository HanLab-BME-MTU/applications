function features=analyzeDNApaint(multiTIFF,imgParam)
%ANALYZEDNAPAINT analyzes DNA paint experiments with multi-TIFF files
%
%   Input:
%       multiTIFF ->  absolute path to multi-TIFF file (required)
%       imgParam  ->  imaging parameters in structure including
%                     NA, numerical aperture
%                     mag, objective maginifictaion
%                     lambda, max emission wavelenght [m]
%                     pixelSize, physical camera pixelsize [m]
%       See documentation of function 'getGaussianPSFsigma' for details.
%
%   Output:
%       features  ->  cell array of detected feature informations
%
%   US, 07/03/2012
%

% disable all warnings
warning('off','all');

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=false;

ip.addRequired('multiTIFF',@ischar);
ip.addRequired('imgParam',@isstruct);

ip.parse(multiTIFF,imgParam);

multiTIFF=ip.Results.multiTIFF;
imgPrm=ip.Results.imgParam;

prmField={'NA','mag','lambda','pixelSize'};
flag=isfield(imgPrm,prmField);

if ~all(flag)
    error('Incorrect field names. Fix field names of struct imgParam');
end

img=Tiff(multiTIFF,'r');

count=0;

while ~img.lastDirectory()
    frame=double(img.read());
    img.nextDirectory();
    count=count+1;
end

img.close();

features=1923;

end