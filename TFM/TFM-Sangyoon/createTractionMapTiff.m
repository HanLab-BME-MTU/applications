function []=createTractionMapTiff(movieData,varargin)
% function []=createTractionMap(MD) creates a heatmap of traction magnitude
% in the coordinates of undeformed position (so that they can be used for
% other applications that does not need 
% without exaggerating it with a constant. It stores the map TIF image to
% /TFMPackage/TractionMapTiff.
% Sangyoon Han 2019 January

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('OutputTFMProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(OutputTFMProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
outputTFMProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(outputTFMProc,paramsIn);
% p.lastToFirst = false;

%% Load TFMPackage
nFrames = movieData.nFrames_;
% Get TFM package
TFMPackage = movieData.getPackage(movieData.getPackageIndex('TFMPackage'));

%% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
pathFF=TFMPackage.outputDirectory_;

padZeros=floor(log10(nFrames))+1;
if p.useCellConfig
    tractionImgFolder=[pathFF filesep 'TractionMapTiff_rawChannelCoord'];
    tMap=forceFieldProc.loadChannelOutput('output','tMapUnshifted'); %'tMap');% in Pa per pixel (1pix x 1pix)
elseif p.useRefConfig
    tractionImgFolder=[pathFF filesep 'TractionMapTiff_SDCCoord'];
    tMap=forceFieldProc.loadChannelOutput('output','tMap'); %'tMap');% in Pa per pixel (1pix x 1pix)
end
%% Image writing
if ~exist(tractionImgFolder,'dir')
%     system(['mkdir -p ' tractionImgFolder]);
    mkdir(tractionImgFolder);
end
progressText(0,'creating traction mag tiff image') % Update text

for ii=1:nFrames
    %imwrite(uint16(Mblue),[target_dir,filesep,'Force_Magnitude',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none') 
    curTmap = tMap(:,:,ii);
    imwrite(uint16(curTmap),[tractionImgFolder,filesep,'ForceMag',num2str(ii,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none');
    progressText(ii/nFrames,'creating traction mag tiff image') % Update text
end
    
disp('Tiff TFM map creation done!')
disp(['Saved in ' tractionImgFolder])