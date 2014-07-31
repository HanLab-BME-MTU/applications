function iceConn = movieViewerImaris(MD,varargin)
%MOVIEVIEWERIMARIS views the 3D dataset described by the input movieData in imaris
%
% movieViewerImaris
% movieViewerImaris(MD)
% iceConn = movieViewerImaris(MD,'OptionName1',optionValue1,...)
%
%   Input: MD - MovieData object describing a 3D dataset
%
%   Output: iceConn - IceImarisConnector object for interacting with imaris instance
%
%Hunter Elliott
%7/2014


%% ------- Input --------- %%

if nargin < 1 || isempty(MD)
    [filePath,fileName] = optionalFileInput('','*.mat','Select the MovieData file to view:');
    MD = MovieData.load([filePath fileName]);        
elseif ~isa(MD,'MovieData') || ~MD.is3D
    error('First input must be a MovieData object describing a 3D dataset!')
end


ip = inputParser;
ip.addParamValue('IceConnector',[],@(x)(isa(x,'IceImarisConnector')));%Optinally input an IceConnector objects
ip.parse(varargin{:});
p = ip.Results;

%% --------- Init ------------ %%

imSize = [MD.imSize_ MD.zSize_];
pixSize = [MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_] / 1e3;%imaris defaults to microns
nFrames = MD.nFrames_;
nChan = numel(MD.channels_);
if ~isempty(MD.timeInterval_)    
    tInt = MD.timeInterval_;
else
    tInt = 1;
end

if ~isempty(MD.camBitdepth_)
    bitDepth = MD.camBitdepth_;
else
    %Just assume 16
    bitDepth = 16;
end
imClass = ['uint' num2str(bitDepth)];

if isempty(p.IceConnector)
    iceConn = IceImarisConnector;
end

%Start imaris
iceConn.startImaris;

%Create the dataset
dataSet = iceConn.createDataset(imClass,imSize(1),imSize(2),imSize(3),nChan,nFrames,...
                                       pixSize(1),pixSize(2),pixSize(3),tInt);

%Set image properties. 
chanCols = jet(nChan);%Default colors
for iChan = 1:nChan    
    if ~isempty(MD.channels_(iChan).emissionWavelength_)
        chanCols(iChan,:)= MD.channels_(iChan).getColor;
        dataSet.SetChannelName(iChan-1,num2str(MD.channels_(iChan).emissionWavelength_));        
    end
    dataSet.SetChannelColorRGBA(iChan-1,iceConn.mapRgbaVectorToScalar([chanCols(iChan,:) 0]));
end

dispMin = Inf(1,nChan);
dispMax = zeros(1,nChan);
%Set the name so we can check which dataset is loaded
dataSet.SetParameter('Image','Name',[MD.outputDirectory_ filesep MD.movieDataFileName_]);

%% ------ Image Loading / Passing ----- %%

wtBar = waitbar(0,'Please wait, loading all image data...');

for iFrame = 1:nFrames    
    for iChan = 1:nChan                            
            %NOTE - for very large images this has lots of overhead and may
            %cause problems passing to imaris. May need to switch to
            %single-plane loading/passing
            currIm = MD.channels_(iChan).loadStack(iFrame); 
            if ~strcmpi(class(currIm),imClass)
                warning('movViewImaris:dataTypeMismatch',['Image data class mismatch: expected ' imClass ' was loaded as ' class(currIm) ' : casting!']);
            end
            %Pass the image to imaris
            iceConn.setDataVolume(cast(currIm,imClass),iChan-1,iFrame-1);
            
            dispMin(iChan) = min(dispMin(iChan),min(currIm(:)));
            dispMax(iChan) = max(dispMax(iChan),max(currIm(:)));
    end  
    
    waitbar(iFrame/nFrames,wtBar)
    
end

if ishandle(wtBar)
    close(wtBar)
end


%% -------- Display Configuration ------ %%

for iChan = 1:nChan
    
   dataSet.SetChannelRange(iChan-1,dispMin(iChan),dispMax(iChan)); 
    
end


%% ------- Process Output Display ------ %%

nProc = numel(MD.processes_);

for iProc = 1:nProc
    
    if ismethod(MD.processes_{iProc},'drawImaris')        
       MD.processes_{iProc}.drawImaris(iceConn); 
    end
    
end


