function [ filoBranch ] = GCASubRegionalExtractFilopodiaMovie( movieData,varargin)
% GCASubRegionalGetWindowsMovie: finds windows in a subroi that last
% throughout the time interval chosen 
%
%%
%
% OUTPUT: 
% idxWindFinal 
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
%ip.addParameter('TSOverlays',true);
defaultSubRoiDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep  'GC' filesep 'masks'];
ip.addParameter('subRoiMaskDirIn',defaultSubRoiDir); % where the subRoimasks are stored.


ip.addParameter('OutputDirectory', []);

defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits']; 
ip.addParameter('InputDirectory', defaultInDir); 
 
ip.addParameter('EndFrame',movieData.nFrames_-1);
ip.addParameter('StartFrame',1);

ip.addParameter('plotOverlays',true,@(x) islogical(x)); 

ip.parse(varargin{:});

% 
if isempty(ip.Results.OutputDirectory)
    x = upDirectory(ip.Results.subRoiMaskDirIn,1);
    %defaultOutDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep 'GC' filesep 'filopodia'];
    outDir = [x filesep 'filopodia']; 
else
    outDir = ip.Results.OutputDirectory;
end
%% Initiate
% load the filoInfo with the fitted data only problem is this is quite a
% large file 
load([ip.Results.InputDirectory filesep 'Channel_1' filesep 'filoBranch.mat']); % NOTE fix channel input here.  
    % rename 
    filoBranchWhole = filoBranch; 
    clear filoBranch
% make folders 
if ~isdir(outDir)
    mkdir(outDir)
end 
%% % Get SubRegional Windows 

if ~isempty(ip.Results.subRoiMaskDirIn)
    
 
for iFrame = ip.Results.StartFrame:ip.Results.EndFrame
    filoInfo = filoBranchWhole(iFrame).filoInfo; 
    roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
    
    % load the subRoiMask
    [filoInfoSub] = GCAextractFilopodiaFromSubRegions(filoInfo,roiMask);
    
    filoBranch(iFrame).filoInfo = filoInfoSub; 
    filoBranch(iFrame).timeStampSub = clock; 
    filoBranch(iFrame).timeStampWhole = filoBranchWhole(iFrame).timeStamp; 
    clear filoInfoSub 
    save([outDir filesep 'filoBranch.mat'  ],'filoBranch'); 
end
close gcf

if ip.Results.plotOverlays
    visualDir = [outDir filesep 'Visualization' ]; 
    if ~isdir(visualDir) 
        mkdir(visualDir)
    end 
    
    for  iFrame = ip.Results.StartFrame:ip.Results.EndFrame
       
        img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}]));
        [ny,nx] = size(img);
        setFigure(nx,ny,'off')
        imshow(-img,[]);
        filoInfoSub = filoBranch(iFrame).filoInfo; 
        hold on
       
        roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
        roiYX = bwboundaries(roiMask); 
        cellfun(@(x) plot(x(:,2),x(:,1),'color','r'),roiYX); 
    
        GCAVisualsFilopodiaMeasurementOverlays(filoInfoSub,[ny,nx]);
        saveas(gcf,[visualDir filesep num2str(iFrame,'%03d') '.png']); 
        close gcf 
    end
end

% find start and end points


else 
    filoInfoSub = []; 
end  
  %% end subRegionalSpecific   
  
  
  


