function [ output_args ] = gcaMakeVeilStemMaskFolder(movieData,varargin)
%gcaMakeVeilStemMaskFolder

%% Check input
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

ip = inputParser;

ip.CaseSensitive = false;


defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
    'III_veilStem_reconstruction' filesep 'Channel_1'];

defaultOutDir = [defaultInDir filesep 'finalMask']; 

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('InputDirectory', defaultInDir,@(x) ischar(x));

ip.addParameter('addExternalSegProcess',true);
defaultOutDirMD = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'masks_VeilStemFinal'];
ip.addParameter('OutputDirectoryMD',defaultOutDirMD,@(x) ischar(x)); 

ip.addParameter('Overlays',true);
ip.parse(varargin{:});
%% Initialize
if ~isdir([ip.Results.OutputDirectory])
    mkdir([ip.Results.OutputDirectory]); 
end 

fileToLoad = [ip.Results.InputDirectory filesep 'veilStem.mat'];
load(fileToLoad);
nImTot = numel(veilStem);
fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];
%%
for iFrame = 1:nImTot
    mask  = veilStem(iFrame).finalMask;
    
    % should have already been done but just to make sure
    % get largest CC
    mask = double(getLargestCC(mask));
    
    % set the border to zero - this just helps make sure that the windows do not
    % error
    dims = size(mask);
    mask(1:dims(1),1) =0;
    mask(1:dims(1),dims(2))=0;
    mask(1,1:dims(2))= 0;
    mask(dims(1),1:dims(2)) =0;
    
    % should have been filled but again just in case
    
    mask = imfill(mask,'holes');
    imwrite(mask,[ip.Results.OutputDirectory filesep 'veilStemMask' num2str(iFrame,fmt) '.tif']);
    if ip.Results.Overlays
        
        img =  movieData.channels_(1).loadImage(iFrame);
        img = double(img);
        [ny,nx] = size(img);
        setFigure(nx,ny,'off');
        imshow(-img,[]);
        hold on
        roiYX = bwboundaries(mask);
        text(5,10,'VeilStem Outline', 'Color','g');
        cellfun(@(x) plot(x(:,2),x(:,1),'color','g'),roiYX);
        if ~isdir([ip.Results.OutputDirectory filesep 'Overlays']);
            mkdir([ip.Results.OutputDirectory filesep 'Overlays']) ;
        end
        saveas(gcf,[ip.Results.OutputDirectory filesep 'Overlays' filesep 'VeilStemOverlay' num2str(iFrame,fmt) '.png']);
        close gcf
    end
    clear mask
end
% add to process via external segmentation
if ip.Results.addExternalSegProcess
    
    extProc = ExternalSegmentationProcess(movieData);
    % % % Save the process in the movie object
    movieData.addProcess(extProc);
    params = movieData.processes_{end}.funParams_;
    params.OutputDirectory = ip.Results.OutputDirectoryMD;
    params.InputData{1} = ip.Results.OutputDirectory; 
    parseProcessParams(movieData.processes_{end},params);
    movieData.processes_{end}.run;
    movieData.save;   
end
end  %function
