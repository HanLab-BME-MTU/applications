function makeWindowWithColorCodedProtrusionState(movieData,...
                                                 iBands,...
                                                 iSectors,...
                                                 iChan,...
                                                 movieName)
                                                  
if nargin < 5 || isempty(movieName)
    movieName = 'windowTestingMovie';
end
    
if nargin < 4 || isempty(iChan)
    iChan = 1;
end
                                                      
%Load the windows
disp('Loading Windows...')
load([movieData.windows.directory filesep movieData.windows.fileName])

%Determine number of bands, windows and images
[nBandsTot,nSectorsTot,nFrames] = size(allWinPoly);

%Load protrusion samples
disp('Loading protrusion samples...');
load(fullfile(movieData.protrusion.directory, movieData.protrusion.samples.fileName));
if ~isfield(protrusionSamples,'stateNames') || ...
        ~isfield(protrusionSamples,'states') || ...
        ~isfield(protrusionSamples,'statePersistence')
    error('Classification of edge velocity has not been performed.');
end

states = protrusionSamples.states;

if exist('allWindowSamples','var') && isfield(allWindowSamples,'nPixels')
    winStati = ~isnan(allWindowSamples.nPixels);
else
    disp('Cannot find window samples, displaying all windows!');
    winStati = ones(nBandsTot,nSectorsTot,nImages); %TEMP TEMP!
end

%Default is to use all windows
if nargin < 3 || isempty(iSectors)
    iSectors = 1:nWindowsTot;
end
%Default is to use all bands
if nargin < 2 || isempty(iBands)
    iBands = 1:nBandsTot;
end

%Determine number of windows/bands requested for tracking
nBands = length(iBands);
nSectors = length(iSectors);

disp('Checking image files...')
imageFileNames = dir([movieData.imageDirectory filesep movieData.channelDirectory{iChan} filesep '*.tif']); %TEMP TEMP

try
    disp('Checking protrusion vectors....')
    load([movieData.protrusion.directory filesep movieData.protrusion.fileName])
    useProt = true;
catch errMsg
    disp(['Warning: ' errMsg]);
    useProt = false;
end

MakeQTMovie('start',[movieData.analysisDirectory filesep movieName '.mov']);
MakeQTMovie('quality',.75);
        
for iFrame = 1:nFrames-1
    clf
    hold on
    
    currImage = imread([movieData.imageDirectory filesep movieData.channelDirectory{iChan} filesep imageFileNames(iImage).name]); %TEMP TEMP!!
    
    currImage = cast(currImage,'double');
    imagesc(currImage), colormap gray, axis image, axis off;
    
    if iImage == 1
        %Find the mean and standard deviation of first image values
        imMean = nanmean(nonzeros(cast(currImage,'double')));
        imStd = nanstd(nonzeros(cast(currImage,'double')));        
    end
    %Use this to set the color axis on all frames
    caxis([ (imMean - 1.7*imStd) (imMean + 2.5*imStd)]); 
    set(gca,'color','k')

    set(gca,'YDir','reverse')
    %Overlay each requested window on the image

    % Color definition for pause, protrusion and retraction states
    colors = {'y','b','r'};
    
    for iSector = 1:nSectors
        for iBand = 1:nBands            
            if ~isempty(allWinPoly(iBands(iBand),iSectors(iSector),iFrame).outerBorder) &&  ~isempty(allWinPoly(iBands(iBand),iSectors(iSector),iFrame).innerBorder)
                xTmp = [ allWinPoly(iBands(iBand),iSectors(iSector),iFrame).outerBorder(1,:) ...
                         allWinPoly(iBands(iBand),iSectors(iSector),iFrame).innerBorder(1,end:-1:1)];
                yTmp = [ allWinPoly(iBands(iBand),iSectors(iSector),iFrame).outerBorder(2,:) ...
                         allWinPoly(iBands(iBand),iSectors(iSector),iFrame).innerBorder(2,end:-1:1)];
                     
                 if winStati(iBands(iBand),iSectors(iSector),iFrame) > 0
                    fill(xTmp,yTmp,colors{states(iSector,iFrame)},'FaceAlpha',.5,'EdgeColor','k','EdgeAlpha',.5);
                 end
            end            
        end
    end
    
    if useProt && iFrame <= length(protrusion) %#ok<USENS>
        %Draw the protrusion vectors
        quiver(protrusion{iFrame}(:,1),protrusion{iFrame}(:,2),protrusion{iFrame}(:,3),protrusion{iFrame}(:,4),0,'Color','m')             
     end
    
    %Draw the time
    text(10,20,[num2str((iFrame-1)*movieData.timeInterval_s) ' s'],'color','w','FontSize',16)
    
    axis fill,axis image
    set(gca,'color','w');

    MakeQTMovie('addaxes'); 
end

MakeQTMovie('finish');

