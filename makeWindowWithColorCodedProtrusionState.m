function makeWindowWithColorCodedProtrusionState(movieData,...
                                                 iBands,...
                                                 iWindows,...
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
[nBandsTot,nWindowsTot,nImages] = size(allWinPoly);

%Load protrusion samples
disp('Loading protrusion samples...');
load(fullfile(movieData.protrusion.directory, movieData.protrusion.samples.fileName));
if ~isfield(protrusionSamples,'classNames') || ~isfield(protrusionSamples,'classes')
    error('Classification of edge velocity has not been performed.');
end

classes = protrusionSamples.classes;

if exist('allWindowSamples','var') && isfield(allWindowSamples,'nPixels')
    winStati = ~isnan(allWindowSamples.nPixels);
else
    disp('Cannot find window samples, displaying all windows!');
    winStati = ones(nBandsTot,nWindowsTot,nImages); %TEMP TEMP!
end

%Default is to use all windows
if nargin < 3 || isempty(iWindows)
    iWindows = 1:nWindowsTot;
end
%Default is to use all bands
if nargin < 2 || isempty(iBands)
    iBands = 1:nBandsTot;
end

%Determine number of windows/bands requested for tracking
nBands = length(iBands);
nWindows = length(iWindows);

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

for iImage = 1:nImages-1
    
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
    for iWindow = 1:nWindows
        switch classes(iWindow,iImage)
            case 0
                color = 'y'; % pause
            case 1
                color = 'b'; % protrusion
            case 2
                color = 'r'; % retraction
        end
        
        for iBand = 1:nBands            
            if ~isempty(allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder) &&  ~isempty(allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder)
                xTmp = [ allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder(1,:) ...
                         allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder(1,end:-1:1)];
                yTmp = [ allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder(2,:) ...
                         allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder(2,end:-1:1)];
                     
                 if winStati(iBands(iBand),iWindows(iWindow),iImage) > 0
                    fill(xTmp,yTmp,color,'FaceAlpha',.5,'EdgeColor','k','EdgeAlpha',.5);
                 end
            end            
        end
    end
    
    if useProt && iImage <= length(protrusion) %#ok<USENS>
        %Draw the protrusion vectors
        quiver(protrusion{iImage}(:,1),protrusion{iImage}(:,2),protrusion{iImage}(:,3),protrusion{iImage}(:,4),0,'Color','m')
                         
     end
    
    %Draw the time
    text(10,20,[num2str((iImage-1)*movieData.timeInterval_s) ' s'],'color','w','FontSize',16)
    
    axis fill,axis image
    set(gca,'color','w');
    
    if iImage == 1
        MakeQTMovie('start',[movieData.analysisDirectory filesep movieName '.mov']);
        MakeQTMovie('quality',.75);
    end    
    MakeQTMovie('addaxes');
    
end

MakeQTMovie('finish');

