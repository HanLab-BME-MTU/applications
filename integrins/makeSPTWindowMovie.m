function makeSPTWindowMovie(movieData,sptPropInWindow,propName,iBands,iWindows,iChan,movieName) %#ok<INUSL>

if nargin < 7 || isempty(movieName)
    movieName = ['windows_' propName];
end

if nargin < 6 || isempty(iChan)
    iChan = 1;
end

%Load the windows
disp('Loading Windows...')
load([movieData.windows.directory filesep movieData.windows.fileName])

%Determine number of bands, windows and images
[nBandsTot,nWindowsTot,nImages] = size(allWinPoly);

%load the samples to see if window passed quality control
disp('Checking activity samples...')
try
    load([movieData.activity.directory filesep movieData.activity.fileName]);
catch
    disp('Couldnt load activity samples, using all windows')
end

if exist('allWindowSamples','var') && isfield(allWindowSamples,'nPixels')
    winStati = ~isnan(allWindowSamples.nPixels);
else
    disp('Cannot find window samples, displaying all windows!');
    winStati = ones(nBandsTot,nWindowsTot,nImages); %TEMP TEMP!
end

%Default is to use all windows
if nargin < 5 || isempty(iWindows)
    iWindows = 1:nWindowsTot;
end
%Default is to use all bands
if nargin < 4 || isempty(iBands)
    iBands = 1:nBandsTot;
end

%Determine number of windows/bands requested for tracking
nBands = length(iBands);
nWindows = length(iWindows);

disp('Checking image files...')
imageFileNames = dir([movieData.imageDirectory filesep movieData.channelDirectory{iChan} filesep '*.tif']); %TEMP TEMP

maskFileNames = dir([movieData.masks.directory filesep movieData.masks.channelDirectory{iChan} filesep '*.tif']);

%Load the first image
currImage = imread([movieData.imageDirectory filesep movieData.channelDirectory{iChan} filesep imageFileNames(1).name]);

%Determine size of the image
[imageM,imageN] = size(currImage);

winColors = rand(nBands,nWindows,3);

try
    disp('Checking protrusion vectors....')
    load([movieData.protrusion.directory filesep movieData.protrusion.fileName])
    useProt = true;
catch
    disp('Couldnt load protrusion vectors - not displaying.')
    useProt = false;
end

%% SPT property color-coding

%get current property
eval(['propCurrent = sptPropInWindow.' propName '.values;']);

switch propName
    
    case {'spDensity','f2fDisp','diffCoef','confRad'}
        
        %get values of property to be plotted
        valuesNoNaNs = propCurrent(iBands,iWindows,:); %#ok<NODEF>
        valuesNoNaNs = valuesNoNaNs(:);
        valuesNoNaNs = valuesNoNaNs(~isnan(valuesNoNaNs));
        minValue = prctile(valuesNoNaNs,2);
        maxValue = prctile(valuesNoNaNs,98);
        
    case {'fracUnclass','fracConf','fracBrown','fracDir'}
        
        minValue = 0;
        maxValue = 1;
        
end

%divide range of values into 100 segments
stepSize = (maxValue - minValue) / 100;
thresholdVal = minValue : stepSize: maxValue;
numThresh = length(thresholdVal);

%classify each window in each frame based on its property value
windowClass = zeros(size(allWinPoly));
for iThresh = 1 : numThresh
    windowClass(propCurrent>=thresholdVal(iThresh)) = windowClass(propCurrent>=thresholdVal(iThresh)) + 1;
end

cmap1 = [1 1 1; colormap(jet(101))];
close

%% Make movie

%open a figure window that fills the screen
figure('units','normalized','position',[0 0 1 1])

for iImage = 1:nImages
    
    clf
    
    % START LEFT PANEL
    
    axes('Position',[0 0.1 0.495 0.9]);
    
    hold on
    
    currImage = imread([movieData.imageDirectory filesep movieData.channelDirectory{iChan} filesep imageFileNames(iImage).name]); %TEMP TEMP!!
    
    currImage = cast(currImage,'double');
    
    %Draw the image
    imagesc(currImage),colormap gray, axis image,axis off
    
    if iImage == 1
        
        %Find the mean and standard deviation of first image values
        imMeanL = nanmean(nonzeros(cast(currImage,'double')));
        imStdL = nanstd(nonzeros(cast(currImage,'double')));
    end
    
    %Use this to set the color axis on all frames
    caxis([ (imMeanL - 1.7*imStdL) (imMeanL + 2.5*imStdL)]);
    set(gca,'color','k')
    
    set(gca,'YDir','reverse')
%     %Overlay each requested window on the image
%     for iBand = 1:nBands
%         for iWindow = 1:nWindows
%             
%             if ~isempty(allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder) &&  ~isempty(allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder)
%                 xTmp = [ allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder(1,:) ...
%                     allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder(1,end:-1:1)];
%                 yTmp = [ allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder(2,:) ...
%                     allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder(2,end:-1:1)];
%                 
%                 if winStati(iBands(iBand),iWindows(iWindow),iImage) > 0
%                     fill(xTmp,yTmp,'r','FaceAlpha',0,'EdgeColor','g','EdgeAlpha',0.5);
%                 end
%             end
%         end
%     end
    
    if useProt && iImage <= length(protrusion)
        %Draw the protrusion vectors
        quiver(protrusion{iImage}(:,1),protrusion{iImage}(:,2),protrusion{iImage}(:,3),protrusion{iImage}(:,4),0,'Color','m')
        
    end
    
    %Draw the time
    text(10,20,[num2str((iImage-1)*movieData.timeInterval_s) ' s'],'color','w','FontSize',32)
    
    axis fill,axis image
    set(gca,'color','w');
    
    % END LEFT PANEL
    
    % START RIGHT PANEL
    
    axes('Position',[0.505 0.1 0.495 0.9]);
    
    hold on
    
    currImage = imread([movieData.masks.directory filesep movieData.masks.channelDirectory{iChan} filesep maskFileNames(iImage).name]); %TEMP TEMP!!
    
    currImage = cast(currImage,'double');
    
    %Draw the image
    imagesc(currImage),colormap gray, axis image,axis off
    
    %Use this to set the color axis on all frames
    caxis([0 1]);
    set(gca,'color','k')
    
    set(gca,'YDir','reverse')
    %Overlay each requested window on the image
    for iBand = 1:nBands
        for iWindow = 1:nWindows
            
            windowColor = cmap1(windowClass(iBands(iBand),iWindows(iWindow),iImage)+1,:);
            
            if ~isempty(allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder) &&  ~isempty(allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder)
                xTmp = [ allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder(1,:) ...
                    allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder(1,end:-1:1)];
                yTmp = [ allWinPoly(iBands(iBand),iWindows(iWindow),iImage).outerBorder(2,:) ...
                    allWinPoly(iBands(iBand),iWindows(iWindow),iImage).innerBorder(2,end:-1:1)];
                
                if winStati(iBands(iBand),iWindows(iWindow),iImage) > 0
                    fill(xTmp,yTmp,windowColor,'FaceAlpha',1,'EdgeColor',[0 0 0],'EdgeAlpha',0.5);
                end
            end
        end
    end
    
    if useProt && iImage <= length(protrusion)
        %Draw the protrusion vectors
        quiver(protrusion{iImage}(:,1),protrusion{iImage}(:,2),protrusion{iImage}(:,3),protrusion{iImage}(:,4),0,'Color','m')
        
    end
    
    axis fill,axis image
    set(gca,'color','w');
    
    % END RIGHT PANEL
    
    % START COLOR LEGEND
    
    axes('Position',[0.55 0.04 0.4 0.015],'Ytick',zeros(1,0),...
        'XTick',[thresholdVal(1) thresholdVal(end)],'Fontsize',24);
    
    hold on
    
    for iThresh = 1 : 101
        plot(thresholdVal(iThresh)*[1 1],[1 2],'color',cmap1(iThresh+1,:),'LineWidth',4)
    end
    xlim([thresholdVal(1) thresholdVal(end)])
    title(['Color legend for ' propName],'Fontsize',24)

    % END COLOR LEGEND
    
    if iImage == 1
        MakeQTMovie('start',[movieData.analysisDirectory filesep movieName '.mov'])
        MakeQTMovie('quality',.75)
    end
    MakeQTMovie('addfigure')
    
end

MakeQTMovie('finish')
%movie2avi(windowTrackMovie,[movieData.analysisDirectory filesep movieName],'compression','Cinepak');
%close(h);

