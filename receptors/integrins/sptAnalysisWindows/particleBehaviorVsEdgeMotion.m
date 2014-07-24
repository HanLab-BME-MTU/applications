function particleBehaviorVsEdgeMotion(protPerWindow,sptPropInWindow,...
    prop2analyze,band2analyze,window2analyze,frameRange) %#ok<INUSL>
%PARTICLEBEHAVIORVSEDGEMOTION looks for a relationship between particle behavior and cell edge protrusion activity
%
%SYNOPSIS
%
%INPUT  protPerWindow  : Average protrusion vector per window (from
%                        Hunter).
%       sptPropInWindow: Output of particleBehaviorInWindows.
%       prop2analyze   : Cell array indicating properties to analyze. The
%                        names must exactly match the field names in
%                        sptPropInWindow.
%                        Optional. Default: All properties.
%       band2analyze   : Vector indicating which bands to analyze. Bands
%                        arise from the division of a cell into windows and
%                        run parallel to the cell edge. Band #1 is at the
%                        cell edge.
%                        Optional. Default: 1.
%       window2analyze : Vector indicating which windows to analyze. 
%                        Optional. Default: all windows.
%       frameRange     : Row vector with first and last frame to include
%                        in analysis.
%                        Optional. Default: all frames.
%
%OUTPUT
%
%REMARKS NEEDS UPDATING because of classifyEdgeMotion
%
%
%Khuloud Jaqaman, September 2010

%% Input

if nargin < 3 || isempty(prop2analyze)
    prop2analyze = {'spDensity','f2fDispMag2D','angleMean','angleStd',...
        'f2fDispMagParaDir','f2fDispMagPerpDir','f2fDispSignParaDir',...
        'f2fDispSignPerpDir','f2fDispMagParaProt','f2fDispMagPerpProt',...
        'f2fDispSignParaProt','f2fDispSignPerpProt','ratioDispMagDir',...
        'ratioDispSignDir','ratioDispMagProt','ratioDispSignProt',...
        'asymParam','f2fDispMag2DLin','f2fDispMagParaDirLin',...
        'f2fDispMagPerpDirLin','f2fDispSignParaDirLin',...
        'f2fDispSignPerpDirLin','f2fDispMagParaProtLin',...
        'f2fDispMagPerpProtLin','f2fDispSignParaProtLin',...
        'f2fDispSignPerpProtLin','ratioDispMagDirLin',...
        'ratioDispSignDirLin','ratioDispMagProtLin',...
        'ratioDispSignProtLin','asymParamLin',...
        'fracUnclass','fracLin','fracIso','fracIsoUnclass',...
        'fracConf','fracBrown','fracDir','diffCoef','confRad'};
end

%get bands to analyze
if nargin < 4 || isempty(band2analyze)
    band2analyze = 1;
end

%calculate number of properties and number of bands
numProp2analyze = length(prop2analyze);
numBands = length(band2analyze);

%get protrusion vector matrix
protNormVecMag = protPerWindow.avgNormal(:,1:end-1);

%get number of windows and time points
[numWindows,numTP] = size(protNormVecMag);

%get which windows to analyze and which to ignore
if nargin < 5 || isempty(window2analyze)
    window2analyze = 1 : numWindows;
end
window2ignore = setdiff(1:numWindows,window2analyze);

%get which frames to include in analysis
if nargin < 6 || isempty(frameRange)
    frameRange = [1 numTP];
else
    frameRange(1) = max(frameRange(1),1);
    frameRange(2) = min(frameRange(2),numTP);
end
numFramesAnalysis = diff(frameRange) + 1;

%define extended color map
cmapExt = [[0 0 0]; colormap];
close

%% Analysis

%classify edge motion activity per window
windowMotionType = classifyEdgeMotion(protPerWindow);
windowMotionType = windowMotionType(:,1:end-1);

%get the indices of windows with the different motion types
indxProtrude = find(windowMotionType==1);
indxRetract = find(windowMotionType==-1);
indxPause = find(windowMotionType==0);

%go over requested particle properties
for iProp = 1 : numProp2analyze
    
    %get particle property name
    propName = prop2analyze{iProp};
    
    %get current particle property
    eval(['propCurrent = sptPropInWindow.' propName '.values;']);
    
    for iBand = band2analyze
        
        propCurrentBand = squeeze(propCurrent(iBand,:,:));
        
        %make a scatter plot of property vs. protrusion normal
        %color coding: green = protrusion, red = retraction, blue = pause
        figure, hold on
        plot(protNormVecMag(indxProtrude),propCurrentBand(indxProtrude),'g+')
        plot(protNormVecMag(indxRetract),propCurrentBand(indxRetract),'rx')
        plot(protNormVecMag(indxPause),propCurrentBand(indxPause),'b.')
        
    end
    
end %(for iProp = 1 : numProp2analyze)












