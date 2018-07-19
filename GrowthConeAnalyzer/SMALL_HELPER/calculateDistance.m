
function [dist,pixIdxTrunc,deltC] = calculateDistance(pixIdx,size,varargin)
% calculateDistance: Small function to calculate the distance along a line of pixels 
% Note this function currently does not perform any smoothing. 
% 
% INPUT 
% pixIdx : rx1 
%
%

% OUTPUT 
% pixIdxTrunc: the input pixels cut-off from a certain length from position 1 
% this can then be used to collect the medial axis transform lengths 

% OUTPUT
% pixIdxTrunc: the input pixels cut-off from a certain length from position 1
% this can then be used to collect the medial axis transform lengths


%% Check Input
ip = inputParser;
ip.addRequired('pixIdx',@(x) isvector(x) || isempty(x) );
ip.addRequired('size',@isvector);
ip.addOptional('mkPlot',0,@isscalar); % 
ip.addParameter('distCutOff',[],@(x) isscalar(x) || isempty(x));
ip.addParameter('pixelSizeMic',0.216, @isscalar); % in microns
ip.parse(pixIdx,size,varargin{:});

distCutOff = ip.Results.distCutOff;
pixelSize = ip.Results.pixelSizeMic;
mkPlot = ip.Results.mkPlot;

[y,x]  = ind2sub(size,pixIdx);

deltX = diff(x);
deltY = diff(y);

delt =  arrayfun(@(i) sqrt(deltX(i)^2+deltY(i)^2),1:length(deltX));
deltC = cumsum(delt);
deltCMic = deltC.*pixelSize;
if ~isempty(distCutOff)
    
    
    pixIdxTrunc= pixIdx(deltCMic<distCutOff);
else
    pixIdxTrunc = [];
    
end



if ~isempty(delt) % ie not a singleton    
    dist = sum(delt);
else
    dist = 0.001 ; % just set singleton distances to 0 for now : note the adjMat will not include any
    % values of 0 so have to set to 0.001
end

dist = dist.*pixelSize; % convert from pixels to um

%% SANITY : Plot the Pixels and the Distances Along
if mkPlot ==1    
    figure
    plot(x,y,'color','r');
    scatter(x,y,10,'r','filled');
    text(x(1),y(1),num2str(dist),'color','g');
end
end






