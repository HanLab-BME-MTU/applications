function [ output_args ] = GCAVisualsColorCodeByCurvature(filoInfo,varargin)
%GCAVisualsColorCodeByCurvature: 
% Small visualization function to plot curvature. 
%
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
ip.addRequired('filoInfo');
ip.addParameter('cMapLimits',[]); 
ip.addParameter('filoFilterSet',[]); 
ip.addParameter('pix2Micron',[]); 

ip.parse(filoInfo,varargin{:});
%%
%set the defaults for filoFilterSet
if isempty(ip.Results.filoFilterSet);
    filoFilterSet = true(length(filoInfo),1);
else
    filoFilterSet = ip.Results.filoFilterSet;
end


filoInfo = filoInfo(filoFilterSet(:,1)); 

filoCurvsAll =   abs(vertcat(filoInfo(:).Ext_FiloCurvIndVals));
xy = vertcat(filoInfo(:).Ext_coordsXY_SplineFit);

if ~isempty(ip.Results.pix2Micron) 
   filoCurvsAll = filoCurvsAll*(1/ip.Results.pix2Micron); % convert from 1/pixels to 1/um  
end 


if isempty(ip.Results.cMapLimits);
    cMapLimits(1) = min(filoCurvsAll);
    cMapLimits(2) = max(filoCurvsAll);
else 
    cMapLimits(1) = ip.Results.cMapLimits(1); 
    cMapLimits(2) = ip.Results.cMapLimits(2); 
end


% create distance mapper
cMapLength=128; cMap=jet(cMapLength);

mapper=linspace(cMapLimits(1),cMapLimits(2),cMapLength)';

D=createDistanceMatrix(filoCurvsAll,mapper);
[sD,idxCMap]=sort(abs(D),2);

for k=1:cMapLength
    %                 idxCand = EPreFilt(idxCMap(:,1) == k,2);
    %                 idxSeed = EPreFilt(idxCMap(:,1)==k,1);
    %                 for iEdge = 1:length(idxCand) % some can have the same color
    
    scatter(xy(idxCMap(:,1)==k,1),xy(idxCMap(:,1) == k,2),5,...
        cMap(k,:),'filled');
    %end
end

end

