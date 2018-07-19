function [thicknessValues,TSFig] = GCAAnalysisExtract_VeilStemThickness(neuriteLongPathIndicesC,veilStemMaskC,varargin)
%% 

%% Output: 
% thicknessValues: nx1 vector of distance transform values along neurite 
% longest path: from tip of
% neurite to distCutOffForThickness (in um) 
%% CheckInput
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('neuriteLongPathIndicesC'); 
ip.addRequired('veilStemMaskC'); 

ip.addParameter('distCutOffForThickness',20); % in um
ip.addParameter('pixelSizeMic',0.216); % pixel size in um
ip.addParameter('visual',false); 

ip.parse(neuriteLongPathIndicesC,veilStemMaskC,varargin{:});
p = ip.Results;

%% 
[ny,nx]= size(veilStemMaskC); 

distCutOff = ip.Results.distCutOffForThickness;  

[~,pixMeasure] = calculateDistance(neuriteLongPathIndicesC,[ny,nx],'distCutOff',distCutOff,'pixelSizeMic',ip.Results.pixelSizeMic); 

distTrans = bwdist(~veilStemMaskC); 
distTransMic = distTrans.*ip.Results.pixelSizeMic; 
thicknessValues = distTransMic(pixMeasure); 
%thicknessValues.*

%% Make the Visuals : Eventually want to make this such that read in image and overlays the colormap 
if ip.Results.visual == true 
    
    TSFig(1).h = setFigure(nx,ny,'off');
    pixMask = false([ny,nx]);
    pixMask(pixMeasure) = true;
    roiYX =  bwboundaries(veilStemMaskC);
    
    imagesc(distTransMic.*pixMask,[0,5]);
    hold on
    cellfun(@(x) plot(x(:,2),x(:,1),'w'),roiYX);
    thickestPt = pixMeasure(thicknessValues == max(thicknessValues));
    [yMax,xMax] = ind2sub([ny,nx],thickestPt);
    distAtMax = thicknessValues(thicknessValues == max(thicknessValues)); 
    distAtMaxPix = distAtMax/ip.Results.pixelSizeMic; % convert back to pixels
    % find the max point
    gcaCircles(xMax,yMax,distAtMaxPix,'facecolor','none','edgeColor','w');
    scatter(xMax,yMax,10,'w','filled')
    
    text(5,10,'Values used for thickness measurement','color','w','FontSize',10); 
    text(5,20,'Circle plots the radius of the thickest point','color','w','FontSize',10); 
    text(5,30,['Distance From Leading Protrusion = ' num2str(p.distCutOffForThickness) 'um'],'color','w','FontSize',10); 
    forScaleBar = 10/ip.Results.pixelSizeMic; % 10 um scale bar 
    plotScaleBar(forScaleBar,'height', 2, 'Color',[1 1 1],'Location', 'SouthEast');
    colorbar
    
end 


%%
end

