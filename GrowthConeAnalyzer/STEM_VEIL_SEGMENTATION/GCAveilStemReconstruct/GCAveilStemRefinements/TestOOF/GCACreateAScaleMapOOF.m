function [scalemap] = GCACreateAScaleMapOOF(OOFDir,MD,frame)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



files = searchFiles('.tif',[],OOFDir,0,'all',1); 


  ires =   cellfun(@(x) double(imread(x)),files,'uniformoutput',0); 
 
%theta = itheta{1};

scaleMap = ones(size(ires{1}));

res = ires{1}; 
for si = 2:numel(ires) 
    idx = ires{si}>res;
    res(idx) = ires{si}(idx);
    
    scaleMap(idx) = si;
end

% for now just load the NMS from previous frame - see fo you can call some
% of these in from C etc. 

GCARefineThinVeilStemWithActiveContoursMovie(MD); 

% load the currrent img 


% just reintegrate here. 
 toLoad  = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'IIIa_veilStem_refine_ActiveContour' filesep 'veilStem.mat'];
load(toLoad); 

% get the NMS for now from the old steerable filters for now (should be relatively similar) 
NMSSteerFilt = veilStem.maxNMSInt; 

%% PLOT the SCALES FROM THE OOF 
scales = (NMSSteerFilt>0).*scaleMap;
[ny,nx] = size(img); 
setFigure(ny,nx,'on'); 
cmap = colormap(jet(10));
imshow(-img,[]) ;
hold on
for iScale = 1:length(ip.Results.RidgeScalesAC)
    idx = find(scales == iScale);
    [cy,cx]  = ind2sub([ny,nx],idx);
    scatter(cx,cy,10,cmap(iScale,:),'filled');
end

end 



