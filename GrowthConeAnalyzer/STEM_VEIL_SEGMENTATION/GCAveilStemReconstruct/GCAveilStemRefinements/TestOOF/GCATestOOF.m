function [scalemap] = GCATestOOF(OOFDir,MD,frame)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

OOFFigDir = [OOFDir filesep 'OOFFigs']; 
if ~isdir(OOFFigDir)
    mkdir(OOFFigDir)
end 
%% load the OOF Filter Response and integrate
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

outputDir = [OOFDir filesep 'SteerableFilterFigs']; 
if ~isdir(outputDir) 
    mkdir(outputDir) 
end 

[NMSSteerFilt, img ]= GCAGetSteerableFilterScaleOverlaysMovie(MD,'OutputDirectory',outputDir,'frame',frame); 

% load the currrent img 

close all 


%% PLOT the SCALES FROM THE OOF 
scales = (NMSSteerFilt>0).*scaleMap;
[ny,nx] = size(img); 
setFigure(ny,nx,'on'); 
cmap = colormap(jet(10));
imshow(-img,[]) ;
hold on
for iScale = 1:length(1:10)
    idx = find(scales == iScale);
    [cy,cx]  = ind2sub([ny,nx],idx);
    scatter(cx,cy,10,cmap(iScale,:),'filled');
end

% setFigure(nx,ny,'on');
% % get the values of the response in the NMS
%  values = res(NMSSteerFilt~=0);
%  NMSApproxOOF = res.*NMSSteerFilt; 
% 
%  cutoff = prctile(values,95); 
%  imshow(res.*NMSSteerFilt,[cutoff,max(values(:))]);
%  
%  setAxis('on')
% 
% [count  ] = hist(values,100);
%     hist(values,100);
%     hold on
%     line([cutoff cutoff],[0 max(count)],'color','r','Linewidth',2);
%     xlabel('Ridge Response NMS');
%     ylabel('Count');
%     title(['Red Line = The 75th percentile']);
%     set(gcf,'Color',[1,1,1]);
% 
% 
%  ridgeCand = (NMSApproxOOF>cutoff) & scaleMap>3;
%  setFigure(nx,ny,'on'); 
%  imshow(res.*ridgeCand);  
%   [respNMSMean,respNMSSTD]   = fitGaussianModeToPDF(values); 
% cutoffTrueResponse = respNMSMean+3*respNMSSTD;
saveas(gcf,[OOFFigDir filesep num2str(frame,'%03d') '.png']); 
saveas(gcf,[OOFFigDir filesep num2str(frame,'%03d'),'.fig']); 
end 



