function quickView3DMovieROI(movieData,iROI)


nROIMax = numel(movieData.rois_);

if nargin < 2 || isempty(iROI)
    iROI = 1:nROIMax;
end

nROI = numel(iROI);

p.ChannelIndex = 1;
p.nSnapshots = 1;

if  any(iROI>nROIMax)
    error('Invalid ROI number or movie has no ROIs!')
end


if p.nSnapshots > 1
    iSnap = round(linspace(1,movieData.nFrames_,p.nSnapshots));
else
    iSnap = 1;
end

imNames = movieData.getImageFileNames(p.ChannelIndex);
imDir = movieData.getChannelPaths(p.ChannelIndex);

fsFigure(.75);

roiCols = jet(nROI);

for j = 1:p.nSnapshots
    
    
    subplot(1,p.nSnapshots,j);
    
    currIm = stackRead([imDir{1} filesep imNames{1}{iSnap(j)}]);    
    currIm = max(currIm,[],3);
    
    imagesc(currIm);axis image,colormap gray,axis off    
    hold on,saturateImageColormap(gca,3);
    
    for k = 1:nROI
    
        roiMask = movieData.rois_(iROI(k)).getROIMask;
        roiMaskBord = bwboundaries(sum(roiMask,3)>0);
        roiMaskBord = roiMaskBord{1};
        plot(roiMaskBord(:,2),roiMaskBord(:,1),'color',roiCols(k,:));
        
    end
    
    
    
    
end
    
    

