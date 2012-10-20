function quickView3DMovieROI(movieData,iROI)


if nargin < 2 || isempty(iROI)
    iROI = 1;
end

p.ChannelIndex = 1;
p.nSnapshots = 3;

if numel(movieData.rois_) < iROI
    error('Invalid ROI number or movie has no ROIs!')
end

roiMask = movieData.rois_(iROI).getROIMask;
roiMaskBord = bwboundaries(sum(roiMask,3)>0);
roiMaskBord = roiMaskBord{1};

iSnap = round(linspace(1,movieData.nFrames_,p.nSnapshots));

imNames = movieData.getImageFileNames(p.ChannelIndex);
imDir = movieData.getChannelPaths(p.ChannelIndex);

fsFigure(.75);
for j = 1:p.nSnapshots
    
    subplot(1,p.nSnapshots,j);
    
    currIm = stackRead([imDir{1} filesep imNames{1}{iSnap(j)}]);    
    currIm = max(currIm,[],3);
    
    imagesc(currIm);axis image,colormap gray,axis off    
    hold on,saturateImageColormap(gca,3);
    
    plot(roiMaskBord(:,2),roiMaskBord(:,1),'r');
    
end
    
    

