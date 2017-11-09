function movieInfo= labelToMovieInfo(label,vol)
[feats,nFeats] = bwlabeln(label);
featsProp = regionprops(feats,vol,'Area','WeightedCentroid','MeanIntensity','MaxIntensity','PixelValues');

% centroid coordinates with 0.5 uncertainties
tmp = vertcat(featsProp.WeightedCentroid)-1;
xCoord = [tmp(:,1) 0.5*ones(nFeats,1)]; 
yCoord = [tmp(:,2) 0.5*ones(nFeats,1)]; 
zCoord = [tmp(:,3) 0.5*ones(nFeats,1)];
amp=[vertcat(featsProp.MaxIntensity) 0.5*ones(nFeats,1)];

% u-track formating
movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
movieInfo.xCoord= xCoord;movieInfo.yCoord=yCoord;movieInfo.zCoord=zCoord;
movieInfo.amp=amp;
movieInfo.int=amp;
