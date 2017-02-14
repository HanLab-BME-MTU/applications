function [detections,sphericalCoord,inliers]=tracks2detections(tracks,detections,poleDetection,scales)

% Translate these changes in the detection structure
for kIdx=1:length(tracks)
    tr=tracks(kIdx);
    for tIdx=1:tr.lifetime
        f=tr.f(tIdx);
        if(tr.tracksFeatIndxCG(tIdx)>0)
            detections(f).xCoord(tr.tracksFeatIndxCG(tIdx),1)=tr.x(tIdx);
            detections(f).yCoord(tr.tracksFeatIndxCG(tIdx),1)=tr.y(tIdx);
            detections(f).zCoord(tr.tracksFeatIndxCG(tIdx),1)=tr.z(tIdx);
        end
   end
end

% Build the associated polar coordinate
[~,sphericalCoord,~,inliers]=poleDist(poleDetection,detections,'anisotropy',scales,'angleRef','poles');
