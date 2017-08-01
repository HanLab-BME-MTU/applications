function densities=estimateDensity(detections,radius)
%%
densities=cell(1,length(detections));
volume=pi*(radius^3)*4/3;
oneMicron=pi*(10^3)*4/3;
for fIdx=1:length(detections)
    pos=detections(fIdx).getPosMatrix();
    closePos=KDTreeBallQuery(pos,pos,radius);
    densities{fIdx}=(cellfun(@(c) oneMicron*length(c)/volume, closePos));
end
        