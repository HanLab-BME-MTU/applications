function [densityMeasurement]=perfDensityMeasurement(cellDistFromEdge,tracksMatxCord,tracksMatyCord,binPix,toDoList)
if nargin < 5 || isempty(toDoList);
    toDoList=1:length(cellDistFromEdge);
end

tracksMatxCord=round(tracksMatxCord);
tracksMatyCord=round(tracksMatyCord);

for frame=toDoList
    k=1;
    maxC2ED=max(cellDistFromEdge(frame).mat(:));
    while (k-1)*binPix<=maxC2ED
        binMask=(cellDistFromEdge(frame).mat>(k-1)*binPix & cellDistFromEdge(frame).mat<=k*binPix);
        binArea=sum(binMask(:));
        % imagesc(binMask);
        % count the number of cells in each bin:
        xvec=tracksMatxCord(:,frame);
        yvec=tracksMatyCord(:,frame);
    
        % find the ones that are located within the mask:
        linInd=sub2ind(size(binMask),yvec,xvec);
        % sort out possible NaN values:
        linInd(isnan(linInd))=[];
        inNuc =linInd(binMask(linInd)==1);
       
        densityMeasurement.cells(k,frame)  =numel(inNuc);
        densityMeasurement.area(k,frame)   =binArea;
        densityMeasurement.density(k,frame)=numel(inNuc)/binArea;
        densityMeasurement.binPix          =binPix;
        k=k+1;
    end
    densityMeasurement.density(densityMeasurement.cells==0)=NaN;
end