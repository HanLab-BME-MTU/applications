function [PL] = CellSelectorRect(PointList)
%CellSelector takes in clusterInfo and prompts you to select the region
%containing a cell

n= numel(PointList);
pnts=[];
for i=1:n
    pnts=vertcat(pnts,PointList{i}.pnts);
end

img = hist3(pnts(:,1:2),'Edges',{[0:max(pnts(:,1))],[0:max(pnts(:,2))]});
h=figure;
imshow(img,[0,max(img(:)*.01)]);
colormap('hot');
rect = getrect;

%convert rect into a polygon
xi = [rect(1),rect(1),rect(1)+rect(3),rect(1)+rect(3),rect(1)];
yi = [rect(2),rect(2)+rect(4),rect(2)+rect(4),rect(2),rect(2)];



for i=1:n
    pnts = PointList{i}.pnts;
    in = inpolygon(pnts(:,1),pnts(:,2),yi,xi);
    PointList{i}.pnts = pnts(in,:);
    PointList{i}.fullData = PointList{i}.fullData(in,:); 
end

close(h);

PL = PointList;

end