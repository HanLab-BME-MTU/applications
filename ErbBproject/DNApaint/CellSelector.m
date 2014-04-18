function [CI] = CellSelector(clusterInfo)
%CellSelector takes in clusterInfo and prompts you to select the region
%containing a cell

pnts = vertcat(clusterInfo.pnts);
img = hist3(pnts(:,1:2),'Edges',{[0:max(pnts(:,1))],[0:max(pnts(:,2))]});
pnts = vertcat(clusterInfo.ptClusterCenter);
imshow(log(img),[]);
[BW,xi,yi]= roipoly;

in = inpolygon(pnts(:,1),pnts(:,2),yi,xi);

CI= clusterInfo(in);

end