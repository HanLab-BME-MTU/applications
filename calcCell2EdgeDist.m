function [cell2edgeDist]=calcCell2EdgeDist(tracksMatxCord,tracksMatyCord,cellDistFromEdge,toDoList)
if nargin < 5 || isempty(toDoList);
    toDoList=1:length(cellDistFromEdge);
end

for frame=toDoList
    D=cellDistFromEdge(frame).mat;
    xvec=tracksMatxCord(:,frame);
    yvec=tracksMatyCord(:,frame);
   
    linInd=sub2ind(size(D),yvec,xvec);
    % sort out possible NaN values:
    badInd=find(isnan(linInd));
    goodInd=find(~isnan(linInd));
    
    cell2edgeDist(goodInd,frame)=D(linInd(goodInd));
    cell2edgeDist(badInd,frame)=NaN; 
end