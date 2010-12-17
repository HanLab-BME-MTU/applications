function [cellDistGradx,cellDistGrady]=calcCellDistGrad(tracksMatxCord,tracksMatyCord,distGrad,toDoList)
if nargin < 4 || isempty(toDoList);
    toDoList=1:length(distGrad);
end

for frame=toDoList
    DX=distGrad(frame).xmat;
    DY=distGrad(frame).ymat;
    xvec=tracksMatxCord(:,frame);
    yvec=tracksMatyCord(:,frame);
    
    linInd=sub2ind(size(DX),yvec,xvec);
    % sort out possible NaN values:
    badInd=find(isnan(linInd));
    goodInd=find(~isnan(linInd));
    
    cellDistGradx(goodInd,frame)=-DX(linInd(goodInd));
    cellDistGrady(goodInd,frame)=-DY(linInd(goodInd));
    
    cellDistGradx(badInd,frame)=NaN;
    cellDistGrady(badInd,frame)=NaN;
    
    % normalize to 1:
    nVec=sqrt(cellDistGradx(:,frame).^2+cellDistGrady(:,frame).^2);
    cellDistGradx(:,frame)=cellDistGradx(:,frame)./nVec;
    cellDistGrady(:,frame)=cellDistGrady(:,frame)./nVec;
end