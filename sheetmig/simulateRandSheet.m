function [xPos yPos]=simulateRandSheet(xLimit,yLimit,nCells,numTpts,Dx,Dy,xDrift,yDrift)
initPos=[rand(nCells,1)*(xLimit(2)-xLimit(1))+xLimit(1),rand(nCells,1)*(yLimit(2)-yLimit(1))+yLimit(1)];
xPos=zeros(nCells,numTpts);
yPos=zeros(nCells,numTpts);
for iCell=1:nCells
    traj=randomWalk(numTpts-1,initPos(iCell,:),Dx,Dy,xDrift,yDrift);
    xPos(iCell,:)=traj(:,1);
    yPos(iCell,:)=traj(:,2);
end

% plot(xCord',yCord')