function [pMean,pSEM95,Ntot,pxMean,pxSEM95,pyMean,pySEM95]=calcPersistence(xPos,yPos)
[~,numFrames]=size(xPos);

for frame=1:numFrames-1
    % Length of each single segment
    dx(:,frame) = xPos(:,frame+1)-xPos(:,frame);
    dy(:,frame) = yPos(:,frame+1)-yPos(:,frame);    
    dL(:,frame) = sqrt(dx(:,frame).^2+dy(:,frame).^2);
    
    % Length of the direct connection:
    dirX(:,frame) = xPos(:,frame+1)-xPos(:,1);
    dirY(:,frame) = yPos(:,frame+1)-yPos(:,1);
    dirL(:,frame) = sqrt(dirX(:,frame).^2+dirY(:,frame).^2);
end

% The persistence of the full movement (x&y): 
cumL = cumsum(dL,2);
p    =  dirL./cumL;

pMean  = nanmean(p,1);
Ntot   = sum(~isnan(cumL),1);
pSEM95 = facSEMtoSEM95*nanstd(p,[],1)./sqrt(Ntot);

if nargout>3
    % The persistence of the movement in x: 
    cumLx = cumsum(abs(dx),2);
    px    =  abs(dirX)./cumLx;

    pxMean  = nanmean(px,1);
    pxSEM95 = facSEMtoSEM95*nanstd(px,[],1)./sqrt(Ntot);

    % The persistence of the movement in y: 
    cumLy = cumsum(abs(dy),2);
    py    =  abs(dirY)./cumLy;

    pyMean  = nanmean(py,1);
    pySEM95 = facSEMtoSEM95*nanstd(py,[],1)./sqrt(Ntot);
end