function [corr]=velCorrFunc(xPos,yPos,xVel,yVel,dR,ptsId)
% Calculates the velocity correlation function in 2D.

if nargin < 6 || isempty(ptsId)
    ptsId=ones(size(xPos));
end

% make sure that all entries are not NaN, that is take only those lines
% into account where all values exist.
nanCheck = ~(isnan(xPos) | isnan(yPos) | isnan(xVel) | isnan(yVel) | isnan(ptsId));
xPos = xPos(nanCheck);
yPos = yPos(nanCheck);
xVel = xVel(nanCheck);
yVel = yVel(nanCheck);
ptsId= ptsId(nanCheck);

% now convert the 0/1 list to indices:
ptsId=find(ptsId==1);
% convert to row vector:
ptsId=(ptsId(:))';

% calculate the cos(alpha) between the velocity and the x-axis:
magVel   = sqrt(xVel.^2+yVel.^2);
cosAlpha = xVel./magVel;
sinAlpha = yVel./magVel;

% initialize R:
% maxDist=sqrt((max(xPos)-min(xPos))^2+(max(yPos)-min(yPos))^2);
% maxBin =ceil(maxDist/dR);
R=[];
% calculate the correlation function for all points identified by ptsId.
k=1;
for pt=ptsId   
%     text=[' Point: ',num2str(k,'%05.f'),' of ',num2str(numel(ptsId))];
%     progressText(k/numel(ptsId),text);    
    px = xPos(pt);
    py = yPos(pt);
    vx = xVel(pt);
    vy = yVel(pt);
    ca = cosAlpha(pt);
    sa = sinAlpha(pt);
    
    % calculate the distance from all points in the set to the current
    % central point pt:
    pts2ptR = sqrt((px-xPos).^2+(py-yPos).^2);
    
    % calculate cos(DeltaAlpha) where DeltaAlpha is the angle between the
    % velocities and the current pt. Here we use
    % cos(a1-a2)=cos(a1)*cos(a2)+sin(a1)*sin(a2);
    % Where a1 and a2 are the angles with respect to the x-axis:
    cosda = cosAlpha*ca+sinAlpha*sa;
    
    
    % Bin the points in the set according to the distance to the current
    % central point pt:
    maxR      = max(pts2ptR);
    % make sure that the autocorrelation gets an own bin! This is important
    % if one wants to normalize the correlation function, such that the
    % autocorrelation == 1.
    [n,binIDs] = histc(pts2ptR,[0,eps:dR:maxR+dR]);
    
    for bin=1:max(binIDs)
        % remove NaNs?!
        checkVec=(binIDs==bin);
        
        if n(bin)~=nansum(checkVec)
            display(['Sum not correct: ',num2str(bin)]);
        end
        
        if length(R)<bin
            R(bin).Ntot=0;
        end
        
        ind=(R(bin).Ntot+1):(R(bin).Ntot+n(bin));
        % store as column vectors:
        R(bin).pts2ptR(ind,1) =  pts2ptR(checkVec);
        R(bin).xVelPts(ind,1) =     xVel(checkVec);
        R(bin).yVelPts(ind,1) =     yVel(checkVec);   
        R(bin).xVelCtr(ind,1) = repmat(vx,n(bin),1);
        R(bin).yVelCtr(ind,1) = repmat(vy,n(bin),1);
        R(bin).cosda(ind,1)   = cosda(checkVec);
        % increase the total number of correlation pairs:
        R(bin).Ntot           = R(bin).Ntot+n(bin);
    end
    
    %remove the point from the list! (otherwise we count each pair twice)
    xPos(pt)    = NaN;
    yPos(pt)    = NaN;
    xVel(pt)    = NaN;
    yVel(pt)    = NaN;
    magVel(pt)  = NaN;
    cosAlpha(pt)= NaN;
    sinAlpha(pt)= NaN;

    k=k+1;
end
% calculate the correlations:
for bin=1:length(R)
    corr.RMean(bin)    = nanmean(R(bin).pts2ptR);
    corr.RSTD(bin)     =  nanstd(R(bin).pts2ptR);
    corr.funcVel(bin)  = corrFunc(R(bin).xVelCtr,R(bin).xVelPts,R(bin).yVelCtr,R(bin).yVelPts);
    corr.funcVelx(bin) = corrFunc(R(bin).xVelCtr,R(bin).xVelPts);
    corr.funcVely(bin) = corrFunc(R(bin).yVelCtr,R(bin).yVelPts);
    corr.Ntot(bin)     = R(bin).Ntot;
    % The cosa needs a special treatment, since here the caclulation could
    % introduce additional NaNs (if vel=0). These positions have to be
    % canceled out in the distance vector:
    corr.RMean4CosaMean(bin) = mean(R(bin).pts2ptR(~isnan(R(bin).cosda)));
    corr.RMean4CosaSTD(bin)  =  std(R(bin).pts2ptR(~isnan(R(bin).cosda)));
    corr.cosaMean(bin) = nanmean(R(bin).cosda);
    corr.cosaSTD(bin)  =  nanstd(R(bin).cosda);    
end
corr.R=R;

% Normalize such that the autocorrelation is 1:
% corr.funcVel  = corr.funcVel /corr.funcVel(1);
% corr.funcVelx = corr.funcVelx/corr.funcVelx(1);
% corr.funcVely = corr.funcVely/corr.funcVely(1);

return;


% test this function:
xvec=linspace(-10,10,1001);
yvec=linspace(-10,10,1001);
dR = 0.1;

lx=0.5;
ly=2;

[xmat,ymat]=meshgrid(xvec,yvec);

xVel=exp(-(1/lx)*sqrt((xmat.^2+ymat.^2)));
yVel=exp(-(1/ly)*sqrt((xmat.^2+ymat.^2)));

xPos=xmat(:);
yPos=ymat(:);
xVel=xVel(:);
yVel=yVel(:);

% add random value:
%xVel=xVel+0.05*randn(size(xVel));
%yVel=yVel+0.05*randn(size(yVel)); 

ptsId=(xPos==0 & yPos==0);

%ptsId=[];
[corr]=velCorrFunc(xmat(:),ymat(:),xVel(:),yVel(:),dR,ptsId);

%figure()
%quiver(xPos,yPos,xVel,yVel)


figure()
plot(corr.RMean,corr.funcVelx,'or');
hold on
plot(corr.RMean,corr.funcVely,'sg');
plot(corr.RMean,corr.funcVel ,'*b');
hold off

figure()
errorbar(corr.RMean4CosaMean,corr.cosaMean,corr.cosaSTD,'sk');


[ux,sigmaUx,goodRowsx]=robustExponentialFit(corr.RMean, corr.funcVelx, 1,1);
[uy,sigmaUy,goodRowsy]=robustExponentialFit(corr.RMean, corr.funcVely, 1,1);
[u,sigmaU,goodRows]   =robustExponentialFit(corr.RMean, corr.funcVel, 1,1);

ux
uy
u