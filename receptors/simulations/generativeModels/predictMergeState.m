function trajCorr = predictMergeState(rTraj, gTraj, winSize, lagFrames)

rTraj = squeeze(rTraj);
gTraj = squeeze(gTraj);

[numFrames, dim] = size(rTraj);
if numFrames < dim
    rTraj = rTraj';
    gTraj = gTraj';
    [numFrames, ~] = size(rTraj);
end

boxSize = 10; %in micrometers

diffConst = 0.1;
dT = 0.01;
separations = squeeze(sqrt(sum((rTraj-gTraj).^2,2)));
trajCorr = zeros(1,numFrames);
distProb = zeros(1,numFrames);

trajCorr((winSize+lagFrames):(numFrames-winSize-lagFrames)) = ...
    trajCorrelation(rTraj(:,:),gTraj(:,:),(winSize+lagFrames):(numFrames-winSize-lagFrames),winSize,lagFrames);
distProb((lagFrames+1):numFrames) = ...
    probMergedByDist(boxSize,diffConst,dT,separations((lagFrames+1):numFrames),separations(1:(numFrames-lagFrames)));

status = (trajCorr>.3);

%{
for boundT = 1:numFrames
    shape = ones(winSize+boundT+1,1);
    winLag = winSize+lagT-1;
    if boundT > winLag
        for i = (1+winLag):-1:(ceil(winLag/2)+2)
            shape(i) = shape(i+1) - 0.5/floor(winLag/2);
        end
        for i = (boundT+2):(boundT+1+ceil(winLag/2))
            shape(i) = shape(i-1) - 0.5/ceil(winSize/2);
        end
        for i = (ceil(winLag/2)+1):2
            shape(i) = shape(i+1)*5/7;
        end
        for i = (boundT+2+ceil(winLag/2)):(length(shape)-1)
            shape(i+ = shape(i-1)*5/7;
        end
        shape(1) = 0;
        shape(end) = 0;
    end
end
%}


%{
minWidth = winSize;

for i = round(minWidth/2):(numFrames-minWidth+round(minWidth/2))
    if sum(status((i-round(minWidth/2)+1):(i-round(minWidth/2)+minWidth))) < (minWidth*.9)
        newStatus(i) = 0;
    end
end

status = newStatus;
%}