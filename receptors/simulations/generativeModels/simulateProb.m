function corrMerged = simulateProb(rTraj,gTraj,locError,winSize,lagT)

%% Set up
N = 1000; %number of simulations to run

%get trajectories into correct orientation
rTraj = squeeze(rTraj);
gTraj = squeeze(gTraj);
[numFrames, dim] = size(rTraj);
if numFrames < dim
    rTraj = rTraj';
    gTraj = gTraj';
    [numFrames, ~] = size(rTraj);
end

%calculate separations between two trajectories
separations = sqrt(sum((rTraj-gTraj).^2,2));

%% Simulations

%gets N different tracks of possible green trajectories if green and red
%were bound for the whole time
tracksMerged = repmat(rTraj,1,1,N) + randn(numFrames,2,N)*locError*sqrt(2);

%% Correlation

corrMerged = zeros(N,numFrames);

%for each of the trajectories calculate correlation statistic
for n = 1:N
    range = (1+round((winSize+lagT)/2)):(numFrames-round(winSize/2)-lagT-1);
    range = range(separations(range)<4*locError);
    corrMerged(n,range) = trajCorrelation(gTraj,tracksMerged(:,:,n),range,winSize,lagT);
end

%I wanted to massage the data; first I take only the best half of the
%simulations
corrMerged = sort(corrMerged,1);
corrMerged = mean(corrMerged(1:round(N/2),:),1);

%this helps to get rid of the noise by removing entries with negative
%correlations and then applying a median filter
corrMerged = corrMerged.*(corrMerged>0);
corrMerged = medianFilter(corrMerged,3);


function newCorr = medianFilter(corr,win)

[~,numFrames] = size(corr);
newCorr = zeros(size(corr));
for i = (floor(win/2)+1):(numFrames-ceil(win/2))
    newCorr(:,i) = median(corr(:,(i-floor(win/2)):(i-floor(win/2)+win)),2);
end