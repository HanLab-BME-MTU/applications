function receptorInfo = getClusterTwoColors(receptTrajR,receptTrajG,plotR,plotG,param)

%% Inputs

%must have 2 inputs
if nargin < 2
    disp('Please enter two data sets.');
    return;
end

if nargin == 4
    param = [];
end

%inputs must be of same time length
if  length(receptTrajR) ~= length(receptTrajG)
    disp('Please use data with the same number of frames.')
    return;
end

%check number of dimensions in inputs
if ndims(receptTrajR)~=3 || ndims(receptTrajG)~=3
   disp('Please give inputs with 3 dimensions.');
   return;
end

%time between frames, default 0.01
dT = .01;
if isfield(param,'dT')
    dT = param.dT;
end

%localization error, default 0
locError = 0;
if isfield(param,'locError')
    locError = param.locError;
end

%standard deviation of separation of two receptors in dimer, default 0.1
sigmaDimer = 0.1;
if isfield(param,'sigmaDimer')
    sigmaDimer = param.sigmaDimer;
end

%diffusion constant of free receptor, default 0.1
D = 0.1;
if isfield(param,'diffConst')
    D = param.diffConst;
end

%side length of field of observation; if not present then it guesses based
%on the largest location
if isfield(param,'fieldSize')
    fieldSize = param.fieldSize;
else
    fieldSize = max(max(receptTrajR(:)),max(receptTrajG(:)));
end
    
%isolate x- and y-coordinates of the red and green trajectories
rXCoord = squeeze(receptTrajR(:,1,:));
rYCoord = squeeze(receptTrajR(:,2,:));
gXCoord = squeeze(receptTrajG(:,1,:));
gYCoord = squeeze(receptTrajG(:,2,:));

%get number of frames and green/red receptors
[numG,numFrames] = size(gXCoord);
[numR,~] = size(rXCoord);

%% Main Body

%get displacement vectors

rDispVectors = zeros(numR,2,numFrames);
gDispVectors = zeros(numG,2,numFrames);

for t = 2:numFrames
    rDispVectors(:,:,t-1) = receptTrajR(:,:,t) - receptTrajR(:,:,t-1);
    gDispVectors(:,:,t-1) = receptTrajG(:,:,t) - receptTrajG(:,:,t-1);
end

%get displacements
receptSep = zeros(numR,numG,numFrames);
for r = 1:numR
    for g = 1:numG
        receptSep(r,g,:) = sqrt((rXCoord(r,:)-gXCoord(g,:)).^2 + (rYCoord(r,:)-gYCoord(g,:)).^2);
    end
end


winSizeMax = 25;
winSizeMin = 10;
lagTMax = 10;
dCorr = zeros(numR,numG,winSizeMax,lagTMax,numFrames);
pCorr = zeros(size(dCorr));

subplot(winSizeMax-winSizeMin+1,lagTMax,1);
for r = plotR
    for g = plotG
        for win = winSizeMin:winSizeMax
            winSize = win-winSizeMin+1;
            disp(winSize)
            for lagT = 1:lagTMax
                buffer = floor((win+lagT)/2);
                for t = 1:numFrames
                    if(receptSep(r,g,t)<.2)
                        %dCorr(r,g,winSize,lagT,t) = distanceCorrelation(receptTrajR(r,:,:),receptTrajG(g,:,:),t,win);
                        pCorr(r,g,winSize,lagT,t) = trajCorrelation(receptTrajR(r,:,:),receptTrajG(g,:,:),t,win,lagT);
                    else
                        dCorr(r,g,winSize,lagT,t) = 0;
                        pCorr(r,g,winSize,lagT,t) = 0;
                    end
                end
                subplot(winSizeMax-winSizeMin+1,lagTMax,(winSize-1)*lagTMax+lagT);
                plot(1:numFrames,squeeze(pCorr(r,g,winSize,lagT,:)),[475 475 525 525], [0 -1 -1 0]);
                xlim([1 numFrames]);
                ylim([-1.1 1.1]);
            end
        end
    end
end

receptorInfo = squeeze(pCorr(plotR,plotG,:,:,:));