function x = hiddenMarkovModel(rTraj,gTraj,locError)

%% Constants
kOff = 1;
kOn = .5;
dT = 0.01;
diffConst = 0.1;
boxSize = 10;
var1 = locError^2;
var2 = 2*locError^2 + 2*diffConst*dT;
pOff = kOff*dT;
pOn = kOn*dT;

%% Fix trajectories
rTraj = squeeze(rTraj);
gTraj = squeeze(gTraj);

if size(rTraj) == size(gTraj')
    gTraj = gTraj';
end

[numFrames, dim] = size(rTraj);
if numFrames < dim
    rTraj = rTraj';
    gTraj = gTraj';
    [numFrames, ~] = size(rTraj);
end

%% Viterbi algorithm with hidden Markov Model
obs = sqrt(sum((rTraj-gTraj).^2,2));
startP = [kOn/kOff (1-kOn/kOff)]';
transP = [(1-pOn) pOn; pOff (1-pOff)]';

V = zeros(2,numFrames);
path = V;
emitP = zeros(2,1);

V(1,1) = startP(1)*(2*obs(1)/boxSize^2);
V(2,1) = startP(2)*obs(1)/var1*exp(-obs(1)^2/(2*var1));

for t=2:numFrames
    emitP(1) = 0;
    for theta = 0:pi/100:2*pi
        emitP(1) = emitP(1) + exp(-(obs(t)^2 + obs(t-1)^2 - 2*obs(t)*obs(t-1)*cos(theta))/2/var2)*pi/100;
    end
    emitP(1) = emitP(1)*obs(t)/2/pi/var2;
    emitP(2) = obs(t)/var1*exp(-obs(t)^2/(2*var1));
    
    [V(1,t), path(1,t)] = max(V(:,t-1).*transP(:,1).*emitP(1));
    [V(2,t), path(2,t)] = max(V(:,t-1).*transP(:,2).*emitP(2));
end

x = zeros(1,numFrames);

[~,x(end)] = max(V(:,numFrames));

for t = numFrames:-1:2
    x(t-1) = path(x(t),t);
end

x = x - 1;