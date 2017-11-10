function states = hiddenMarkovModel(rTraj,gTraj,locError,pOn,pOff,dT,diffConst,boxSize)

%HIDDENMARKOVMODEL returns the path of states for a given set of trajectories and parameters
%
%SYNPOSIS x = hiddenMarkovModel(rTraj,gTraj,locError,pOn,pOff,dT,diffConst,boxSize)
%
%INPUT  rTraj              : The 2D position of the red particle over time
%       gTraj              : The 2D position of the green particle over
%                            time
%       locError           : The one-dimensional localization error for all
%                            the particles. More work is needed to change
%                            this for each time point and particle from the
%                            observed data.
%       pOn                : The estimated probability of binding for each
%                            time step
%       pOff               : The estimated probability of splitting for
%                            each time step
%       dT                 : The time between each subsequent frame
%       diffConst          : Diffusion constant (um^2/s)
%       boxSize            : The length of each side of the square (in um)
%                            for the given data set
%
%OUTPUT states             : The path of states predicted by a hidden
%                            Markov model (1 = monomer, 2 = dimer)
%
%Code by Paul Blazek June 2015.

%% Constants
%variances for the two different probability distributions used for each
%state
var1 = 2*locError^2;
var2 = 8*locError^2 + 2*diffConst*dT;

%% Fix trajectories
%squeeze them and make sure they are in the proper orientation.

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
%The Viterbi algorithm is used to find which path of states is most likely
%given the markov model and the observed state.

%Note for below: I switched the L matrix for a logL matrix because my
%values were going below realmin

separations = sqrt(sum((rTraj-gTraj).^2,2));

%likelihoods at first time point
startP = ones(2,1);
startP(2) = separations(1)/var1*exp(-separations(1)^2/(2*var1));
startP(1) = 2*separations(1)/boxSize^2;

%set up likelihood matrix and path
logL = zeros(2,numFrames);
possiblePaths = zeros(2,numFrames);
P = zeros(2,1);

if startP(1) == 0
    logL(1,1) = log(realmin);
else
    logL(1,1) = log(startP(1));
end
logL(2,1) = log(startP(2));

%T is the transition matrix, or the probability of going from one state to
%another
T = zeros(2);
T(2,1) = pOff;
T(2,2) = 1 - pOff;
for t=2:numFrames %loop across the whole time
    %Calculate the likelihoods of being in state 1 or 2; for state 1 there
    %is a cutoff after which the likelihood is 0; the cutoff is the
    %likelihood of being beyond the end of the matrix used to store the
    %probabilities of being within the cutoff distance
    P(1) = 0;
    for theta = 0:pi/100:2*pi*.9999
        P(1) = P(1) + exp(-(separations(t)^2 + separations(t-1)^2 - 2*separations(t)*separations(t-1)*cos(theta))/2/var2)*pi/100;
    end
    P(1) = P(1)*separations(t)/2/pi/var2;
    P(2) = separations(t)/var1*exp(-(separations(t)^2)/(2*var1));
    P(2) = P(2)*(P(2)>1e-10);
    
    %The below was blotted out because it was found to be less accurate;
    %basically it was scaling pOn by the probability that the separation
    %difference between the two particles was less than some cutOff
    %parameter
    
    %create the rest of the T matrix, which is dependent on the separation
    %if (separations(t) < lessThanCutoffProb(end,1))
    %    [~,sepApprox] = min(abs(separations(t)-lessThanCutoffProb(:,1)));
    %    T(1,2) = pOn*lessThanCutoffProb(sepApprox,2);
    %else
    %    T(1,2) = 0;
    %end
    T(1,2) = pOn;
    T(1,1) = 1 - T(2,1);
    
    %find which end point is more probable for each starting condition, and  
    if P(2) == 0
        [logL(2,t), possiblePaths(2,t)] = max(logL(:,t-1) + 15*log(realmin));
    else
        [logL(2,t), possiblePaths(2,t)] = max(logL(:,t-1) + log(T(:,2).*P(2)));
    end
    [logL(1,t), possiblePaths(1,t)] = max(logL(:,t-1) + log(T(:,1).*P(1)));
    
end

%states is the true path (i.e. the answer of the status)
states = zeros(1,numFrames);

%choose which of the start conditions based on overall likelihood
[~,states(end)] = max(logL(:,numFrames));

%now that you have that, find the actual path and store it in x
for t = numFrames:-1:2
    states(t-1) = possiblePaths(states(t),t);
end

%median filter to remove lifetimes of only 1 frame
for i = 2:(numFrames-1)
    if and(states(i)~=states(i-1),states(i)~=states(i+1))
        states(i) = states(i-1);
    end
end