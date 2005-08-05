function [simTraj,errFlag] = mtGammaTdSd(modelParam,length0,totalTime)
%MTGAMMATDSD simulates an MT trajectory assuming that phase duration is gamma-distributed and speed is normally-distributed
%
%SYNOPSIS [simTraj,errFlag] = mtGammaTdSd(modelParam,length0,totalTime)
%
%INPUT  modelParam  : Structure containing model parameters:
%           .growthSpeed   : 2-element row vector containing mean and
%                            standard deviation of growth rate [microns/min].
%           .shrinkageSpeed: 2-element row vector containing mean and
%                            standard deviation of shrinkage rate [microns/min].
%           .growthTime    : 2-element row vector containing mean and
%                            standard deviation of growth time [s].
%           .shrinkageTime : 2-element row vector containing mean and
%                            standard deviation of shrinkage time [s].
%
%       length0     : Initial length of microtubule [microns].
%       totalTime   : Total time of simulation [s].
%       
%OUTPUT simTraj     : 2 column vector. 1st column: Time at transitions
%                     (rescue and catastrophe). 2nd column: Corresponding
%                     MT length.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simTraj = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin ~= 3 

    disp('--mtGammaTimeDistr: Wrong number of input variables!');
    errFlag = 1;
    return

end

%get model parameters
vG = modelParam.growthSpeed/60;
vS = modelParam.shrinkageSpeed/60;
tG = modelParam.growthTime;
tS = modelParam.shrinkageTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate approximate number of iterations needed
totItno = round(totalTime/mean([tG tS]));

%get totItno uniformly distributed random numbers
rand('state',sum(100*clock));
randNum = rand(totItno,1);

%calculate phase duration gamma distribution parameters
rG = (tG(1)/tG(2))^2;
invThetaG = tG(2)^2/tG(1);
rS = (tS(1)/tS(2))^2;
invThetaS = tS(2)^2/tS(1);

%calculate growth and shrinkage times
gsTimes = zeros(totItno,1);
gsTimes(1:2:end) = gaminv(randNum(1:2:end),rG,invThetaG);
gsTimes(2:2:end) = gaminv(randNum(2:2:end),rS,invThetaS);

%get totItno uniformly distributed random numbers
randNum = rand(totItno,1);

%calculate speed gamma distribution parameters
rG = (vG(1)/vG(2))^2;
invThetaG = vG(2)^2/vG(1);
rS = (vS(1)/vS(2))^2;
invThetaS = vS(2)^2/vS(1);

%determine growth and shrinkage speeds in each growth and shrinkage
%intervel, respectively
gsSpeeds = zeros(totItno,1);
gsSpeeds(1:2:end) = gaminv(randNum(1:2:end),rG,invThetaG);
gsSpeeds(2:2:end) = -gaminv(randNum(2:2:end),rS,invThetaS);

%calculate length change in each interval
lengthChange = zeros(totItno,1);
lengthChange(1:2:end) = gsTimes(1:2:end).*gsSpeeds(1:2:end);
lengthChange(2:2:end) = gsTimes(2:2:end).*gsSpeeds(2:2:end);

%calculate simTraj
simTraj = zeros(totItno,2);
simTraj(1,:) = [0 length0]; %initial length
for i=1:totItno
    simTraj(i+1,:) = simTraj(i,:) + [gsTimes(i) lengthChange(i)];
end


%%%%% ~~ the end ~~ %%%%%

