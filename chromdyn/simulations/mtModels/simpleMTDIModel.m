function [simTraj,errFlag] = simpleMTDIModel(modelParam,length0,totalTime)
%SIMPLEMTDIMODEL simulates an MT trajectory using the 4 traditional MTDI descriptors
%
%Simulation assumes a gamma distribution for growth and shrinkage times, as
%in D. J. Odde and H. M. Buettner 1995, Annals of Biomedical Engineering
%23:268-286.
%
%SYNOPSIS [simTraj,errFlag] = simpleMTDIModel(modelParam,length0,totalTime)
%
%INPUT  modelParam  : Structure containing model parameters:
%           .growthSpeed   : MT growth rate [microns/min].
%           .shrinkageSpeed: MT shrinkage rate [microns/min].
%           .growthTime    : 2-element row vector containing mean and
%                            standard deviation of growth time 
%                            distribution [s].
%           .shrinkageTime : 2-element row vector containing mean and
%                            standard deviation of shrinkage time 
%                            distribution [s].
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

    disp('--simpleMTDIModel: Wrong number of input variables!');
    errFlag = 1;
    return

end

%get model parameters
vG = modelParam.growthSpeed/60;
vS = -modelParam.shrinkageSpeed/60;
tG = modelParam.growthTime;
tS = modelParam.shrinkageTime;

%calculate gamma distribution parameters
rG = (tG(1)/tG(2))^2;
invThetaG = tG(2)^2/tG(1);
rS = (tS(1)/tS(2))^2;
invThetaS = tS(2)^2/tS(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate approximate number of iterations needed
totItno = round(totalTime/mean([tG tS]));

%get totItno uniformly distributed random numbers
rand('state',sum(100*clock));
randNum = rand(totItno,1);

%calculate growth and shrinkage times (which are gamma distributed)
gsTimes = zeros(totItno,1);
gsTimes(1:2:end) = gaminv(randNum(1:2:end),rG,invThetaG);
gsTimes(2:2:end) = gaminv(randNum(2:2:end),rS,invThetaS);

%calculate length change in each interval
lengthChange = zeros(totItno,1);
lengthChange(1:2:end) = gsTimes(1:2:end)*vG;
lengthChange(2:2:end) = gsTimes(2:2:end)*vS;

%calculate simTraj
simTraj = zeros(totItno,2);
simTraj(1,:) = [0 length0]; %initial length
for i=1:totItno
    simTraj(i+1,:) = simTraj(i,:) + [gsTimes(i) lengthChange(i)];
end


%%%%% ~~ the end ~~ %%%%%

