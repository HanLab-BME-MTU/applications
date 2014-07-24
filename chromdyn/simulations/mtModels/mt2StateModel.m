function [traj,errFlag] = mtKmcModel(modelParam,initialState,totalTime)
%MTKMCMODEL simulates MTDI using the KMC method
%
%INPUT  modelParam  : Structure containing model parameters:
%           .kGtpAdd  : Rate constant of GTP-unit addition to growing 
%                       plus-end [(micromolar*s)^-1]. Row vector with 3
%                       entries: fast rate, slow rate and probability of
%                       fast rate.
%           .kCapLoss : Rate constant of cap loss (via hydrolysis from MT 
%                       plus-end [s^-1]. Row vector with 2 entries: fast
%                       rate and slow rate.
%           .kGdpRem  : Rate constant of GDP-unit loss from shrinking 
%                       minus-end [s^-1]. Row vector with 3 entries: fast 
%                       rate, slow rate and probability of fast rate.
%           .kCapGain : Rate constant of cap gain GTP-unit addition to GDP-end
%                       [(micromolar*s)^-1]. Row vector with 2 entries: fast
%                       rate and slow rate.
%       initialState: Structure containing information on initial state of system:
%           .mtLength0: Initial length of microtubule [microns].
%           .unitConc : Concentration of free GTP-tubulin units [micromolar].
%                       Assumed constant throughout the simulation.
%       totalTime   : Total simulation time [s].
%
%OUTPUT traj        : 2 column vector of time [s] and MT length [microns]
%                     at each time point.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Model assumes coupled hydrolysis
%
%Khuloud Jaqaman, 6/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traj = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin ~= 3
    disp('--mtKmcModel: Incorrect number of input arguments!');
    errFlag = 1;
    return;
end

%get simulation parameters
mtLength0 = initialState.mtLength0;
unitConc  = initialState.unitConc;
kGtpAdd = modelParam.kGtpAdd(1:2)*unitConc;
pGtpAdd = modelParam.kGtpAdd(3);
kCapLoss = modelParam.kCapLoss;
kGdpRem = modelParam.kGdpRem(1:2);
pGdpRem = modelParam.kGdpRem(3);
kCapGain = modelParam.kCapGain*unitConc;

%assign unit length [microns]
unitLength = 0.008;

%get initial number of units in MT and assign mtLength0 making sure that it
%has an integer number of units
mtIndex = round(mtLength0/unitLength);
mtLength0 = unitLength*mtIndex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%length of mtLength
vecLength = 1.5*round(totalTime*max(kGtpAdd(1)+kCapLoss(1),kGdpRem(1)+kCapGain(1)));

%allocate memory for mtLength
traj = NaN*ones(vecLength,2);

%initialize variables
iter = 1;
time = 0;
traj(iter,:) = [time mtLength0];
mtLengthT = mtLength0;
capSizeT = 1;
mType = 2;

%choose random numbers for Monte Carlo simulation
rand('state',sum(100*clock)); %initialization
randNum = rand(vecLength,3);

%iterate until totalTime is reached
while time < totalTime

    %update iter by 1
    iter = iter + 1;

    if capSizeT == 0 %shrinkage state

        %calculate "total" rate
        kTot = kGdpRem(mType) + kCapGain(mType);
        
        %get time step
        dt = -log(randNum(iter,1))/kTot;

        %put event probabilities in a line from 0 to 1
        pEvent(1) = kGdpRem(mType)/kTot; %losing a GDP unit
        pEvent(2) = 1; %gaining a GTP cap

        %decide which event takes place
        event = find(pEvent>randNum(iter,2),1,'first');
        
        %take action
        switch event
            case 1 %lose a GDP unit
                mtLengthT = mtLengthT - unitLength;
            case 2 %gain cap
                mtLengthT = mtLengthT + unitLength;
                capSizeT = 1;
                mType = randNum(iter,3)>pGtpAdd;
                mType = mType + 1;
        end

    else %growth state

        %calculate "total" rate
        kTot = kGtpAdd(mType) + kCapLoss(mType);

        %get time step
        dt = -log(randNum(iter,1))/kTot;

        %put event probabilities in a line from 0 to 1
        pEvent(1) = kGtpAdd(mType)/kTot; %adding a GTP unit
        pEvent(2) = 1; %losing the GTP cap

        %decide which event takes place
        event = find(pEvent>randNum(iter,2),1,'first');

        %take action
        switch event
            case 1 %add GTP unit
                mtLengthT = mtLengthT + unitLength;
            case 2 %lose cap
                mtLengthT = mtLengthT - unitLength;
                capSizeT = 0;
                mType = randNum(iter,3)>pGdpRem;
                mType = mType + 1;
        end

    end

    %update time
    time = time + dt;

    %update MT length in vector
    traj(iter,:) = [time mtLengthT];
    
end

%remove extra time points in the end
traj = traj(find(~isnan(traj(:,1))),:);


%%%%% ~~ the end ~~ %%%%%
