function [mtLength,capSize,errFlag] = mtGTPCapLDepK(modelParam,initialState,...
    totalTime,dt,timeEps,saveTraj)
%MTGTPCAPLDEPK: GTP-cap model with length-dependent k's of microtubule DI
%
%It uses a simple model where a microtubule is treated as a sequence of "units",
%which can be thought of as rings of 13 tubulin dimers. These units attach to and 
%detach from its plus end, while the minus end is fixed. "units" attached to the 
%microtubule may get hydrolyzed, if the "unit" preceding them is. Generally, the 
%microtubule grows when there is a GTP cap and shrinks when the GDP cap is
%lost. In order to limit MT lengths to a certain range, the rate of GTP-"unit" 
%addition to GTP-end is reduced as MT grows and the rate of GDP-"unit" loss
%is reduced as MT shrinks. 
%
%SYNOPSIS [mtLength,capSize,errFlag] = mtGTPCapLDepK(modelParam,initialState,...
%    totalTime,dt,timeEps,saveToFile)
%
%INPUT  modelParam  : Structure containing model parameters:
%           .minLength   : minimum length of microtubules, in micrometers.
%           .maxLength   : maximum length of microtubules, in micrometers.
%           .kTOnT = kTOnTFree - addAmpT*{tanh[addWidT*(mtLength-addLenT)]+1}/2.
%                        : Rate constant of GTP-"unit" addition to GTP-"unit"
%                          at + end, in (micromolar*s)^-1. Units:
%                            [kTOnTFree],[addAmpT] = (micromolar*s)^-1.
%                            [addWidT] = (micrometers)^-1.
%                            [addLenT] = micrometers.
%           .kTOff       : Rate constant of GTP-"unit" loss from + end, in s^-1.
%           .kTOnD       : Rate constant of GTP-"unit" addition to GDP-"unit"
%                          at + end, in (micromolar*s)^-1.
%           .kDOff = kDOffFree + addAmpD*{tanh[addWidD*(mtLength-addLenD)]-1}/2.
%                        : Rate constant of GDP-"unit" loss from + end, in s^-1.
%                            [kDOffFree],[addAmpD] = s^-1.
%                            [addWidD] = (micrometers)^-1.
%                            [addLenD] = micrometers.
%           .kHydrolysis : Rate constant of hydrolysis of GTP-"units" to
%                          GDP-"units" in microtubule.
%
%       initialState: Structure containing information on initial state of system:
%           .mtLength0   : Initial length of microtubule, in micrometers.
%           .capSize0    : Initial number of "units" forming GTP-cap.
%           .unitConc    : Concentration of free GTP-tubulin "units", in micromolar,
%                          assumed constant throughout the simulation.
%
%       totalTime   : Total time of simulation in seconds.
%       dt          : Time step used for time discretization. 
%       timeEps     : Value of the product of time step and maximum rate
%                     constant.
%
%       saveTraj    : Structure defining whether and where results will be
%                     saved.
%           .saveOrNot   : 1 if user wants to save, 0 if not.
%           .fileName    : name (including location) of file where results 
%                          will be saved. If empty and saveOrNot is 1, the name
%                          is chosen automatically to be
%                          "mtTraj-day-month-year-hour-minute-second",
%                          and the data is saved in directory where
%                          function is called from
%       
%OUTPUT mtLength    : Length of the microtubule throughout the simulation.
%       capSize     : Number of units forming GTP-cap throughout the simulation.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, 10/2003

errFlag = 0;                    %logical variable used to indicate whether 
                                %there is an error in input data

%check if correct number of arguments were used when function was called
if nargin ~= nargin('mtGTPCapLDepK');
    disp('--mtGTPCapLDepK: Incorrect number of input arguments!');
    errFlag = 1;
    return;
end

mtLength0 = initialState.mtLength0;
capSize0  = initialState.capSize0;
unitConc  = initialState.unitConc;
minLength = modelParam.minLength;
maxLength = modelParam.maxLength;
kTOnTFree = modelParam.kTOnTFree;
addAmpT = modelParam.addAmpT;
addWidT = modelParam.addWidT;
addLenT = modelParam.addLenT;
kTOff = modelParam.kTOff;
kTOnD = modelParam.kTOnD;
kDOffFree = modelParam.kDOffFree;
addAmpD = modelParam.addAmpD;
addWidD = modelParam.addWidD;
addLenD = modelParam.addLenD;
kHydrolysis = modelParam.kHydrolysis;

%check for error in input
if minLength < 0
    disp('--mtGTPCapLDepK: Minimum length should be positive!');
    errFlag = 1;
    return;
end
if maxLength <= minLength
    disp('--mtGTPCapLDepK: Maximum length should be positive and larger than minimum length!');
    errFlag = 1;
    return;
end
if addWidT < 6/(maxLength-minLength)
    disp('--mtGTPCapLDepK: addWidT is too small!');
    errFlag = 1;
    return;
end
if addWidD < 6/(maxLength-minLength)
    disp('--mtGTPCapLDepK: addWidD is too small!');
    errFlag = 1;
    return;
end
if mtLength0 <= minLength || mtLength0 >= maxLength 
    disp('--mtGTPCapLDepK: Initial length of microtubule should be');
    disp(sprintf('  between %f and %f micrometers!',minLength,maxLength));
    errFlag = 1;
end
if capSize0 < 0
    disp('--mtGTPCapLDepK: capSize0 cannot be negative!');
    errFlag = 1;
end
if unitConc <= 0
    disp('--mtGTPCapLDepK: Tubulin concentration should be positive!');
    errFlag = 1;
end
if kTOnTFree < 0
    disp('--mtGTPCapLDepK: kTOnTFree cannot be negative!');
    errFlag = 1;
end
if addAmpT < 0 || addAmpT >= kTOnTFree
    disp('--mtGTPCapLDepK: addAmpT should be positive and smaller than kTOnTFree!');
    errFlag = 1;
end
if addLenT < minLength + 3/addWidT || addLenT > maxLength - 3/addWidT
    disp('--mtGTPCapLDepK: addLenT is not within acceptable range!');
    errFlag = 1;
end
if kTOff < 0
    disp('--mtGTPCapLDepK: kTOff cannot be negative!');
    errFlag = 1;
end
if kTOnD < 0
    disp('--mtGTPCapLDepK: kTOnD cannot be negative!');
    errFlag = 1;
end
if kDOffFree < 0
    disp('--mtGTPCapLDepK: kDOff cannot be negative!');
    errFlag = 1;
end
if addAmpD < 0 || addAmpD >= kDOffFree
    disp('--mtGTPCapLDepK: addAmpD should be positive and smaller than kDOffFree!');
    errFlag = 1;
end
if kHydrolysis < 0
    disp('--mtGTPCapLDepK: kHydrolysis cannot be negative!');
    errFlag = 1;
end 
if totalTime <= 0
    disp('--mtGTPCapLDepK: Total time should be positive!');
    errFlag = 1;
end
if timeEps > 1
    disp('--mtGTPCapLDepK: The product of the time step and maximum rate');
    disp('  constant should be smaller than 1!');
    errFlag = 1;
end
kTOnTFreeEff = kTOnTFree*unitConc;   %effective rate constants of "unit" addition
addAmpTEff = addAmpT*unitConc;       %(s^-1)
kTOnDEff = kTOnD*unitConc;
if max([kTOnTFreeEff kTOff kTOnDEff kDOffFree kHydrolysis])*dt > timeEps
    disp('--mtGTPCapLDepK: Time step too large!');
    errFlag = 1;
end
if saveTraj.saveOrNot ~= 0 && saveTraj.saveOrNot ~= 1
    disp('--analyzeMtTrajectory: "saveTraj.saveOrNot" should be 0 or 1!');
    errFlag = 1;
end
if errFlag
    disp('--mtGTPCapLDepK: Please fix input data!');
    mtLength = [];
    capSize = [];
    return;
end

%initialize random number generator
rand('state',sum(100*clock));

unitLength = 0.008;   %length, in micrometers, of a tubulin "unit"

mtIndex = round(mtLength0/unitLength); %number of "units" in microtubule
currentLength = mtIndex*unitLength;    %length of microtubule
if capSize0 > mtIndex                  %size of cap cannot be larger than number 
    capSize0 = mtIndex;                %of units in MT
end
capSizeT = capSize0;                   %initial number of units in GTP-cap

%first run a 2000 second simulation to get rid of any artifacts in the
%initial data, namely the size of the GTP-cap.

time = 0;
while time < 2000        %iterate for required "equilibration time"
    
    if capSizeT == 0     %mtState = +1 if there is GTP cap
        mtState = -1;    %and -1 if there isn't
    else
        mtState = +1;
    end
    
    time = time + dt;    %update time
    
    %determine whether a "unit" will be added or removed from the
    %microtubule + end, or whether no action will take place. 
    %also determine whether the "unit" at the edge of the GTP-cap, if it
    %exists, will get hydrolyzed.
    switch mtState
        case +1                           %use parameters when there is GTP-cap
            kTOnTEff = kTOnTFreeEff - addAmpTEff*...
                (tanh(addWidT*(currentLength-addLenT))+1)/2;
            pOn = kTOnTEff*dt;
            pOff = kTOff*dt;
            pHyd = kHydrolysis*dt;
        case -1                           %use parameters when there is no GTP-cap
            kDOff = kDOffFree + addAmpD*...
                (tanh(addWidD*(currentLength-addLenD))-1)/2;
            pOn = kTOnDEff*dt;
            pOff = kDOff*dt;
            pHyd = NaN;
        otherwise
            disp('--mtGTPCapLDepK: mtState is neither +1 nor -1!');
            errFlag = 1;
            return;
    end
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pOn
    smallerOn = randomNumber < pOn;
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pOff
    smallerOff = randomNumber < pOff;
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pHyd
    smallerHyd = randomNumber < pHyd;
    
    increment = smallerOn - smallerOff;   %determine whether a unit is added or
                                          %removed, or whether there's a pause,
    
    mtIndex = mtIndex + increment;                      %update number of units in MT
    capSizeT = max(0,capSizeT-smallerHyd+increment);    %and in cap

    currentLength = mtIndex*unitLength;       %store new length of microtubule
        
end

%Now do the actual simulation

vecLength = floor(totalTime/dt);   %approximate length of mtLength and capSize 
mtLength = zeros(vecLength,1);     %allocate memory for arrays mtLength
capSize = zeros(vecLength,1);      %and capSize

mtLength(1) = currentLength;       %initial length of microtubule
capSize(1) = capSizeT;             %initial number of units in GTP-cap

time = 0;                          %time of simulation, in seconds
iter = 1;                          %number of iterations

while time < totalTime             %iterate until totalTime is reached
    
    if capSize(iter) == 0          %mtState = +1 if there is GTP cap
        mtState = -1;              %and -1 if there isn't
    else
        mtState = +1;
    end
    currentLength = mtLength(iter);
    
    iter = iter + 1;               %update iteration numer
    time = time + dt;              %update time
    
    %determine whether a "unit" will be added or removed from the
    %microtubule + end, or whether no action will take place. 
    %also determine whether the "unit" at the edge of the GTP-cap, if it
    %exists, will get hydrolyzed.
    switch mtState
        case +1                           %use parameters when there is GTP-cap
            kTOnTEff = kTOnTFreeEff - addAmpTEff*...
                (tanh(addWidT*(currentLength-addLenT))+1)/2;
            pOn = kTOnTEff*dt;
            pOff = kTOff*dt;
            pHyd = kHydrolysis*dt;
        case -1                           %use parameters when there is no GTP-cap
            kDOff = kDOffFree + addAmpD*...
                (tanh(addWidD*(currentLength-addLenD))-1)/2;
            pOn = kTOnDEff*dt;
            pOff = kDOff*dt;
            pHyd = NaN;
        otherwise
            disp('--mtGTPCapLDepK: mtState is neither +1 nor -1!');
            errFlag = 1;
            return;
    end
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pOn
    smallerOn = randomNumber < pOn;
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pOff
    smallerOff = randomNumber < pOff;
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pHyd
    smallerHyd = randomNumber < pHyd;
    
    increment = smallerOn - smallerOff;   %determine whether a unit is added or
                                          %removed, or whether there's a pause,
    
    mtIndex = mtIndex + increment;                %update number of units in MT
    capSize(iter) = max(0,capSize(iter-1)-smallerHyd+increment);    %and in cap

    mtLength(iter) = mtIndex*unitLength;       %store new length of microtubule
    
end

%save data if user wants to
if saveTraj.saveOrNot
    if isempty(saveTraj.fileName)
        save(['mtTraj-',nowString],'mtLength','capSize'); %save in file
    else
        save(saveTraj.fileName,'mtLength','capSize'); %save in file (directory specified through name)
    end
end
