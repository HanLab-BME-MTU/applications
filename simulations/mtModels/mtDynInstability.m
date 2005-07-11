function [traj,jump,errFlag] = mtDynInstability(modelParam,initialState,...
                totalTime,dt,timeEps,saveTraj)
%MTDYNINSTABILITY sim. the dyn. instabil. of a microtubule, free or in G1 yeast
%
%It uses a simple model where a microtubule is treated as a sequence of "units",
%which can be thought of as rings of 13 tubulin dimers. These units attach to and 
%detach from its plus end, while the minus end is fixed. The length of the
%microtubule must be always greater than a minimum length, below which it
%either (1) comletely disassembles if it's free, or (2) is rescued "manually" 
%if it's in G1 yeast. In the free case, a new microtubule is created and the 
%simulation proceeds.
%
%SYNOPSIS [traj,jump,errFlag] = mtDynInstability(modelParam,initialState,...
%               totalTime,dt,timeEps,saveTraj)
%
%INPUT  modelParam  : Structure containing model parameters:
%           .minLength   : minimum length of microtubule, in micrometers.
%           .kOnELong    : Rate constant of "unit" addition to elongating
%                          microtubule, in (micromolar*s)^-1.
%           .kOffELong   : Rate constant of "unit" loss in an elongating
%                          microtubule, in s^-1.
%           .kOnShrink   : Rate constant of "unit" addition to shrinking
%                          microtubule, in (micromolar*s)^-1.
%           .kOffShrink  : Rate constant of "unit" loss in a shrinking
%                          microtubule, in s^-1.
%
%       initialState: Structure containing information on initial state of system:
%           .mtLength0   : Initial length of microtubule, in micrometers.
%           .mtState0    : State of microtubule: -1 if shrinking or in a pause after
%                          shrinking, +1 if growing or in pause after growing.
%           .free        : 1 if microtubule is free, 0 if in G1 yeast.
%           .unitConc    : Concentration of free GTP-tubulin "units", in micromolar,
%                          assumed constant throughout the simulation.
%
%       totalTime   : Total time of simulation in seconds.
%       dt          : Time step used for time discretization. 
%       timeEps     : Value of the product of time step and maximum rate
%                     constant.
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
%OUTPUT traj        : 2 column vector, where first column is time (s) and
%                     second column is MT length (microns)
%       jump        : instances throughout the simulation where MT dies or is rescued.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, 9/2003

minLength = modelParam.minLength;
kOnElong = modelParam.kOnElong;
kOffElong = modelParam.kOffElong;
kOnShrink = modelParam.kOnShrink;
kOffShrink = modelParam.kOffShrink;

mtLength0 = initialState.mtLength0;
mtState0 = initialState.mtState0;
free = initialState.free;
unitConc = initialState.unitConc;

%effective rate constants of unit addition while growing and shrinking
kOnElongEff = kOnElong*unitConc;
kOnShrinkEff = kOnShrink*unitConc;

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('mtDynInstability')
    disp('--mtDynInstability: Incorrect number of input arguments!');
    errFlag  = 1;
    return;
end

%exception handling
% if minLength <= 0
%     disp('--mtDynInstability: Minimum length should be positive!');
%     errFlag = 1;
%     return;
% end
if mtLength0 <= minLength 
    disp('--mtDynInstability: Initial length of microtubule should be');
    disp(sprintf('  larger than %f micrometers!',minLength));
    errFlag = 1;
end
if mtState0 ~= -1 && mtState0 ~= +1
    disp('--mtDynInstability: mtState0 should be either +1 or -1!');
    errFlag = 1;
end 
if free ~= 0 && free~=1
    disp('--mtDynInstability: The variable "free" should be either 0 or 1!');
    errFlag = 1;
end    
if unitConc <= 0
    disp('--mtDynInstability: Tubulin concentration should be positive!');
    errFlag = 1;
end
if kOnElong < 0
    disp('--mtDynInstability: kOnElong should be positive!');
    errFlag = 1;
end
if kOffElong < 0
    disp('--mtDynInstability: kOffElong should be positive!');
    errFlag = 1;
end
if kOnShrink < 0
    disp('--mtDynInstability: kOnShrink should be positive!');
    errFlag = 1;
end
if kOffShrink < 0
    disp('--mtDynInstability: kOffShrink should be positive!');
    errFlag = 1;
end
if totalTime <= 0
    disp('--mtDynInstability: Total time should be positive!');
    errFlag = 1;
end
if timeEps > 1
    disp('--mtDynInstability: The product of the time step and maximum rate constant should be');
    disp('  smaller than 1');
    errFlag = 1;
end
if max([kOnElongEff kOffElong kOnShrinkEff kOffShrink])*dt > 1.001*timeEps
    disp('--mtDynInstability: Time step too large!');
    errFlag = 1;
end
if saveTraj.saveOrNot ~= 0 && saveTraj.saveOrNot ~= 1
    disp('--analyzeMtTrajectory: "saveTraj.saveOrNot" should be 0 or 1!');
    errFlag = 1;
end
if errFlag
    disp('--mtDynInstability: Please fix input data!');
    return;
end

%initialize random number generator
rand('state',sum(100*clock));

time = 0;                          %time of simulation, in seconds
iter = 1;                          %number of iterations
counter = 1;                       %number of times a MT dies or is rescued

vecLength = floor(totalTime/dt);   %approximate length of traj
traj = zeros(vecLength,2);         %allocate memory for array traj 
mtStateTemp = zeros(vecLength,1); 

unitLength = 0.008;                %length, in micrometers, of a tubulin "unit"
mtIndex0 = round(mtLength0/unitLength); %initial number of "units" in microtubule

traj(1,:) = [0 mtIndex0*unitLength]; %initial time and MT length

mtIndex = mtIndex0;                %current # of units = initial # of units 
mtState = mtState0;                %current state = initial state

jump(1) = 0;                       %i.e. MT was shorter than min. before sim. started

while time < totalTime             %iterate until totalTime is reached
    
    mtStateTemp(iter) = mtState;
    
    iter = iter + 1;               %update iteration numer
    time = time + dt;              %update time
    
    %determine whether a "unit" will be added or removed from the
    %microtubule, or whether no action will take place. Thus get the new
    %length and state of the system.
    switch mtState
        case +1                        %use elongation parameters
            pOn = kOnElongEff*dt;
            pOff = kOffElong*dt;
        case -1                        %use shrinkage parameters
            pOn = kOnShrinkEff*dt;
            pOff = kOffShrink*dt;
        otherwise
            disp('--mtDynInstability: mtState is neither +1 nor -1!');
            errFlag = 1;
            return;
    end
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pOn
    smallerOn = randomNumber < pOn;
    
    randomNumber = rand;   %get a random number between 0 and 1 and compare it to pOff
    smallerOff = randomNumber < pOff;
    
    increment = smallerOn - smallerOff;  %determine whether a unit is added or
                                         %removed, or whether there's a pause,
    if increment == 1  %determine new state
        mtState = +1;
    elseif increment == -1
        mtState = -1;
    end
    
    mtIndex = mtIndex + increment;           %get new number of "units" in MT
    traj(iter,:) = [time mtIndex*unitLength];%store new length of microtubule
    
    if traj(iter,2) <= minLength            %if microtubule reaches min length
        
        counter = counter + 1;              %increment counter by 1
        jump(counter) = iter;               %iter. # where min. length was reached
        
        iter = iter + 1;                    %new step
        time = time + dt; 
        
        if free
            %start new microtubule
            mtIndex = mtIndex0;
            traj(iter,2) = mtIndex*unitLength;
            mtState = mtState0;
        else
            %force a rescue
            mtIndex = mtIndex + 1;
            traj(iter,2) = mtIndex*unitLength;
            mtState = +1;
        end
        
    end
    
end

mtStateTemp(iter) = mtState;

counter = counter + 1;                      %end of simulation
jump(counter) = size(traj,1);               %total number of iterations in simulation
                                            %store even if MT did not reach
                                            %min. length in the end.
%save data if user wants to
if saveTraj.saveOrNot
    if isempty(saveTraj.fileName)
        save(['mtTraj-',nowString],'traj'); %save in file
    else
        save(saveTraj.fileName,'traj'); %save in file (directory specified through name)
    end
end
