function [traj,errFlag] = mtGTPCapLDepK(modelParam,initialState,...
    totalTime,dt,timeEps,saveTraj)
%MTGTPCAPLDEPK: GTP-cap model of MTDI with length-dependent k's
%
%It uses a simple model where a microtubule is treated as a sequence of units,
%which can be thought of as rings of 13 tubulin dimers. These units attach to and
%detach from its plus end, while the minus end is fixed. A Unit attached to the
%microtubule may get hydrolyzed if the unit preceding it is. Generally, the
%microtubule grows when there is a GTP cap and shrinks when the GTP cap is
%lost due to hydrolysis. All rates can be length-dependent.
%
%SYNOPSIS [traj,errFlag] = mtGTPCapLDepK(modelParam,initialState,...
%    totalTime,dt,timeEps,saveToFile)
%
%INPUT  modelParam  : Structure containing model parameters:
%           .kTOnT    : Rate constant of GTP-unit addition to GTP-end
%                       [(micromolar*s)^-1].
%           .kHydr    : Rate constant of hydrolysis of GTP-units in MT [s^-1].
%           .kTOff    : Rate constant of GTP-unit loss from end [s^-1].
%           .kTOnD    : Rate constant of GTP-unit addition to GDP-end
%                       [(micromolar*s)^-1].
%           .kDOff    : Rate constant of GDP-unit loss from end [s^-1].
%
%           All rate constants have the following subfields:
%               .lDepend : +1 if rate increases with MT length;
%                          -1 if rate decreases with MT length;
%                           0 if rate does not depend on length.
%                           If rate is not length dependent, only kmax
%                           must be specified.
%               .kmax    : Maximum value that rate constant can have.
%               .kmin    : Minimum value that rate constant can have.
%               .lmin    : MT length where changes "start" [microns].
%               .lmax    : MT length where changes "stop" [microns].
%           In addition, kHydr has the subfield
%               .coupled : 1 if hydrolysis is coupled to unit addition 
%                          (only the last added unit might have GTP); 
%                          0 if hydrolysis is not directly coupled to
%                          addition (more than one unit at the + end might 
%                          have GTP).
%
%       initialState: Structure containing information on initial state of system:
%           .mtLength0: Initial length of microtubule [microns].
%           .capSize0 : Initial number of "units" forming GTP-cap [units].
%           .unitConc : Concentration of free GTP-tubulin units [microns].
%                          Assumed constant throughout the simulation.
%
%       totalTime   : Total simulation time [s].
%       dt          : Time step used for time discretization [s].
%       timeEps     : Value of the product of time step and maximum rate
%                     constant. Optional. Default: 0.5. Use [] for default.
%
%       saveTraj    : Structure defining whether and where results are
%                     saved. Optional. Default: no saving.
%           .saveOrNot: 1 if user wants to save, 0 if not.
%           .fileName : name (including location) of file where results
%                       will be saved. If empty and saveOrNot is 1, the name
%                       is chosen automatically to be
%                       "mtTraj-day-month-year-hour-minute-second",
%                       and the data is saved in directory where function
%                       is called from.
%
%OUTPUT traj        : 3-column vector:
%                       1st column: time [s];
%                       2nd column: MT length [microns];
%                       3rd column: cap size [rings].
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, 10/2003. Major updates: 6/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traj = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < 4
    disp('--mtGTPCapLDepK: Incorrect number of input arguments!');
    errFlag = 1;
    return;
end

%get simulation parameters
mtLength0 = initialState.mtLength0;
capSize0  = initialState.capSize0;
unitConc  = initialState.unitConc;
kTOnT = modelParam.kTOnT;
kHydr = modelParam.kHydr;
kTOff = modelParam.kTOff;
kTOnD = modelParam.kTOnD;
kDOff = modelParam.kDOff;

%use default if timeEps is not given
if nargin < 5 | isempty(timeEps)
    timeEps = 0.5;
end

%use default if saveTraj is not given
if nargin < 6
    saveTraj.saveOrNot = 0;
end

%make sure that simulation time step is not too large
kMax = max([kTOnT.kmax*unitConc kHydr.kmax kTOff.kmax kTOnD.kmax*unitConc ...
    kDOff.kmax]);
if kMax*dt > timeEps
    disp('--mtGTPCapLDepK: Time step too large!');
    errFlag = 1;
end

%exit if there are problems in input
if errFlag
    disp('--mtGTPCapLDepK: Please fix input data!');
    return;
end

%assign unit length [microns]
unitLength = 0.008;

%get initial number of units in MT and assign mtLength0 making sure that it
%has an integer number of units
mtIndex = round(mtLength0/unitLength);
mtLength0 = unitLength*mtIndex;

%make sure that # units in cap <= # units in MT
if capSize0 > mtIndex
    capSize0 = mtIndex;
end

%calculate all parameters specifying length dependence of rate constants
kTOnT = deriveParam(kTOnT);
kHydr = deriveParam(kHydr);
kTOff = deriveParam(kTOff);
kTOnD = deriveParam(kTOnD);
kDOff = deriveParam(kDOff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulation is done for "totalTime + 100s", where in the end the first 100 s
%are removed to eliminate any artifacts in the initial data

%length of traj
vecLengthI = round(100/dt);
vecLength = vecLengthI + round(totalTime/dt);

%allocate memory for traj
traj = zeros(vecLength,3);

%assign initial time, MT length and cap size
time = 0;
mtLengthT = mtLength0;
capSizeT = capSize0;
traj(1,:) = [time mtLengthT capSizeT];

%choose random numbers for Monte Carlo simulation
rand('state',sum(100*clock)); %initialization
randNum = rand(vecLength,3);

%iterate until "totalTime + 100" is reached
for iter = 2:vecLength
    
    time = time + dt;
    
    if capSizeT == 0 %shrinkage state

        %get rates based on current MT length        
        kTOnDEff = kTOnD.kmin + (kTOnD.kmax-kTOnD.kmin)*...
            (kTOnD.lDepend*(tanh(kTOnD.coef*(mtLengthT-kTOnD.shift)))+1)/2;
        kDOffEff = kDOff.kmin + (kDOff.kmax-kDOff.kmin)*...
            (kDOff.lDepend*(tanh(kDOff.coef*(mtLengthT-kDOff.shift)))+1)/2;

        %assign probability of each event taking place
        pOn = kTOnDEff*unitConc*dt;
        pOff = kDOffEff*dt;
        pHyd = NaN;

    else %growth state

        %get rates based on current MT length
        kTOnTEff = kTOnT.kmin + (kTOnT.kmax-kTOnT.kmin)*...
            (kTOnT.lDepend*(tanh(kTOnT.coef*(mtLengthT-kTOnT.shift)))+1)/2;
        kHydrEff = kHydr.kmin + (kHydr.kmax-kHydr.kmin)*...
            (kHydr.lDepend*(tanh(kHydr.coef*(mtLengthT-kHydr.shift)))+1)/2;
        kTOffEff = kTOff.kmin + (kTOff.kmax-kTOff.kmin)*...
            (kTOff.lDepend*(tanh(kTOff.coef*(mtLengthT-kTOff.shift)))+1)/2;

        %assign probability of each event taking place
        pOn = kTOnTEff*unitConc*dt;
        pHyd = kHydrEff*dt;
        pOff = kTOffEff*dt;

    end

    %using these 3 random numbers, decide which events take place
    addUnit = randNum(iter,1) < pOn; %unit addition
    remUnit = randNum(iter,2) < pOff; %unit loss
    hydUnit = randNum(iter,3) < pHyd; %unit hydrolysis

    %find net unit addition or loss
    increment = addUnit - remUnit;

    %update MT length
    mtLengthT = mtLengthT + unitLength*increment;
    
    %update GTP cap size
    if ~kHydr.coupled %hydrolysis not coupled to addition
        capSizeT = max(0,capSizeT-hydUnit+increment);
    else %hydrolysis coupled to addition
        if capSizeT == 0 %when in shrinkage state
            capSizeT = max(0,increment);
        else %when in growth state
            if increment == -1
                capSizeT = 0;
            elseif increment == 0
                capSizeT = 1 - hydUnit;
            end
        end
    end
    
    %save MT length and cap size in arrays
    traj(iter,:) = [time mtLengthT capSizeT];

end

%remove first 100 seconds
traj = traj(vecLengthI+1:end,:);
traj(:,1) = traj(:,1) - traj(1,1);

%save data if requested
if saveTraj.saveOrNot
    if isempty(saveTraj.fileName)
        save(['mtTraj-',nowString],'traj'); %save in file
    else
        save(saveTraj.fileName,'traj'); %save in file (directory specified through name)
    end
end


%%%%% ~~ the end ~~ %%%%%

%subfunction called by mtGTPCapLDepK

function y = deriveParam(x);
%derives a few parameters

y = x;

if y.lDepend == 0
    y.kmin = y.kmax;
    y.lmin = -1;
    y.lmax = 1;
end
y.shift = mean([y.lmin y.lmax]);
y.coef = 4.6/(y.lmax-y.lmin);


%%%%% ~~ the end ~~ %%%%%

