function [tagPos,errFlag] = hingeModel(hingeParam,coordInit,totalTime,dt);
%HINGEMODEL performs the free/diffusive rotation of tag about MT tip.
%
%SYNOPSIS [tagPos,errFlag] = hingeModel(hingeParam,coordInit,totalTime,dt);
%
%INPUT  hingeParam : Parameters used in the model:
%           free        : 0 if diffusive rotation, 1 if free rotation.
%           chromL      : Long semi-axis of chromosome (microns).
%           chromS      : Short semi-axis of chromosome (microns).
%           viscosity   : viscosity of medium inside the nucleus (Pa*s).
%           temperature : Absolute temperature of system (K).
%       coordInit  : Initial position of hinge tip, in Cartesian
%                    coordinates, entered as a row vector (microns).
%       totalTime  : Total time of simulation in seconds.
%       dt         : Time step used for time discretization. 
%
%OUTPUT tagPos  : position of hinge tip at every time step.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%COMMENTS For now, radius is assumed constant. So the hinge only rotates. 
%
%Khuloud Jaqaman, 12/03

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('hingeModel');
    disp('--hingeModel: Incorrect number of input arguments!');
    errFlag = 1;
    return;
end

free = hingeParam.free;
if ~free
    %     chromL = hingeParam.chromL;
    %     chromS = hingeParam.chromS;
    %     viscosity = hingeParam.viscosity;
    %     temperature = hingeParam.temperature;
    diffConst = hingeParam.diffConst;
end

%check for error in input
if free ~= 0 && free ~= 1
    disp('--hingeModel: The variable free should be wither 0 or 1');
    tagPos = [];
    errFlag = 1;
    return;
end
% if ~free
%     %     if isempty(chromL) || chromL < 0
%     %         disp('--hingeModel: Large radius of chromosome should be positive!');
%     %         errFlag = 1;
%     %     end
%     %     if isempty(chromS) || chromS < 0
%     %         disp('--hingeModel: Small radius of chromosome should be positive!');
%     %         errFlag = 1;
%     %     end
%     %     if isempty(viscosity) || viscosity < 0
%     %         disp('--hingeModel: Viscosity should be positive!');
%     %         errFlag = 1;
%     %     end
%     %     if isempty(temperature) || temperature < 0
%     %         disp('--hingeModel: temperature should be positive!');
%     %         errFlag = 1;
%     %     end
% end
if totalTime <= 0
    disp('--hingeModel: Total time should be positive!');
    errFlag = 1;
end
if dt <= 0
    disp('--hingeModel: Time step should be positive!');
    errFlag = 1;
end
if errFlag
    disp('--hingeModel: Please fix input data!');
    tagPos = [];
    return;
end

if ~free
    
%     %rotational diffusion constant of a chromosome about its two short semi-axes, 
%     %assuming it is a long prolate ellipsoid of revolution with semi-axes chromL 
%     %and chromS (chromL > 5*chromS). The formula for the rotational friction coefficient
%     %is taken from [C. Tanford, Physical Chemistry of Macromolecules (1961), p.436]. 
%     %The formula also uses the viscosity of the medium.
%     chromL = chromL*1e-6; %in meters
%     chromS = chromS*1e-6; %in meters
%     frictionCoef = 16*pi*viscosity*chromL^3/(-1+2*log(2*chromL/chromS))/3; %kg*m^2/s
%     diffConst = 1.38e-23*temperature/frictionCoef; %s^-1
    
    %maximum anglular displacement per time step based on the relation <dtheta^2> = 2Ddt
    maxAngle = sqrt(2*diffConst*dt);
    
end

%distance between tag and MT tip. Must maintain the same value throughout the simulation
radius = sqrt(sum(coordInit.^2)); 

vecLength = floor(totalTime/dt);   %approximate length of tagPos 
tagPos = zeros(vecLength,3);       %allocate memory for array tagPos

tagPos(1,:) = coordInit;   %initial position of hinge tip in Cartesian coordinates

time = 0;   %initialization of loop
iter = 1;

if free     %free rotation, no diffusion constraints
    while time < totalTime

        iter = iter + 1;    %update iteration number
        time = time + dt;   %update time
        
        theta = acos(2*rand-1); %new theta
        phi = 2*pi*rand;        %new phi
        
        tagPos(iter,1) = radius*sin(theta)*cos(phi);    %new position in Cartesian coordinates
        tagPos(iter,2) = radius*sin(theta)*sin(phi);
        tagPos(iter,3) = radius*cos(theta);
        
    end
else       %diffusion limited rotation
    while time < totalTime

        x = tagPos(iter,1);
        y = tagPos(iter,2);
        z = tagPos(iter,3);
                
        iter = iter + 1;               %update iteration number
        time = time + dt;              %update time
        
        theta = (2*rand-1)*maxAngle; %rotation in time step about one short axis
        phi   = (2*rand-1)*maxAngle; %rotation in time step about other short axis
        
        cosTheta = cos(theta); %sines and cosines of angles 
        sinTheta = sin(theta);
        cosPhi = cos(phi);
        sinPhi = sin(phi);
        
        dummy = -sinTheta*y + cosTheta*z;
        tagPos(iter,1) = cosPhi*x - sinPhi*dummy; %get new Cartesian coordinates
        tagPos(iter,2) = cosTheta*y + sinTheta*z;
        tagPos(iter,3) = sinPhi*x + cosPhi*dummy;
        
    end
end
