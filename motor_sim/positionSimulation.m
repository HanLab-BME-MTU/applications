function [param,tracks, motorInfo] = positionSimulation(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%positionSimulation takes param and outputs tracks structure from
%simulation
%
% This function carries out the simulation of the motors walking on a
% microtubule and outputs the tracks structure for each vesicle, which
% includes x-positions, y-positions, and intensities for each time point.
% The function takes the structure param as input, which contains all the
% simulation parameters (see loadParam.m). The motor can be in one of three
% states: walking, paused, or unbound. The vesicle itself never unbinds from the
% microtubule track, and therefore there is no diffusion.
%
% INPUT
%               param             : simulation parameter structure
%
% OUTPUT        This function saves the currentParam and tracks
%               structures for each simulation iteration, the param (input
%               into the function), the parameterData (parameter values,
%               mean of average run velocities, and std of average run
%               velocities), and a plot of the mean of the average run
%               velocites plotted against value of the parameter scanned
%               through.
%
% DEPENDENCES
%               positionSimulation is used by parameterScan
%
% REMARKS
%
% Created November 9, 2007 by DA Nunez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%calculate probabilities not included in param structure
param.kinesin_probability_matrix(1,1) = 1 - param.kinesin_probability_matrix(1,2) - param.kinesin_probability_matrix(1,3);
param.kinesin_probability_matrix(2,1) = 1 - param.kinesin_probability_matrix(2,2) - param.kinesin_probability_matrix(2,3);
param.kinesin_probability_matrix(3,3) = 1 - param.kinesin_probability_matrix(3,1) - param.kinesin_probability_matrix(3,2);
param.dynein_probability_matrix(1,1) = 1 - param.dynein_probability_matrix(1,2) - param.dynein_probability_matrix(1,3);
param.dynein_probability_matrix(2,1) = 1 - param.dynein_probability_matrix(2,2) - param.dynein_probability_matrix(2,3);
param.dynein_probability_matrix(3,3) = 1 - param.dynein_probability_matrix(3,1) - param.dynein_probability_matrix(3,2);

%save local copy of probability matrix that will be edited multiple times
%during looping
kinesin_probability_matrix = param.kinesin_probability_matrix;
dynein_probability_matrix = param.dynein_probability_matrix;

%initialize array of structures tracks as taken by Ge's programs, with the
%state counters for kinesin and dynein added (which are used to ensure that
%state transition probabilities are consistent with input probabilities
tracks(1:param.vesicle_number) = struct('startID',1,'len',param.time_max/param.sampling_rate,'points',nan(param.time_max/param.sampling_rate,3));

%initiate variables and allocate space
iter = 2; %iteration counter start; iter = 1 are initial conditions, so looping really starts on second loop
vesiclePosition = nan(1,param.vesicle_number,param.time_max/param.time_step); %this matrix will store the positions of all the vesicles for each time point
vesiclePosition(1,1:param.vesicle_number,1) = param.microtubule_length*rand(param.vesicle_number,1); %starting position in nm is chosen randomly along microtubule length

%start the kinesin and dynein current state variables that store the state
%of each kinesin and dynein on each vesicle (the second dimension) for each
%time step (the third dimention)
kinesin_current_state_number = zeros(param.total_kinesin,param.vesicle_number,param.time_max/param.time_step);
dynein_current_state_number = zeros(param.total_dynein,param.vesicle_number,param.time_max/param.time_step);
kinesin_current_state_number(:,:,1) = min(ceil(3*rand(param.total_kinesin,param.vesicle_number)+0.0000001),repmat(3,param.total_kinesin,param.vesicle_number)); %starting state of each kinesin is randomly chosen
dynein_current_state_number(:,:,1) = min(ceil(3*rand(param.total_dynein,param.vesicle_number)+0.0000001),repmat(3,param.total_kinesin,param.vesicle_number)); %starting state of each dynein is randomly chosen

kinesin_state_decision_number = zeros(param.total_kinesin,param.vesicle_number); %this matrix will hold the values of the states that each kinesin from each vesicle decides to change to
dynein_state_decision_number = zeros(param.total_dynein,param.vesicle_number); %this matrix will hold the values of the states that each dynein from each vesicle decides to change to
kinesinPositionChange = zeros(param.total_kinesin,1); %this column vector will hold the change in tether length for each kinesin on a given vesicle; since the length values are stored on the ThetherLength matrix, this particular vector is erased for each vesicle during looping.
dyneinPositionChange = zeros(param.total_dynein,1); %this column vector will hold the change in tether length for each dynein on a given vesicle; since the length values are stored on the ThetherLength matrix, this particular vector is erased for each vesicle during looping.

%tether lengths will be stored for each unit time (3rd dimension) for each
%vesicle(2nd dimension) and for each motor on these vesicles (1st
%dimension)
kinesinTetherLength = zeros(param.total_kinesin,param.vesicle_number,param.time_max/param.time_step);
dyneinTetherLength = zeros(param.total_dynein,param.vesicle_number,param.time_max/param.time_step);
%starting tether lengths will be the resting length at time t = 0 (iter=1)
kinesinTetherLength(1:param.total_kinesin,1:param.vesicle_number,1) = param.kinesinRestingLength; %kinesin tether starts at resting length
dyneinTetherLength(1:param.total_dynein,1:param.vesicle_number,1) = param.dyneinRestingLength; %dynein tether starts at resting length
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main simulation loop over time
while param.time_step*(iter - 1) <= param.time_max

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This loop carries out simulation for each vesicle
    for vesicle = 1:param.vesicle_number

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This loop decides state and position for each kinesin on vesicle
        if param.total_kinesin ~= 0
            for kinesin = 1:param.total_kinesin

                %Recalculate Force on kinesin based on F =
                %springConstant*displacementFromEquilibrium
                if kinesinTetherLength(kinesin,vesicle,iter-1) >= param.kinesinRestingLength
                    kinesinForce = param.kinesinTetherStiffness*(kinesinTetherLength(kinesin,vesicle,iter-1) - param.kinesinRestingLength)/param.total_kinesin;
                else
                    kinesinForce = 0;
                end

                %Recalculate kinesin probability matrix values that change with
                %changing force
                kinesin_probability_matrix(1,3) = param.kinesin_probability_matrix(1,3)*exp(kinesinForce/param.Force_stall_kinesin);
                kinesin_probability_matrix(1,2) = param.kinesin_probability_matrix(1,2)*exp(kinesinForce/param.Force_stall_kinesin);
                kinesin_probability_matrix(2,3) = param.kinesin_probability_matrix(2,3)*exp(kinesinForce/param.Force_stall_kinesin);
                kinesin_probability_matrix(2,2) = param.kinesin_probability_matrix(2,2)*exp(kinesinForce/param.Force_stall_kinesin);
                kinesin_probability_matrix(1,1) = max(1 - kinesin_probability_matrix(1,2) - kinesin_probability_matrix(1,3),0);
                kinesin_probability_matrix(2,1) = max(1 - kinesin_probability_matrix(2,2) - kinesin_probability_matrix(2,3),0);
                kinesin_probability_matrix(3,3) = max(1 - kinesin_probability_matrix(3,1) - kinesin_probability_matrix(3,2),0);

                %CHECK:prevents probabilities from being negative

                if any(kinesin_probability_matrix(:) + 0.00001 < 0)
                    error('probability less than zero');
                    return
                end

                %decide to which state kinesin going to change based on current state
                kinesin_decision_rand = rand(1);
                if kinesin_decision_rand <= kinesin_probability_matrix(kinesin_current_state_number(kinesin,vesicle,iter-1),1)
                    kinesin_state_decision_number(kinesin,vesicle) = 1;
                elseif kinesin_decision_rand <= kinesin_probability_matrix(kinesin_current_state_number(kinesin,vesicle,iter-1),1) + kinesin_probability_matrix(kinesin_current_state_number(kinesin,vesicle,iter-1),2)
                    kinesin_state_decision_number(kinesin,vesicle) = 2;
                else
                    kinesin_state_decision_number(kinesin,vesicle) = 3;
                end %kinesin decision loop

                %set position as resting length from vesicle position if
                %dissociating or determine change in kinesin position (no
                %change if pausing; change in position determined by
                %force-velocity curve for given time step if walking)
                if kinesin_current_state_number(kinesin,vesicle,iter-1) == 1 | kinesin_current_state_number(kinesin,vesicle,iter-1) == 2
                    if kinesin_state_decision_number(kinesin,vesicle) == 3
                        kinesinPositionChange(kinesin) = -kinesinTetherLength(kinesin,vesicle,iter-1)+param.kinesinRestingLength;
                    elseif kinesin_state_decision_number(kinesin,vesicle) == 1
                        kinesinVelocity = param.velocity_max_kinesin;
                        kinesinPositionChange(kinesin) = kinesinVelocity*param.time_step;
                    else
                        kinesinPositionChange(kinesin) = 0;
                    end %kinesin tether length loop
                %if binding to microtubule pick a random position to bind;
                %should maybe allow binding further than resting length
                %with some sort of decreasing probability distribution
                elseif kinesin_current_state_number(kinesin,vesicle,iter-1) == 3
                    if kinesin_state_decision_number(kinesin,vesicle) == 1 | kinesin_state_decision_number(kinesin,vesicle) == 2
                        kinesinPositionChange(kinesin) = -param.kinesinRestingLength*rand(1);
                    else
                        kinesinPositionChange(kinesin) = 0;
                    end
                end %of if attached
                kinesinTetherLength(kinesin,vesicle,iter) = kinesinTetherLength(kinesin,vesicle,iter - 1) + kinesinPositionChange(kinesin);
            end %kinesin for loop
        end %kinesin if loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This loop decides state and position for each dynenin on vesicle
      if param.total_dynein ~= 0
            for dynein = 1:param.total_dynein

                %Recalculate Force on dynein based on F =
                %springConstant*displacementFromEquilibrium
                if dyneinTetherLength(dynein,vesicle,iter-1) >= param.dyneinRestingLength
                    dyneinForce = param.dyneinTetherStiffness*(dyneinTetherLength(dynein,vesicle,iter-1) - param.dyneinRestingLength)/param.total_dynein;
                else
                    dyneinForce = 0;
                end

                %Recalculate dynein probability matrix values that change with
                %changing force
                dynein_probability_matrix(1,3) = param.dynein_probability_matrix(1,3)*exp(dyneinForce/param.Force_stall_dynein);
                dynein_probability_matrix(1,2) = param.dynein_probability_matrix(1,2)*exp(dyneinForce/param.Force_stall_dynein);
                dynein_probability_matrix(2,3) = param.dynein_probability_matrix(2,3)*exp(dyneinForce/param.Force_stall_dynein);
                dynein_probability_matrix(2,2) = param.dynein_probability_matrix(2,2)*exp(dyneinForce/param.Force_stall_dynein);
                dynein_probability_matrix(1,1) = max(1 - dynein_probability_matrix(1,2) - dynein_probability_matrix(1,3),0);
                dynein_probability_matrix(2,1) = max(1 - dynein_probability_matrix(2,2) - dynein_probability_matrix(2,3),0);
                dynein_probability_matrix(3,3) = max(1 - dynein_probability_matrix(3,1) - dynein_probability_matrix(3,2),0);

                %CHECK:prevents probabilities from being negative

                if any(dynein_probability_matrix(:) + 0.00001 < 0)
                    error('probability less than zero');
                    return
                end

                %decide to which state dynein going to change based on current state
                dynein_decision_rand = rand(1);
                if dynein_decision_rand <= dynein_probability_matrix(dynein_current_state_number(dynein,vesicle,iter-1),1)
                    dynein_state_decision_number(dynein,vesicle) = 1;
                elseif dynein_decision_rand <= dynein_probability_matrix(dynein_current_state_number(dynein,vesicle,iter-1),1) + dynein_probability_matrix(dynein_current_state_number(dynein,vesicle,iter-1),2)
                    dynein_state_decision_number(dynein,vesicle) = 2;
                else
                    dynein_state_decision_number(dynein,vesicle) = 3;
                end %dynein decision loop

                %set position as resting length from vesicle position if
                %dissociating or determine change in dynein position (no
                %change if pausing; change in position determined by
                %force-velocity curve for given time step if walking)
                if dynein_current_state_number(dynein,vesicle,iter-1) == 1 | dynein_current_state_number(dynein,vesicle,iter-1) == 2
                    if dynein_state_decision_number(dynein,vesicle) == 3
                        dyneinPositionChange(dynein) = -dyneinTetherLength(dynein,vesicle,iter-1)+param.dyneinRestingLength;
                    elseif dynein_state_decision_number(dynein,vesicle) == 1
                        dyneinVelocity = param.velocity_max_dynein;
                        dyneinPositionChange(dynein) = dyneinVelocity*param.time_step;
                    else
                        dyneinPositionChange(dynein) = 0;
                    end %dynein tether length loop
                %if binding to microtubule pick a random position to bind;
                %should maybe allow binding further than resting length
                %with some sort of decreasing probability distribution
                elseif dynein_current_state_number(dynein,vesicle,iter-1) == 3
                    if dynein_state_decision_number(dynein,vesicle) == 1 | dynein_state_decision_number(dynein,vesicle) == 2
                        dyneinPositionChange(dynein) = -param.dyneinRestingLength*rand(1);
                    else
                        dyneinPositionChange(dynein) = 0;
                    end
                end %of if attached
                dyneinTetherLength(dynein,vesicle,iter) = dyneinTetherLength(dynein,vesicle,iter - 1) + dyneinPositionChange(dynein);
            end %dynein for loop
        end %dynein if loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %Calculate change of positon of bead based on force balance in
        %between dynein and kinesin tether forces.
        if param.total_kinesin == 0
            param.dyneinTetherStiffness = param.kinesinTetherStiffness*param.tetherStiffnessRatio;
            vesiclePositionChange = -param.Force/param.dyneinTetherStiffness/param.total_dynein - sum(dyneinTetherLength(:,vesicle,iter))/param.total_dynein+param.dyneinRestingLength;
        elseif param.total_dynein == 0
            vesiclePositionChange = -param.Force/param.kinesinTetherStiffness/param.total_kinesin + sum(kinesinTetherLength(:,vesicle,iter))/param.total_kinesin-param.kinesinRestingLength;
        else
            vesiclePositionChange = (sum(kinesinTetherLength(:,vesicle,iter)) - param.total_kinesin*param.kinesinRestingLength...
                - (sum(dyneinTetherLength(:,vesicle,iter)) - param.total_dynein*param.dyneinRestingLength)*param.tetherStiffnessRatio...
                - param.Force/param.kinesinTetherStiffness)/(param.total_kinesin + param.total_dynein*param.tetherStiffnessRatio);
        end
        %score changes in position
        if any(kinesinTetherLength(:,vesicle,iter) > param.kinesinRestingLength)
            vesiclePosition(1,vesicle,iter) = vesiclePosition(1,vesicle,iter-1) + vesiclePositionChange;
        else
            vesiclePosition(1,vesicle,iter) = vesiclePosition(1,vesicle,iter-1);
            vesiclePositionChange = 0;
        end
        %change state
        kinesin_current_state_number(1:param.total_kinesin,vesicle,iter) = kinesin_state_decision_number(1:param.total_kinesin,vesicle);
        dynein_current_state_number(1:param.total_dynein,vesicle,iter) = dynein_state_decision_number(1:param.total_dynein,vesicle);

        %recalculate tether lengths taking into account change in position
        %of vesicle; if the motor is not attached (state #3) then the
        %tether length should remain at the resting length and should not
        %take part in the force calculation at the beginning of the vesicle
        %loop
        for kinesin = 1:param.total_kinesin
            if kinesin_state_decision_number(kinesin,vesicle) == 3
                %do nothing
            else
                kinesinTetherLength(kinesin,vesicle,iter) = kinesinTetherLength(kinesin,vesicle,iter) - vesiclePositionChange;
            end %kinesin tether length recalculation (due to movement of vesicle) conditional loop
        end %kinesin tether length recalculation loop
        for dynein = 1:param.total_dynein
            if dynein_state_decision_number(dynein,vesicle) == 3
                %do nothing
            else
                dyneinTetherLength(dynein,vesicle,iter) = dyneinTetherLength(dynein,vesicle,iter) + vesiclePositionChange;
            end %dynein tether length recalculation (due to movement of vesicle) conditional loop
        end %dynein tether length recalculation loop

    end %vesicle loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %update iteration count
    iter = iter + 1;
end %time point loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%store kinesin positions and states for each time step in a structure
%called kinesinInfo
if param.total_kinesin ~= 0
    motorInfo.kinesinState = kinesin_current_state_number;
    motorInfo.kinesinPosition = (repmat(vesiclePosition,param.total_kinesin,1) + kinesinTetherLength); %in nanometers
end
%store dynein positions and states for each time step in a structure
%called kinesinInfo
if param.total_dynein ~= 0
    motorInfo.dyneinState = dynein_current_state_number;
    motorInfo.dyneinPosition = (repmat(vesiclePosition,param.total_dynein,1) - dyneinTetherLength); %in nanometers
end

%save all vesicle position data for plotting with motor data
motorInfo.vesiclePosition = vesiclePosition; % in nanometers

%take position data according to sampling rate
sampledVesiclePosition(1:param.vesicle_number,1:param.time_max/param.sampling_rate) = vesiclePosition(1,1:param.vesicle_number,1:param.sampling_rate/param.time_step:param.time_max/param.time_step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This loop populates tracks parameter data
for n = 1:param.vesicle_number
    tracks(n).points(:,2) = 0;
    tracks(n).points(:,3) = 0.0000001;
    tracks(n).points(:,1) = sampledVesiclePosition(n,:);%/1000/0.126'; %position data is in pixels
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
