%loadParam outputs parameter structure used by positionSimulation
%
% This function is used to load the param structure required for
% simulations of motor driven vesicular transport.
%
% INPUT         none
%
% OUTPUT        param: Structure of parameters used in simulation
%
% DEPENDENCES   loadParam calls no function
%               oneParameterScan is used by parameterScan 
%
% REMARKS
%
% Created November 7, 2007 by DA Nunez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assign default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.Force = 0;

%General simulation parameters
param.time_step = 0.01; %time step of simulation in seconds
param.time_max = 10; %total time of simulation run
param.sampling_rate = 0.1; %frame sampling rate in seconds
param.vesicle_number = 1;   %number of vessicles in simulation
param.microtubule_length = 32256; %length of microtubule in field of view in nm
param.total_kinesin = 2; %total number of kinesin on vesicle
param.total_dynein = 0; %total number of dynein on vesicle

%Run Tracking Variables
param.result_number = '001';
param.id = 'KinesinFVsimulation';
param.result_id = [param.id '_' num2str(param.vesicle_number) 'vesicles' param.result_number];

%Probabilities for Kinesin state changes
param.kinesin_probability_matrix = nan(3);
param.kinesin_probability_matrix(1,2) = 0.1;
param.kinesin_probability_matrix(1,3) = 0.01;
param.kinesin_probability_matrix(3,1) = 0.35;
param.kinesin_probability_matrix(2,3) = 0.2;
param.kinesin_probability_matrix(3,2) = 0.5;
param.kinesin_probability_matrix(2,2) = 0.1;

% Porbabilities for Dynein state changes
param.dynein_probability_matrix = nan(3);
param.dynein_probability_matrix(1,2) = 0.3;
param.dynein_probability_matrix(1,3) = 0.3;
param.dynein_probability_matrix(3,1) = 0.3;
param.dynein_probability_matrix(2,1) = 0.3;

%Force/Velocity constants for Kinesin
param.velocity_max_kinesin = 1100; %velocity under zero load in nm/sec
param.Force_stall_kinesin = 7; %stall force in pN
param.kinesinRestingLength = 65;  %resting length of kinesin tether in nanometers 65 nm as described in Atzberger and Peskin (2006).
param.kinesinTetherStiffness = 0.2; %the stiffness of the kinesin tether (rough estimate from figure 7 in Atzberger and Peskin (2006).

%Force/Velocity constants for Dynein
param.velocity_max_dynein = 1000; %velocity under zero load in nm/sec
param.Force_stall_dynein = 1; %stall force in pN
param.dyneinRestingLength = 65;
param.tetherStiffnessRatio = 1; %ratio of the stiffness of the dynein tether over the stiffness of the kinesin tether
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%