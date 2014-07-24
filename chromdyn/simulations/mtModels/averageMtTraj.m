function [mtLengthAve,mtLengthSD,errFlag] = averageMtTraj(mtLength,dt,...
    expTimeStep,aveInterval)
%AVERAGEMTTRAJ finds the average and stad. dev. of MT positions over experimental intervals
%
%The function takes the value of MT length at every simulation time step, 
%determines the appropriate intervals that must be averaged over to get
%values corresponding to experimental data, then get the average and standard
%deviation of microtubule length at experimental time points.
%
%Synopsis [mtLengthAve,mtLengthSD,errFlag] = averageMtTraj(mtLength,dt,...
%    expTimeStep,aveInterval)
%
%INPUT  mtLength    : Length of microtubule as a function of time.
%       dt          : Time step used in Monte Carlo simulation.
%       expTimeStep : Time step used in experimental measurements.
%       aveInterval : Interval used for averaging (<=expTimeStep).
%
%OUTPUT mtLengthAve : Average length of microtubule at experimental time
%                     points.
%       mtLengthSD  : Standard deviation of microtubule length in interval
%                     used for averaging about time points.
%       errFlag     : 0 if function runs normally, 1 otherwise.
%
%Khuloud Jaqaman, 09/03

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('averageMtTraj');
    disp('--averageMtTraj: Incorrect number of input arguments!');
    errFlag = 1;
    return;
end

%check for error in input data

% if max(max(mtLength<=0))
%     disp(' ');
%     disp('Microtubule length should be always positive!');
%     errFlag = 1;
% end
if dt <= 0
    disp('--averageMtTraj: Simulation time step should be positive!');
    errFlag = 1;
end
if expTimeStep <= 0
    disp('--averageMtTraj: Experimental time step should be positive!');
    errFlag = 1;
end
if aveInterval <= 0
    disp('--averageMtTraj: Averaging Interval should be positive!');
    errFlag = 1;
end
if aveInterval > expTimeStep
    disp('--averageMtTraj: Averaging interval should not be larger than experimental time step!');
    errFlag = 1;
end
if dt >= aveInterval
    disp('--averageMtTraj: Simulation time step should be smaller than averaging interval!');
    errFlag = 1;
end
if errFlag
    disp('--averageMtTraj: Please fix input data!');
    return;
end
    
ratio = round(expTimeStep/dt);          %ratio of experimental time step to 
                                        %simulation time step, rounded to
                                        %the nearest integer. 

vectorLength = length(mtLength);        %number of entries in mtLength

trash = rem(vectorLength,ratio);        %number of entries to be removed from
                                        %end of mtLength to make its length an 
                                        %integer multiple of ratio 
                                        
mtLength = mtLength(1:vectorLength-trash); %remove extra entries from the end

vectorLength = length(mtLength);        %new number of entries in mtLength

trajLength = floor(vectorLength/ratio); %num. of entries in "experimental" trajectory

trajectories = reshape(mtLength,ratio,trajLength);    %reshape mtLength to get
                                                      %"ratio" trajectories of
                                                      %"trajLength" entries each

ratio = round(aveInterval/dt);          %ratio of averaging interval to 
                                        %simulation time step, rounded to 
                                        %nearest integer.
                                                      
trajectories = trajectories(1:ratio,:); %real interval that is averaged over.
                                                     
mtLengthAve = mean(trajectories)';     %average lengthof MT at each interval
mtLengthSD  = max(1e-15,abs(std(trajectories)')); %stand. dev. of MT length at each interval
