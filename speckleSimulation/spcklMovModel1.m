function [params,fluorData,not2keep]=spcklMovModel1(params,fluorData,not2keep)
% stationary fluorophores with random poly/depoly at kon/koff

k01=params.k01; % kon,  G-actin to F-actin rate
k10=params.k10; % koff, F-actin to G-actin rate

% probabilities for switching from 1 to 0 and 0 to 1
P10=k10*params.dT;
P01=k01*params.dT;

% initialize temporary variables for state and coodinates
state=[fluorData.state(:,1) zeros(size(fluorData.state(:,1)))]; % nFluor x 2
py1=fluorData.py(:,1); % nFluor-vector
px1=fluorData.px(:,1); % nFluor-vector

% loop through every time point in the movie
for timePt=1:not2keep.nTmPtsInMovie
    % t is empty if current time point is not one of the ones to store;
    % otherwise it equals the current timePt.
    t=intersect(timePt,not2keep.pts2Store);  
    if ~isempty(t) % store the state/position info if t isn't empty
        col=find(not2keep.pts2Store==t); % column in which to store data
        fluorData.state(:,col)=state(:,1);
        fluorData.py(:,col)=py1;
        fluorData.px(:,col)=px1;
    end

    % where the state=0, use the prob. for switching to 1 (P01)
    % where the state=1, use the prob. for switching to 0 (P10)
    P=zeros(size(state));
    P(state(:,1)==0,:)=repmat([0 P01],sum(state(:,1)==0),1);
    P(state(:,1)==1,:)=repmat([P10 0],sum(state(:,1)==1),1);

    % initialize random number to compare against probability
    R=rand(size(state));
    
    % state switch occurs if the probability of switching is higher than a
    % random number from (0,1)
    PoverR=P./R; 
    [C I]=max(PoverR,[],2); % C gives value of largest P/R ratio; I gives the column number
    state(:,2)=state(:,1); % assign next time point to have same values
    state(C>1,2)=I(C>1)-1; % only change the ones where the max ratio is > 1
    
    % for this model there is no real need to store the state for both
    % timePt and timePt+1 or have separate columns for the probabilities in
    % P.  these become useful in situations where you have several
    % probabilities, like in model 3.
    
    % calculate next position with diffusion
    D = params.D*(10^9/10^2)^2; % convert D from cm^2/s to nm^2/s
   
    % assuming that the step taken in time step dT is normally distributed
    % with mean zero, get standard deviation of distribution from
    % diffusion constant D
    stepStd = sqrt(2*D*params.dT);
    %get particle's position at this time point
    py2 = py1 + stepStd*randn(length(py1),1);
    px2 = px1 + stepStd*randn(length(px1),1);
    
    
    % finally, assign the new state/position to the temporary variables for
    % calculating them at the next time step.
    px1=px2;
    py1=py2;
    state(:,1)=state(:,2);
end