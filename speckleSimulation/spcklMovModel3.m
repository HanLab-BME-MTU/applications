function [params,fluorData,not2keep]=spcklMovModel3(params,fluorData,not2keep)
% fluorophores move in one direction and change flow based on state

if params.protein == 1 % actin
    
    % actin can only go between 0 (unbound) and 2 (bound)
    k02=params.k02; % kon,  G-actin to F-actin rate
    k20=params.k20; % koff, F-actin to G-actin rate

    % all the other rate constants should be ZERO, since those state
    % changes aren't allowed for actin

    % kOFF from all states to unbound
    k10=0.00;
    k20=k20;
    k30=0.00; 

    % kON from all states to substrate bound
    k01=0.00;
    k21=0.00;
    k31=0.00; 

    % kON from all states to actin bound
    k02=k02;
    k12=0.00;
    k32=0.00; 

    % kON from all states to substrate and actin bound
    k03=0.00; 
    k13=0.00; 
    k23=0.00; 
    
    % rates of staying in present state is 1 minus sum of all the other rates
    k00=1-k01-k02-k03;
    k11=1-k10-k12-k13;
    k22=1-k20-k21-k23;
    k33=1-k30-k31-k32;

else % adhesion
    
    % kOFF from all states to unbound
    k10=params.k10;
    k20=params.k20;
    k30=params.k30;

    % kON from all states to substrate bound
    k01=params.k01;
    k21=params.k21;
    k31=params.k31;

    % kON from all states to actin bound
    k02=params.k02;
    k12=params.k12;
    k32=params.k32;

    % kON from all states to substrate and actin bound
    k03=params.k03;
    k13=params.k13;
    k23=params.k23;
    
    % rates of staying in present state is 1 minus sum of all the other rates
    k00=1-k01-k02-k03;
    k11=1-k10-k12-k13;
    k22=1-k20-k21-k23;
    k33=1-k30-k31-k32;

end

save([params.outputDirMovieInfo filesep 'rateConstants'],'k*');

%==============================================================

% CALCULATE PROBABILITIES

dT=params.dT;

P00=k00*dT;
P01=k01*dT;
P02=k02*dT;
P03=k03*dT;
P10=k10*dT;
P11=k11*dT;
P12=k12*dT;
P13=k13*dT;
P20=k20*dT;
P21=k21*dT;
P22=k22*dT;
P23=k23*dT;
P30=k30*dT;
P31=k31*dT;
P32=k32*dT;
P33=k33*dT;

% -----------------------------------------------------------
nmPerSecFlowSpeed=not2keep.nmPerSecFlowSpeed; % flow speed (nm/s)

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

    R=rand(size(state,1),4);
    P=zeros(size(state,1),4);
    P(state(:,1)==0,:)=repmat([P00 P01 P02 P03],sum(state(:,1)==0),1);
    P(state(:,1)==1,:)=repmat([P10 P11 P12 P13],sum(state(:,1)==1),1);
    P(state(:,1)==2,:)=repmat([P20 P21 P22 P23],sum(state(:,1)==2),1);
    P(state(:,1)==3,:)=repmat([P30 P31 P32 P33],sum(state(:,1)==3),1);

    PoverR=P./R;
    [C I]=max(PoverR,[],2); % C gives value of largest ratio and I the column number
    state(:,2)=state(:,1); % assign next time point to have same values
    state(C>1,2)=I(C>1)-1; % only change the ones where the max ratio is > 1

    v=zeros(size(state,1),1);
    v(state(:,1)==0)=0;                      % unbound
    v(state(:,1)==1)=0;                      % substrate
    v(state(:,1)==2)=nmPerSecFlowSpeed;      % actin
    v(state(:,1)==3)=nmPerSecFlowSpeed/2;    % both

    % calculate next position with diffusion
    D = params.D*(10^9/10^2)^2; % convert D from cm^2/s to nm^2/s
    driftVel = [v params.theta*ones(length(v),1)]; % [drift speed, angle in degrees]
    theta = driftVel(:,2)*pi/180; % convert to radians
    driftVelCart = repmat(driftVel(:,1),[1 2]).*[sin(theta) cos(theta)]; % cartesian [vy vx]
    % assuming that the step taken in time step dT is normally distributed
    % with mean zero, get standard deviation of distribution from
    % diffusion constant D
    stepStd = sqrt(2*D*params.dT);
    %get particle's position at this time point
    py2 = py1 + driftVelCart(:,1)*params.dT + stepStd*randn(length(py1),1);
    px2 = px1 + driftVelCart(:,2)*params.dT + stepStd*randn(length(px1),1);
    
    % finally, assign the new state/position to the temporary variables for
    % calculating them at the next time step.
    px1=px2;
    py1=py2;
    state(:,1)=state(:,2);
end

