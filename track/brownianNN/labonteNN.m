function [M,conflict,options,v,w]=labonteNN(x,y,options)
% labonteNN performs nearest-neighbor matching by neural network algorithms (self-organizing maps)
%
% SYNOPSIS   [M,conflict,options,v,w]=labonteNN(x,y,options)
%
% INPUT      x        :   coordinate matrix [y x]m for frame 1 (m particles)
%            y        :   coordinate matrix [y x]n for frame 2 (n particles)
%            options  :   structure of parameters (optional)
%                            options.r         : initial maximum distance under which a neuron can 
%                                                be activated [ default: max found coordinate in {y,x} ]
%                                                Pass options.r == 0 to have the function select a radius 
%                                                based on the data (3x the mean minimal distance among particles)
%                            options.rf        : difference of the criterion of two consecutive 
%                                                iterations which stops the process [ default: 5e-3 ]
%                            options.alpha     : initial learning rate [ default: 5e-4 ]
%                            options.beta      : value by which .r is decreased and .alpha 
%                                                is increased at each iteration [ default: 0.99 ]
%                            options.threshold : max accepted distances between two particles in the 
%                                                final result (not necessarily equal to the search radius)
%                            options.maxiter   : maximum number of iterations; in case the algorithm gets stuck 
%                                               [default: 1000]
%                                                
% OUTPUT     M        :   matrix of matched positions [y x y x]
%            conflict :   matrix of still unresolved positions [y x y x]
%                         (with repetitions)
%            options  :   used options structure
%            v        :   x coordinates as tranformed during the process (final weights)
%            w        :   y coordinates as tranformed during the process (final weights)
%
% Should some positions fail to be resolved appropriately (this is a rare event if the options are 
% set carefully), labonteNN returns the conflicting positions in a separate matrix ('conflict'), 
% which can be passed for a second iteration to labonteNN. The output of the second iteration can 
% then be appended to the one of the first.
%
% DEPENDENCES   labonteNN uses { createDistanceMatrix ;
%                                createSparseDistanceMatrix;
%                                missingIndices;
%                                repeatingIndices }
%               fsmTrackMain is used by { fsmTrackTrackerBMTNN }
%
% REFERENCES 
%
% (1) Labonté G. (1998) "A SOM Neural Network That Reveals Continuous Displacement Fields" 
%     In the Proceedings of the 1998 International Joint Conference on Neural Networks 
%     at WCCI'98 (the World Congress on Computational Intelligence), held at Anchorage Alaska,
%     published by IEEE, ISBN 0-7803-4862-1.
%
% (2) Labonté G. (1999) "A New Neural Network for Particle-Tracking Velocimetry" Experiments 
%     in Fluids, 26 (4): 340-346.
%
% (3) Labonté G. (2000) "On a neural network that performs an enhanced nearest-neighbour matching"
%     PATTERN ANALYSIS AND APPLICATIONS 3 (3): 267-278 2000
%
% Aaron Ponti, 2002 - Major update 04/16/2004

% FLAGS
DEBUG=0;    % Plots the particles at every iteration to see convergence progress
VERBOSE=1;  % Outputs some information about the current iteration (text only)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT CHECK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default options if none is passed
if nargin==2
    options.r=max([max(x(:)) max(y(:))]);
    options.rf=5e-3; 
    options.alpha=5e-4;
    options.beta=0.99;
    options.threshold=1;
    options.maxiter=1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOME OUTPUT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DEBUG==1
    
    mn=min([min(x(:)) min(y(:))])-0.001;
    mx=max([max(x(:)) max(y(:))])+0.001;
    figure;
    plot(x(:,1),x(:,2),'k.','MarkerSize',20);
    hold on;
    plot(y(:,1),y(:,2),'r.','MarkerSize',20);
    hold on;
    axis([mn mx mn mx]);
    title('Starting positions [black: source, red: target]');
    
    figure;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOME INITIALIZATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a constant
nDim=size(x,2);     % Number of dimensions

% The initial weights are set to the particle positions
v=x;
w=y;

% Initial search radius
r=options.r;
if r==0
    % The user decided to let this function pick the initial search radius
    r=3*mean(min(createDistanceMatrix(v,w),[],2));
    % Check r, in case all particles in frame 1 have a perfect match at the exact same position in frame 2
    if r==0
        r=1;
    end
end

alpha=options.alpha;

% Initialize weight matrices dw and dv
dw=zeros(size(w));
dv=zeros(size(v));

% Initialize some values
crit=1;          % Current criterion (initially set larger than the threshold)
critThreshold=0; % Criterion threshold, set to zero (will be updated after the first iteration)
lastCrit=0;      % Criterion of the previous iteration (initially zero)

% Maximum number of iterations
maxIter=options.maxiter;

% Current iteration number 
nIter=0;

% Iterate
while crit>critThreshold & nIter<maxIter
    
    % Reset current dw and dv
    dw=0*dw;
    dv=0*dv;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % FEED NETWORK 1 WITH THE WEIGHTS FROM NETWORK 2
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the distance of each neuron to all others
    d=createDistanceMatrix(v,w);
    
    % For all neurons from network 1 find the winning neurons in network 2
    [minD,vN]=min(d,[],2);
    
    % The winning neuron is the one (or one of the ones) having the smallest activation (distance)
    distN=v-w(vN,:);

    % Mark those neurons which have to be updated (distance criterion)
    neuronsToUpdate=createSparseDistanceMatrix(w,w(vN,:),r)~=0;
    
    % Calculate total weights for participating neurons in network 2 (it is just a matrix multiplication)
    dw=(alpha*distN'*neuronsToUpdate')';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % FEED NETWORK 2 WITH THE WEIGHTS FROM NETWORK 1
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate the distance of each neuron to all others (just the transpose of the previous)
    d=d';
    
    % For all neurons from network 2 find the winning neurons in network 1
    [minD,vN]=min(d,[],2);
    
    % The winning neuron is the one (or one of the ones) having the smallest activation (distance)
    distN=w-v(vN,:);
    
    % Mark those neurons which have to be updated (distance criterion)
    neuronsToUpdate=createSparseDistanceMatrix(v,v(vN,:),r)~=0;
    
    % Calculate total weights for participating neurons in network 1 (it is just a matrix multiplication)
    dv=(alpha*distN'*neuronsToUpdate')';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % UPDATE WEIGHTS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    w=w+dw;
    v=v+dv;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % UPDATE TIME-DEPENDENT PARAMETERS (AND ITERATION COUNTER)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=r*options.beta;
    alpha=alpha/options.beta;
    
    % Update iteration counter
    nIter=nIter+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CALCULATE CRITERION
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate distance criterion (sum of squared distance differences)
    minDist=min(createDistanceMatrix(v,w),[],2);
    sqErr=sum(minDist.^2);
    crit=abs(sqErr-lastCrit);
    lastCrit=sqErr;
    
    % Use the criterion at iteration 2 to update the threshold
    if nIter==2
        critThreshold=crit*options.rf;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % SOME DEBUG INFO
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if DEBUG==1
        mx=max([max(x(:)) max(y(:))])+0.001;
        plot(v(:,1),v(:,2),'r*','MarkerSize',10);
        %     axis([0 mx 0 mx]);
        hold on;
        plot(w(:,1),w(:,2),'b*','MarkerSize',10);
        %     axis([0 mx 0 mx]);
        title('Convergence');
        pause(0.001);
    end
    
    if VERBOSE==1
        fprintf(1,'Iteration number: %4.0d - Criterion: %5.5f - Threshold: %.5f - Current radius: %3.3f\n',nIter,crit,critThreshold,r);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOOK FOR CORRESPONDANCES IN THE RETURNED CONVERGED POSITIONS (THEY SHOULD BE UNIVOCALLY DETERMINED)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build pairing matrix
D=createDistanceMatrix(v,w); % We don't set a distance constraint here
Dorig=createSparseDistanceMatrix(x,y,options.threshold); % Only in the original distance matrix

% Initialize matrices E, F, and H
E=zeros(size(D));
F=E;
H=E;

% Row
for i=1:size(D,1)
    t=D(i,:);
    t=t==min(t);
    E(i,:)=t;
end
      
% Column
for i=1:size(D,2)
    t=D(:,i);
    t=t==min(t);
    F(:,i)=t;
end

% Thresholding (on the distances between the original positions)
H=Dorig~=0;
      
% Resulting selected distance matrix
if ~isempty(E) % If E is empty, then also F and H are empty
    G=(E & F & H); 
else
    G=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BUILD M FOR UNIVOCALLY RESOLVED PAIRS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract non-paired speckles
[gy gx]=find(G);

if isempty(gy) % And therefore also gx - no matches found
    M=[];
    return
end

m1=missingIndices(gy,size(G,1));
m2=missingIndices(gx,size(G,2));
nonPairedI=x(m1,:);      % Non paired 1->2
nonPairedJ=y(m2,:);      % Non paired 2->1

matches=[gy gx];
[r1,n1]=repeatingIndices(gy);
[r2,n2]=repeatingIndices(gx);

% More than one partner
oneMatch=cat(2,r1,r2);
conflictMatches=matches(oneMatch,:);  % In conflict
matches(oneMatch,:)=[];               % Univocally assigned

% Add univocally assigned positions to M
for c1=1:size(matches,1)
    M(c1,1:4)=[x(matches(c1,1),:) y(matches(c1,2),:)];
end

% Add non-paired positions from frame 1
lenNPI=size(nonPairedI,1);
MNPI=zeros(lenNPI,4);
MNPI(1:lenNPI,1:2)=nonPairedI;
M=cat(1,M,MNPI);

% Add non-paired positions from frame 2
lenNPJ=size(nonPairedJ,1);
MNPJ=zeros(lenNPJ,4);
MNPJ(1:lenNPJ,3:4)=nonPairedJ;
M=cat(1,M,MNPJ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BUILD MATRIX conflict, CONTAINING THE POSITIONS STILL UNRESOLVED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conflict=[];
if ~isempty(conflictMatches)
    for c1=1:size(conflictMatches,1)
        conflict(c1,1:4)=[x(conflictMatches(c1,1),:) y(conflictMatches(c1,2),:)];
    end
end

if DEBUG==1
    
    figure;
    % Plot
    plotCoords(x,y);hold on;
    for i=1:size(M,1)
        if ~(any(M(i,:)==0))
            plot([M(i,1) M(i,nDim+1)],[M(i,2) M(i,nDim+2)],'r-');
            hold on;
            plot(M(i,1),M(i,2),'k.','Markersize',20);
            plot(M(i,nDim+1),M(i,nDim+2),'r.','Markersize',10);
            axis([mn mx mn mx]);
            hold on;
        end
    end
    
    if ~isempty(conflict)
        % Add positions in conflict
        for i=1:size(conflict,1)
            plot([conflict(i,1) conflict(i,nDim+1)],[conflict(i,2) conflict(i,nDim+2)],'c-');
            hold on;
            plot(conflict(i,1),conflict(i,2),'c.','Markersize',20);
            plot(conflict(i,nDim+1),conflict(i,nDim+2),'c.','Markersize',10);
            axis([mn mx mn mx]);
            hold on;
        end
    end

    hold off;
    
end
    
            