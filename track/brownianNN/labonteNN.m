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
%                            options.rf        : difference of the criterion of two consecutive 
%                                                iterations which stops the process [ default: 1e-5 ]
%                            options.alpha     : initial learning rate [ default: 5e-3 ]
%                            options.beta      : value by which .r is decreased and .alpha 
%                                                is increased at each iteration [ default: 0.99 ]
%                            options.threshold : max accepted distances between two particles in the 
%                                                final result
%                                                
% OUTPUT     M        :   matrix of matched positions [y x y x]
%            conflict :   matrix of still unresolved positions [y x y x]
%                         (with repetitions)
%            options  :   used options structure
%            v        :   x coordinates as tranformed during the process
%            w        :   y coordinates as tranformed during the process
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
%               fsmTrackMain is used by { fsmTrackTrackerA }
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
% Aaron Ponti, 2002

DEBUG=0;

% Define a constant
nDim=size(x,2);     % Number of dimensions

if nargin==2
    options.r=max([max(x(:)) max(y(:))]);
    options.rf=1e-5; 
    options.alpha=0.005;
    options.beta=0.99;
    options.threshold=1;
end

%
% Initialization
%
v=x;
w=y;

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

r=options.r;
alpha=options.alpha;

% Initialize dw and dv
dw=zeros(size(w));
dv=zeros(size(v));

crit=options.r;
lastCrit=0;
while crit>options.rf

    % Reset dw and dv
    dw=0*dw;
    dv=0*dv;
    
    %
    % FEED NETWORK 1 WITH THE WEIGHTS FROM NETWORK 2
    %
    
    % Calculate the distance of each neuron to all others
    d=createDistanceMatrix(v,w);
    
    % Now go through all neurons from network 1 and find the winning neuron from 2
    for i=1:size(d,1)
        % Extract distances to all neurons in netowrk 2 corresponding to neuron i from d
        cD=d(i,:);
        % Find the winning neuron in network 2
        vN=find(cD==min(cD)); % & cD<r);
        % The winning neuron is the one (or one of the ones) having the smallest activation
        if ~isempty(vN)
            if length(vN)>1
                vN=vN(1);
            end
            % alpha must be set to zero where the distance || v(i)-w || is > options.r
            dist=v(i,:)-w(vN,:);
            for j=1:nDim
                dwi(1:size(w,1),j)=alpha*dist(j);
            end
            currDist=createDistanceMatrix(w,w(vN,:));
            currDist=(currDist<r);
            for k=2:nDim
                currDist(:,k)=currDist(:,1);
            end
            dwi=dwi.*currDist;
            % Update weight changes
            dw=dw+dwi;
        end
        % No winner, no weight update 
    end
    
    
    %
    % FEED NETWORK 2 WITH THE WEIGHTS FROM NETWORK 1
    %

    % Calculate the distance of each neuron to all others
    d=d';
    
    % Now go through all neurons from network 2 and find the winning neuron from 1
    for j=1:size(d,1)
        % Extract distances to all neurons in network 1 corresponding to neuron j from d
        cD=d(j,:);
        % Find the winning neuron in network 1
        vN=find(cD==min(cD)); % & cD<r);
        % The winning neuron is the one (or one of the ones) having the smallest activation
        if ~isempty(vN)
            if length(vN)>1
                vN=vN(1);
            end
            % alpha must be set to zero where the distance || v(i)-w || is > options.r
            dist=w(j,:)-v(vN,:);
            for i=1:nDim
                dvj(1:size(v,1),i)=alpha*dist(i);
            end
            currDist=createDistanceMatrix(v,v(vN,:));
            currDist=(currDist<r);
            for k=2:nDim
                currDist(:,k)=currDist(:,1);
            end
            dvj=dvj.*currDist;
            % Update weight changes
            dv=dv+dvj;
        end
        % No winner, no weight update 
    end

    %
    % Update weights
    %
 
    if ~isempty(dw)
        w=w+dw;
    end
    if ~isempty(dv)
        v=v+dv;
    end
    
    %
    % Update time-dependent parameters
    %
    
    r=r*options.beta;
    alpha=alpha/options.beta;

    if DEBUG==1
        mx=max([max(x(:)) max(y(:))])+0.001;
        plot(v(:,1),v(:,2),'r*','MarkerSize',10);
        %     axis([0 mx 0 mx]);
        hold on;
        plot(w(:,1),w(:,2),'b*','MarkerSize',10);
        %     axis([0 mx 0 mx]);
        title('Convergence');
    end
    
    % Calculate distance criterion (sum of squared distance differences)
    minDist=min(createDistanceMatrix(v,w),[],2);
    sqErr=sum(minDist.^2);
    crit=abs(sqErr-lastCrit);
    lastCrit=sqErr;
    
    if DEBUG==1
        crit
        pause(0.001);
    end

end

% Build pairing matrix
D=createDistanceMatrix(v,w); % We don't set a distance constraint here
Dorig=createSparseDistanceMatrix(x,y,options.r); % Only in the original
                                                         % distance matrix

% Initialize matrices E, F, and H
E=[]; F=[]; H=[];

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
H=Dorig~=0;%<=options.threshold;
      
% Resulting selected distance matrix
if ~isempty(E) % If E is empty, then also F and H are empty
    G=(E & F & H); 
else
    G=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%
% BUILD M FOR UNIVOCALLY RESOLVED PAIRS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

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
matches(oneMatch,:)=[];        % Univocally assigned

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%
% BUILD MATRIX conflict, CONTAINING THE POSITIONS STILL UNRESOLVED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

conflict=[];
if ~isempty(conflictMatches)
    for c1=1:size(conflictMatches,1)
        conflict(c1,1:4)=[x(conflictMatches(c1,1),:) y(conflictMatches(c1,2),:)];
    end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% %
% % BUILD M FOR UNIVOCALLY RESOLVED PAIRS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% 
% % Extract non-paired speckles
% [gy gx]=find(G);
% m1=missingIndices(gy,size(G,1));
% m2=missingIndices(gx,size(G,2));
% nonPairedI=x(m1,:);      % Non paired 1->2
% nonPairedJ=y(m2,:);      % Non paired 2->1
% 
% matches=[gy gx];
% [r1,n1]=repeatingIndices(gy);
% [r2,n2]=repeatingIndices(gx);
% 
% % More than one partner
% oneMatch=cat(2,r1,r2);
% conflict=matches(oneMatch,:);  % In conflict, to be returned
% matches(oneMatch,:)=[];        % Univocally assigned
% 
% % Add univocally assigned positions to M
% for c1=1:size(matches,1)
%     M(c1,1:4)=[posI(matches(c1,1),:) posJ(matches(c1,2),:)];
% end
% 
% % Add non-paired positions from frame 1
% lenNPI=size(nonPairedI,1);
% MNPI=zeros(lenNPI,4);
% MNPI(1:lenNPI,1:2)=nonPairedI;
% M=cat(1,M,MNPI);
% 
% % Add non-paired positions from frame 2
% lenNPJ=size(nonPairedJ,1);
% MNPJ=zeros(lenNPJ,4);
% MNPJ(1:lenNPJ,3:4)=nonPairedJ;
% M=cat(1,M,MNPJ);
% 
% 
% 
% % % Assign pairings
% % pairsI=[]; c0=0; toBeAdded=[]; cCat=1;
% % for c1=1:size(G,1)
% %     pos=find(G(c1,:)); % Find matches 
% %     if length(pos)==1
% %         c0=c0+1;
% %         pairsI(c0,1:4)=[x(c1,:) y(pos,:)];
% %     elseif length(pos)==0
% %         c0=c0+1;
% %         pairsI(c0,1:4)=[x(c1,:) 0 0];
% %     else
% %         % More than one match found - this is a rare event
% %         % Assign the one which was originarily closer
% %         % and discard the other(s)
% %         %
% %         % Add first match
% %         d=createDistanceMatrix(x(c1,:),y(pos,:));
% %         posCloser=find(d==min(d));
% %         if length(posCloser)>1      % If they were at the same distance 
% %             posCloser=posCloser(1); % also in the beginning, just take one     
% %         end
% %         pairsI(c1,1:4)=[x(c1,:) y(pos(posCloser),:)];
% %         % Remove assigned match
% %         pos(posCloser)=[];
% %         % Add the rest
% %         for c2=0:length(pos)-1
% %             toBeAdded(cCat+c2,1:4)=[0 0 y(pos(c2+1),:)];
% %         end
% %     end
% % end
% % % Check for non-paired speckles in frame 2
% % pairsJ=[]; c0=0;
% % for c1=1:size(G,2)
% %     pos=find(G(:,c1)); % Find matches
% %     if length(pos)==0
% %         c0=c0+1;
% %         pairsJ(c0,1:4)=[0 0 y(c1,:)];
% %     end
% %     % Other cases have already been treated in the check 
% %     % for speckles in frame 1
% % %     if length(pos)>1
% % %         disp('One speckle in frame 2 with two candidates in 1');
% % %     end
% % end
% % % Concatenate
% % M=cat(1,pairsI,pairsJ);
% % M=cat(1,M,toBeAdded);

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
    
            