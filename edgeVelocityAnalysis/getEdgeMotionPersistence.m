function [Protrusion,Retraction] = getEdgeMotionPersistence(TS,varargin)
% This function bootstrap the mean and confidence intervals for the
% protrusion and retraction persistence time using all time series in input TS
%THIS FUNCTION HAS NO TIME SERIES PRE-PROCESSING
%
%Usage: [Protrusion,Retraction] = getEdgeMotionPersistence(TS,'nBoot',1000,'alpha',0.05)
%
%Input:
%       TS - cell array with one time series in each element
%       Optional:
%               nBoot    - scalar  - number of bootstrap samples
%               alpha    - scalar  - confidence level
%               deltaT   - frame rate
%               cluster  - logical - true for cluster the data (fuzzy k-means based on the edge velocity magnitude)
%               nCluster - scalar  - number of clusters
%Output:
%
%       Protrusion.meanValue.persTime - persistence time
%                           .maxVeloc
%                           .Veloc
%                           .minVeloc
%                           .mednVeloc - median velocity (not sure what for)
%                           .cluster.persTime - persistence time
%                                   .maxVeloc
%                                   .Veloc
%                                   .minVeloc
%                                   .mednVeloc - median velocity (not sure what for)
%
%       Protrusion.CI.same fields as in meanValue
%
%       Protrusion.windows     - structure with all measurements per window (see getPersistenceTime for the measurements)
%       
%
%       Same structure for Retraction    
%
%See also: getPersistenceTime, findingProtRetrTime
%
%Marco Vilela, 2012

%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('TS',@(x) iscell(x));
ip.addOptional('nBoot',1e3,@isscalar);
ip.addOptional('alpha',.05,@isscalar);
ip.addOptional('deltaT',1,@isscalar);
ip.addOptional('cluster',false,@islogical);
ip.addOptional('nCluster',2,@isscalar);
ip.parse(TS,varargin{:});

nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
deltaT   = ip.Results.deltaT;
cluster  = ip.Results.cluster;
nCluster = ip.Results.nCluster;

%**************************************************************************


%%
nWin       = length(TS);
Retraction = setupEdgeVelocityStructure(nWin);        
Protrusion = setupEdgeVelocityStructure(nWin);

auxP = struct('persTime',[],'maxVeloc',[],'Veloc',[],'minVeloc',[],'mednVeloc',[]);
auxR = struct('persTime',[],'maxVeloc',[],'Veloc',[],'minVeloc',[],'mednVeloc',[]);

for iWin = 1:nWin
    
    [Protrusion.windows(iWin),Retraction.windows(iWin)] = getPersistenceTime(TS{iWin},deltaT);
    %believe it or not, this is the fastest way to do it
    auxP.persTime  = [auxP.persTime;Protrusion.windows(iWin).persTime];
    auxP.maxVeloc  = [auxP.maxVeloc;Protrusion.windows(iWin).maxVeloc];
    auxP.Veloc     = [auxP.Veloc;Protrusion.windows(iWin).Veloc];
    auxP.minVeloc  = [auxP.minVeloc;Protrusion.windows(iWin).minVeloc];
    auxP.mednVeloc = [auxP.mednVeloc;Protrusion.windows(iWin).mednVeloc];
    
    auxR.persTime  = [auxR.persTime;Protrusion.windows(iWin).persTime];
    auxR.maxVeloc  = [auxR.maxVeloc;Protrusion.windows(iWin).maxVeloc];
    auxR.Veloc     = [auxR.Veloc;Protrusion.windows(iWin).Veloc];
    auxR.minVeloc  = [auxR.minVeloc;Protrusion.windows(iWin).minVeloc];
    auxR.mednVeloc = [auxR.mednVeloc;Protrusion.windows(iWin).mednVeloc];
    
end

%*****************************************************************

%Bootstrapping the average and CI
[Protrusion.CI,Protrusion.meanValue] = structfun(@(x) bootStrapMean(x,alpha,nBoot),auxP,'Unif',0);
[Retraction.CI,Retraction.meanValue] = structfun(@(x) bootStrapMean(x,alpha,nBoot),auxR,'Unif',0);

if cluster
    
    [Protrusion.CI.cluster,Protrusion.meanValue.cluster] = structfun(@(x) clusterWindowsVelocity(x,nBoot,alpha,nCluster),auxP,'Unif',0);
    [Protrusion.CI.cluster,Protrusion.meanValue.cluster] = structfun(@(x) clusterWindowsVelocity(x,nBoot,alpha,nCluster),auxP,'Unif',0);
    
end

end%End of main function

function [cltConfI,cltDisp,index] = clusterWindowsVelocity(Veloc,nBoot,alpha,nCluster)

    cltDisp  = zeros(1,nCluster);
    cltConfI = zeros(2,nCluster);
    index    = cell(1,nCluster);
    [~,U]    = fcm(Veloc, nCluster);
    maxU     = max(U);
    
    for i =1 : nCluster
        index{i} = find(U(i,:) == maxU);
        [cltConfI(:,i),cltDisp(i)] = bootStrapMean( Veloc(index{i}),alpha,nBoot );
    end
    [cltDisp,idx] = sort(cltDisp); 
    cltConfI      = cltConfI(:,idx); 
end    

function outStruc = setupEdgeVelocityStructure(nWin)

outStruc.meanValue        = struct('persTime',[],'maxVeloc',[],'veloc',[],'minVeloc',[],'mednVeloc',[],...
                                'persTimeCluster',[],'maxVelocCluster',[],'velocCluster',[],'minVelocCluster',[],'mednVelocCluster',[]);

outStruc.CI               = struct('persTime',[],'maxVeloc',[],'veloc',[],'minVeloc',[],'mednVeloc',[],...
                                'persTimeCluster',[],'maxVelocCluster',[],'velocCluster',[],'minVelocCluster',[],'mednVelocCluster',[]);
                            
outStruc.windows(1:nWin)  = struct('limit',[],'persTime',[],'blockOut',[],'maxVeloc',[],'Veloc',[],'minVeloc',[],'mednVeloc',[]);

end