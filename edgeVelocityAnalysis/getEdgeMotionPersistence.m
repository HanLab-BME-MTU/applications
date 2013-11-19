function [outP,outR] = getEdgeMotionPersistence(TS,varargin)
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
%       outP.CI.same fields as in meanValue
%
%       outP.windows     - structure with all measurements per window (see getPersistenceTime for the measurements)
%       
%
%       Same structure for outR    
%
%See also: getPersistenceTime, findingProtRetrTime
%
%Marco Vilela, 2012

%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('TS',@(x) iscell(x));
ip.addParamValue('nBoot',1e3,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('deltaT',1,@isscalar);
ip.addParamValue('cluster',false,@islogical);
ip.addParamValue('nCluster',2,@isscalar);
ip.addParamValue('winInterval',{[]},@iscell);
ip.addParamValue('selection',{[]},@iscell);

ip.parse(TS,varargin{:});

nBoot       = ip.Results.nBoot;
alpha       = ip.Results.alpha;
deltaT      = ip.Results.deltaT;
cluster     = ip.Results.cluster;
nCluster    = ip.Results.nCluster;
winInterval = ip.Results.winInterval;
selection   = ip.Results.selection;
%**************************************************************************
nWin = length(TS);

if isempty(winInterval{1})
    
    winInterval = cellfun(@(x) 1:numel(x),TS,'Unif',0);
    
end


firstP      = min(cell2mat(cellfun(@(x) x(1),winInterval,'Unif',0)));
lastP       = max(cell2mat(cellfun(@(x) x(end),winInterval,'Unif',0)));
motionState = nan(lastP-firstP+1,nWin);
newInt      = cellfun(@(x) x-firstP+1,winInterval,'Unif',0);
    


%%

outR  = setupEdgeVelocityStructure(nWin);        
outP  = setupEdgeVelocityStructure(nWin);


protrusion = struct('persTime',[],'totalTime',[],'maxVeloc',[],'Veloc',[],'minVeloc',[],'mednVeloc',[]);
retraction = struct('persTime',[],'totalTime',[],'maxVeloc',[],'Veloc',[],'minVeloc',[],'mednVeloc',[]);

for iWin = 1:nWin
    
    [outP.windows(iWin),outR.windows(iWin),motionState(newInt{iWin},iWin)] = getPersistenceTime(TS{iWin},deltaT);
    %believe it or not, this is the fastest way to do it
    protrusion.persTime  = [protrusion.persTime; outP.windows(iWin).persTime];
    protrusion.totalTime = [protrusion.totalTime;nansum(outP.windows(iWin).persTime)];
    protrusion.maxVeloc  = [protrusion.maxVeloc; outP.windows(iWin).maxVeloc];
    protrusion.Veloc     = [protrusion.Veloc;    outP.windows(iWin).Veloc];
    protrusion.minVeloc  = [protrusion.minVeloc; outP.windows(iWin).minVeloc];
    protrusion.mednVeloc = [protrusion.mednVeloc;outP.windows(iWin).mednVeloc];
    
    retraction.persTime  = [retraction.persTime; outR.windows(iWin).persTime];
    retraction.totalTime = [retraction.totalTime;nansum(outR.windows(iWin).persTime)];
    retraction.maxVeloc  = [retraction.maxVeloc; outR.windows(iWin).maxVeloc];
    retraction.Veloc     = [retraction.Veloc;    outR.windows(iWin).Veloc];
    retraction.minVeloc  = [retraction.minVeloc; outR.windows(iWin).minVeloc];
    retraction.mednVeloc = [retraction.mednVeloc;outR.windows(iWin).mednVeloc];
    
end

%*****************************************************************

%Bootstrapping the average and CI
[outP.CI,outP.meanValue] = structfun(@(x) bootStrapMean(x,alpha,nBoot),protrusion,'Unif',0);
[outR.CI,outR.meanValue] = structfun(@(x) bootStrapMean(x,alpha,nBoot),retraction,'Unif',0);

outP.total      = protrusion;
outR.total      = retraction;

outP.total.time = nansum(protrusion.persTime);
outR.total.time = nansum(retraction.persTime);

outP.total.percentage = sum( motionState > 0, 2)./sum(isfinite(motionState),2);
outR.total.percentage = sum( motionState < 0, 2)./sum(isfinite(motionState),2);

if ~isempty(selection{1})
    
    nSel = numel(selection);
    pCC = 1;
    rCC = 1;
    for iSel = 1:nSel
        out = find(eval(selection{iSel}));
        if strcmp(selection{iSel}(1:10),'protrusion')
        
            outP.selection(pCC).persTime  = protrusion.persTime(out);
            outP.selection(pCC).maxVeloc  = protrusion.maxVeloc(out);
            outP.selection(pCC).Veloc     = protrusion.Veloc(out);
            outP.selection(pCC).mednVeloc = protrusion.mednVeloc(out);
        
            pCC = pCC + 1;
        elseif strcmp(selection{iSel}(1:10),'retraction')
            
            outR.selection(rCC).persTime  = retraction.persTime(out);
            outR.selection(rCC).maxVeloc  = retraction.maxVeloc(out);
            outR.selection(rCC).Veloc     = retraction.Veloc(out);
            outR.selection(rCC).mednVeloc = retraction.mednVeloc(out);

            rCC = rCC + 1;
        else
            
            error('Selection should start with protrusion or retraction')
            
        end
        
    end
    
end


if cluster
    
    [outP.CI.cluster,outP.meanValue.cluster] = structfun(@(x) clusterWindowsVelocity(x,nBoot,alpha,nCluster),protrusion,'Unif',0);
    [outP.CI.cluster,outP.meanValue.cluster] = structfun(@(x) clusterWindowsVelocity(x,nBoot,alpha,nCluster),protrusion,'Unif',0);
    
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