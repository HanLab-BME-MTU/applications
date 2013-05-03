function [protrusion,retraction] = mergingEdgeResults(TS,varargin)
%This function merges results from different modules
%
ip = inputParser;
ip.addRequired('TS',@(x) iscell(x) );
ip.addParamValue('nBoot',1e3,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('cluster',0,@isscalar);
ip.addParamValue('nCluster',2,@isscalar);
ip.addParamValue('deltaT',1,@isscalar);

ip.parse(TS,varargin{:});
nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
cluster  = ip.Results.cluster;
nCluster = ip.Results.nCluster;
deltaT   = ip.Results.deltaT;

%% Module 1 
[protrusion,retraction] = getEdgeMotionPersistence(TS,'cluster',cluster,'nCluster',nCluster,'alpha',alpha,'nBoot',nBoot,'deltaT',deltaT);
%[p1,r1] = getEdgeMotionPersistence(TS,'cluster',cluster,'nCluster',nCluster,'alpha',alpha,'nBoot',nBoot);
%% Module 2
%[p2,r2] = getEdgeMotionAverageVeloc(TS,'cluster',cluster,'nCluster',nCluster,'alpha',alpha,'nBoot',nBoot);
%% Next module

%%
%protrusion = MergeStruct(p1,p2);
%retraction = MergeStruct(r1,r2);

end
