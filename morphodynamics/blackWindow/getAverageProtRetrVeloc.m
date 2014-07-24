function [meanVeloc,confI,cltDisp,cltConfI,index,ProtVeloc,RetrVeloc] = getAverageProtRetrVeloc(Prot,varargin)
% This function calculates the average of the protrusion/retraction instantaneous velocities
%
% Additional functionality: It can cluster regions with similar velocities 
%
%Note: Remove outliers and NaN before using this function
%      Ex: 
%      TimeSeries( detectOutliers(TimeSeries,5) ) = NaN;
%      TimeSeries = removeMeanTrendNaN(TimeSeries,'trendType',-1)
%
%Usage:
%       [meanVeloc,confI,cltDisp,cltConfI,index,ProtVeloc,RetrVeloc] = getAverageProtRetrVeloc(Prot,varargin)
%
%Input:
%       Prot - cell array with the edge motion data - e.g., protSamples.avg from
%       the protrusion sampling process - elements can have different lengths 
%
%       nBoot - # of boostrap samples to be used (default value 1000)  
%
%       alpha - alpha used to generate the bootstrap confidence intervals
%         (default value 0.05)
%
%       cluster - scalar "1" to perform cluster analysis; "0" otherwise       
%
%       nCluster - number of cluster 
%
%   
%Output:
%       meanVeloc - mean velocity of the whole data set
%                   meanVeloc(1) - for protrusion
%                   meanVeloc(2) - for retraction
%
%       confI    - confidence interval for the meanVeloc 
%                  If the data is normal, this is the standard error
%                  Otherwise, this is the percentiles of the distribution generated by the bootstrap algorithm at 
%                  the given alpha value
%
%       cltDisp  - vector containing the mean velocity for each cluster
%
%       cltConfI - matrix (2,nCluster) - confidence interval for each
%       cluster
%
%       index    - cell array with indeces for each cluster
%
%See also: getEdgeMotionPersistence
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('Prot',@(x) iscell(x) );
ip.addParamValue('nBoot',1e3,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('cluster',0,@isscalar);
ip.addParamValue('nCluster',2,@isscalar);

ip.parse(Prot,varargin{:});
nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
cluster  = ip.Results.cluster;
nCluster = ip.Results.nCluster;

meanVeloc = NaN(1,2);
confI     = NaN(2);
cltDisp   = NaN(1,2*nCluster);
cltConfI  = NaN(2,2*nCluster);
index     = cell(1,2*nCluster);


%Mean Protusion and Retraction instantaneous velocities
cellProtVeloc =  cellfun(@(x)  x( x > 0 ) ,Prot,'UniformOutput',0) ;
cellRetrVeloc =  cellfun(@(x)  x( x < 0 ) ,Prot,'UniformOutput',0) ;

ProtVeloc = cell2mat(cellProtVeloc(:));
RetrVeloc = cell2mat(cellRetrVeloc(:));

ProtVeloc(isnan(ProtVeloc)) = [];
RetrVeloc(isnan(RetrVeloc)) = [];

if adtest(ProtVeloc) %If the distribution is not Gaussian
    [confI(:,1),meanVeloc(1)] = bootStrapMean(ProtVeloc,alpha,nBoot);
    
else
    meanVeloc(1) = mean(ProtVeloc);
    quant        = norminv(1-(alpha/2),0,1);
    confI(1,1) = meanVeloc(1) + std(ProtVeloc)*quant/sqrt(length(ProtVeloc));
    confI(2,1) = meanVeloc(1) - std(ProtVeloc)*quant/sqrt(length(ProtVeloc));
end

if adtest(RetrVeloc)
    [confI(:,2),meanVeloc(2)] = bootStrapMean(RetrVeloc,alpha,nBoot);
else
    meanVeloc(2) = mean(RetrVeloc);
    quant        = norminv(1-(alpha/2),0,1);
    confI(1,2)   = meanVeloc(2) + std(RetrVeloc)*quant/sqrt(length(RetrVeloc));
    confI(2,2)   = meanVeloc(2) - std(RetrVeloc)*quant/sqrt(length(RetrVeloc));
end

if cluster
    
    [cltConfI(:,1:nCluster),cltDisp(1:nCluster),index(1:nCluster)] = ...
                                    clusterWindowsVelocity(ProtVeloc,nBoot,alpha,nCluster);
                                
    [cltConfI(:,nCluster+1:2*nCluster),cltDisp(nCluster+1:2*nCluster),index(nCluster+1:2*nCluster)] = ...
                                    clusterWindowsVelocity(RetrVeloc,nBoot,alpha,nCluster);
    
end

end%End of main function


function [conf,meanS] = bootStrapMean(variable,alpha,nBoot)
%This subFunction bootstrap the mean value of the input "variable"
% bootci calculates the confidence interval at alpha level based on the
% percentile (corrected for bias) of the distribution for the speficic
% stats
opt = statset('UseParallel','never');
if matlabpool('size')
    opt = statset('UseParallel','always');
end

[conf,meanSample] = bootci(nBoot,{@mean,variable},'alpha',alpha,...
    'type','bca','Options',opt);

meanS = mean(meanSample);

end

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