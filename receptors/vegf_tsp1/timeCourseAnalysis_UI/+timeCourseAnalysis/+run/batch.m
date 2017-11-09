function [ job ] = batch( p, cluster , maxWorkers)
%timeCourseAnalysis.run.batch Run timeCourseAnalysis in batch mode

    if(nargin < 2 || isempty(cluster))
        cluster = parcluster(p.batchClusterName);
    end
    if(nargin < 3)
        maxWorkers = 24;
    end
    job = batch(cluster,@timeCourseAnalysis, ...
        0, ... % No output
        {p.CML_FullPath, p.outputDir, p}, ...
        'AutoAttachFiles',false, ...
        'CaptureDiary',true, ...
        'Pool',min(cluster.NumWorkers,maxWorkers) - 1 ...
    );

end

