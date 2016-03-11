function [] = pcVisualizeTracking(params,dirs)

trackingResultsVis = [dirs.results dirs.expname '_tracking.jpg'];

fprintf(sprintf('visualize tracking %s\n',dirs.expname));

if exist(trackingResultsVis,'file') && ~params.always
    return;
end


trackingResultsFname = [dirs.results dirs.expname '_tracking.mat'];

load(trackingResultsFname);%'tracksFinal','kalmanInfoLink','errFlag'

tmp = vertcat(tracksFinal.seqOfEvents);
numTimePoints = max(tmp(:,1));

h = plotTracks2D(tracksFinal,[1 numTimePoints],1);

saveas(h,trackingResultsVis,'jpg');

end