function [] = pcTrackingOld(params,dirs)

trackingResultsFname = [dirs.results dirs.expname '_tracking.mat'];

fprintf(sprintf('Tracking %s\n',dirs.expname));

if exist(trackingResultsFname,'file') && ~params.always
    return;
end


cells4TrackingFname = [dirs.results dirs.expname '_cells4tracking.mat'];
load(cells4TrackingFname); % movieInfo

[tracksFinal,kalmanInfoLink,errFlag] = KhuloudTracker(movieInfo,params);

save(trackingResultsFname,'tracksFinal','kalmanInfoLink','errFlag');

end