function [] = pcPrepareForTracking(params,dirs)

cells4TrackingFname = [dirs.results dirs.expname '_cells4tracking.mat'];

if exist(cells4TrackingFname,'file') && ~params.always
    return;
end

movieInfo = [];
for t = 2 : params.nTime - params.frameJump - 1
    cellsFname = [dirs.detectData pad(t,3) '_cells.mat'];    
    
    load(cellsFname);
    ncells = length(cells.xs);
    frame.xCoord = zeros(ncells,2);
    frame.yCoord = zeros(ncells,2);
    frame.amp = zeros(ncells,2);
    
    frame.xCoord(:,2) = 0;
    frame.yCoord(:,2) = 0;
    frame.amp(:,2) = 0;
    
    for c = 1 : ncells
        y = cells.ys(c);
        x = cells.xs(c);
        frame.xCoord(c,1) = x;
        frame.yCoord(c,1) = y;
        %         frame.amp(c,1) = scores(y,x);
    end
    movieInfo = [movieInfo, frame];
end
  save(cells4TrackingFname,'movieInfo');
end
   