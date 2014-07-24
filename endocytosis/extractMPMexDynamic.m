function [extract]=extractMPMexDynamic(data, restvector)
% for the hotspot analysis, extract MPM of points with certain
% characteristics from the data
% input     data field
%           udist  = distance vector 
%           resvector = restriction vector, can have variable length:
%           [stat da minfr minlft maxlft minint maxint minmot maxmot]
% output    extract(i).mpm  = mpm of points
%
% NOTE: This function outputs the trajectories of all suitable points (e.g.
% objects of lifetime > 300s), thus the dimensions of the extracted data is
% (n x nf*2) - this data contains the dynamic time information which can be
% used for motion correlation analysis
%
% last updated: 02/22/2008 by Dinah Loerke
% Francois Aguet, last modified 02/09/2010

for i = 1:length(data)

    if isfield(data,'lftInfo') && ~isempty(data(i).lftInfo)
        currLI = data(i).lftInfo;
    else
        trackInfoPath = [data(i).source 'TrackInfoMatrices' filesep 'trackInfo.mat'];
        if (exist(trackInfoPath, 'file')==2)
            tfile = load(trackInfoPath);
            trackFields = fieldnames(tfile);
            currTrackInfo = tfile.(trackFields{1});
                
            % calculate lft matrices
            [CMlft,CMstat,CMx,CMy,CMda] = findLifetimesStatusSimple(currTrackInfo);
            currLI.Mat_lifetime = CMlft;
            currLI.Mat_status = CMstat;
            currLI.Mat_xcoord = CMx;
            currLI.Mat_ycoord = CMy;
            currLI.Mat_disapp = CMda;
              
        else
            currLI = [];
        end
    end
    
    % status matrix
    mat_stat = full(currLI.Mat_status);
    % frame number
    [nx,nf] = size(mat_stat);
    % status vector
    vec_stat = zeros(nx,1);
    for k=1:nx, vec_stat(k) = min(nonzeros(mat_stat(k,:))); end
    
    % lifetime matrix
    mat_lft = full(currLI.Mat_lifetime);
    % lifetime vector
    vec_lft = zeros(nx,1);
    for k=1:nx, vec_lft(k) = max(mat_lft(k,:)); end
    
    % x-coordinate matrix
    mat_x = currLI.Mat_xcoord;
    % y-coordinate matrix
    mat_y = currLI.Mat_ycoord;
    % framerate
    fr = data(i).framerate;
    % image size
    imsiz = data(i).imagesize;
    % frame number
    [nx,nf] = size(mat_lft, 2);
    
%     % intensity - rescale values to range between 0 and 1, where the value
%     % corresponds to the percentile of all available intensities
%     svec_int = [];
%     svec_int(:,1) = currLI.Vec_int; le = length(currLI.Vec_int);
%     svec_int(:,2) = 1:le;
%     svec_int = sortrows(svec_int,1);
%     svec_int(:,3) = [1:le]/le;
%     svec_int = sortrows(svec_int,2);
%     % mat_int = []; mat_int = repmat(svec_int(:,3),1,nf);
%     vec_int = svec_int(:,3);
%     
%     % motility - rescale values to range between 0 and 1
%     svec_mot = [];
%     svec_mot(:,1) = currLI.Vec_displ; le = length(currLI.Vec_displ);
%     svec_mot(:,2) = 1:le;
%     svec_mot = sortrows(svec_mot,1);
%     svec_mot(:,3) = [1:le]/le;
%     svec_mot = sortrows(svec_mot,2);
%     % mat_mot = []; mat_mot = repmat(svec_mot(:,3),1,nf);
%     vec_mot = svec_mot(:,3);
    
    %=====================================================================
    % DESIRED CONDITIONS
    % for all frames in the movie, collect those locations of points that
    % fulfill a number of requirements
    % e.g. status, minimum lifetime, maximum motility
    
    % desired minimum lifetime in seconds
    dstat       = restvector(1);
%   dapp        = restvector(2);
    dminfr      = restvector(3);
    minlft      = restvector(4);
    maxlft      = restvector(5);
    minlft_fr   = round(minlft/fr);
    maxlft_fr   = round(maxlft/fr);
    
%     if length(restvector)>5
%         minint = restvector(6);
%         maxint = restvector(7);
%         minmot = restvector(8);
%         maxmot = restvector(9);
%     else
%         minint = 0;
%         maxint = 1;
%         minmot = 0;
%         maxmot = 1;
%     end
    
    
    %=====================================================================
    
    % required conditions
%     findpos = find( ( vec_stat==dstat ) &...
%                     ( vec_lft>dminfr ) &...
%                     ( vec_lft>minlft_fr ) &...
%                     ( vec_lft<maxlft_fr ) &...
%                     ( vec_int>=minint ) &...
%                     ( vec_int<=maxint ) &...
%                     ( vec_mot>=minmot ) &...
%                     ( vec_mot<=maxmot ) );
    
    % required conditions
    findpos = find( ( vec_stat==dstat ) &...
                    ( vec_lft>dminfr ) &...
                    ( vec_lft>minlft_fr ) &...
                    ( vec_lft<maxlft_fr ) );
    
    findx = full(mat_x(findpos,:));
    findy = full(mat_y(findpos,:));
    
    currMPM = NaN(length(findpos),2*nf);
    currMPM(:,1:2:2*nf) = findx;
    currMPM(:,2:2:2*nf) = findy;
    
    extract(i).mpm = currMPM;
    extract(i).imagesize = [imsiz(2) imsiz(1)];    
end    