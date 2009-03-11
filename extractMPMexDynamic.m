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


% loop over all movies
%h = waitbar(0,['looping over ',num2str(length(data)),' movies...']);

odir = cd;


% loop over all fields in the data
for i=1:length(data)
    
    %=====================================================================
    % CURRENT lftInfo
    % determine if lftInfo is already available in the data structure of if
    % it has to be read/calculated from trackInfo
    varReadLftInfo = 1;
    if isfield(data,'lftInfo')
        if ~isempty(data(i).lftInfo)
            varReadLftInfo = 0;
        end
    end
    
    % read or load depending on the variable
    if varReadLftInfo==0
        currLI = data(i).lftInfo;
    else
        intDir = cd;    
        % go to appropriate directory and load trackInfo file
        currPath = data(i).source;
        cd(currPath);

        % if TrackInfoMatrices folder exists
        if exist('TrackInfoMatrices')==7

            cd('TrackInfoMatrices');
            % if trackInfo data exists
            if exist('trackInfo.mat')==2

                % load trackInfo file
                loadfile = load('trackInfo.mat');
                if isfield(loadfile,'trackInfo')
                    currTrackInfo = loadfile.trackInfo;
                elseif isfield(loadfile,'trackInfoMat')
                    currTrackInfo = loadfile.trackInfoMat;
                else
                    currTrackInfo = [];
                end
                
                % calculate lft matrices
                [CMlft,CMstat,CMx,CMy,CMda,CVlft]=findLifetimesStatusSimple(currTrackInfo);
                currLI.Mat_lifetime = CMlft;
                currLI.Mat_status = CMstat;
                currLI.Mat_xcoord = CMx;
                currLI.Mat_ycoord = CMy;
                currLI.Mat_disapp = CMda;
              
            else
                currLI = [];
            end % of if trackInfo exists
        else
            currLI = [];
        end % of if TrackInfo folder exists
        
        cd(intDir);
        
    end % of if varReadlftInfo == 0
    
    
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
    % disapp status matrix
    mat_da = currLI.Mat_disapp;
    % framerate
    fr = data(i).framerate;
    % image size
    imsiz = data(i).imagesize;
    msx = imsiz(1);
    msy = imsiz(2);
    % frame number
    [nx,nf] = size(mat_lft);
    
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
    dapp        = restvector(2);
    dminfr      = restvector(3);
    minlft      = restvector(4);
    maxlft      = restvector(5);
    minlft_fr   = round(minlft/fr);
    maxlft_fr   = round(maxlft/fr);
    
    if length(restvector)>5
        minint = restvector(6);
        maxint = restvector(7);
        minmot = restvector(8);
        maxmot = restvector(9);
    else
        minint = 0;
        maxint = 1;
        minmot = 0;
        maxmot = 1;
    end
    
    
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
    
    currMPM = nan*zeros(length(findpos),2*nf);
    currMPM(:,1:2:2*nf) = findx;
    currMPM(:,2:2:2*nf) = findy;
    
    extract(i).mpm = currMPM;
    extract(i).imagesize = [imsiz(2) imsiz(1)];
    
    %waitbar(i/length(data));
    
end

%close(h);

cd(odir);

end % of function
       