function [result] = MUEC_intensity(data, restvector, dvector, tvector, channel); 
% MUEC (Multichannel Universal Event Correlator) Intensity: correlates
% tracking events with intensity
% 
% SYNOPSIS [result] = MUEC_intensity(data, restvector, dvector, tvector, channel); 
%
% INPUT:    data:       data structure
%           restvector: restriction vector
%                   NOTE1: if restvector is a string instead of a vector, 
%                   then the function will try to read the vector from the 
%                   structure field of the corresponding name
%                   NOTE2: standard format for the restriction vector is
%                   [lstat dstat minfr minsec maxsec], where
%                   lstat = lifetime status [nan,1,2,3]
%                   dstat = appearance/disappearance status [-1,0,1]
%                   minfr = minimum lifetime in frames
%                   minsec = minimum lifetime in seconds
%                   maxsec = maximum lifetime in seconds
%                   EXAMPLE: [nan -1 4 60 800] extracts the disappearance
%                   events of objects that live >4 frames and between
%                   60-800 seconds (determined from .framerate), lifetime
%                   status is irrelevant
%           dvector:    distance vector in pixels (e.g. 2 for 2 pix radius)       
%           tvector:    time shift vector (e.g. [-5:10])
%           channel (optional): structure field in data from which to read
%               the intensity data; if no value is specified, the function
%               will try to find image data at the path location specified
%               in the structure field '.channel2'
%
% OUTPUT:   result:  contains the fields:
%                   .parfield = averaged parameter results above background
%                   .parerror = standard error of the mean (all events)
%                   .restrictions = restriction vector (see above)
%                   .distvector = distance vector (see above)
%                   .timevector = tvector (see above)
%                   .numtraj = total number of trajectories      
%
% Dinah Loerke
% last modified:    August 28, 2008
%                   January 20, 2009
%
% NOTE: the assumption is that the entered tvector (e.g. [-5:10]) is in
% frames and corresponds to the fastest framerate present in the analyzed 
% data (fframerate); data with slower framerates are also collected for the
% same number of frames pre- and post-event, but in the final averaging
% analysis, the data are interpolated only for the time range of the
% fastest rate


od = cd;
nr = 10;
refplane = find(tvector==0);

%% select image files from which intensity is to be read; if no specific
% location is specified, the function will try to read the image locations
% from the structure field 'channel2'

for i=1:length(data)
    
    % select first image for the second (intensity) channel
    if nargin>4
        if isfield(data, channel)
            path2 = getfield(data,{1,i},channel);
        else
            error(['no field of the specified name (,'channel,') is found');
        end
    else
        if isfield(data, channel2)
            path2 = data(i).channel2;
        else
            error('function requires the existence of a field called .channel2');
        end
    end
    
    % change to specified directory
    cd(path2); 
    % select first image file
    [imageName, imagePath] = uigetfile('.tif',['Select first original intensity image for movie #',num2str(i)]); 
    % make a list of all complete image file names (including path) for
    % later upload of the images
    completeImageName = strcat(imagePath, imageName);
    imageStackList = getFileStackNames(completeImageName);
    
    Ch2stackList(i).list = imageStackList;
    
end
pause(0.1);


%% loop over all movies
for i=1:length(data)
    
    fprintf('movie #%02d',i);
    fprintf('\n');
    
    % load trackinfo file
    cpath = [data(i).source,filesep,'TrackInfoMatrices'];
    trackinfopath = [data(i).source,filesep,'LifetimeInfo'];
    cd(cpath);
    loadfile = open('trackInfo.mat');
    trackInfo  = loadfile.trackInfo;
    framerate = data(i).framerate;
    
    % load reference image - this is a binary mask that shows the outline
    % of the cell; the mask is required because the function will choose
    % random positions inside the cell as reference positions for the
    % background signal, and it is desirable to choose the reference
    % positions within the area of the cell (as opposed to the empty area 
    % on the coverslip next to it)
    cpath = [data(i).source,filesep,'SubregionsMask'];
    cd(cpath);
    image = imread('mask0001.tif');
    mask = image';
    
    % if necessary, read restriction vector from the specified structure
    % field
    if isstr(restvector)
        restrictions = getfield( data(i),restvector );
    else
        % extract positions with desired restricition properties
       restrictions = restvector;
    end

    % extract the positions with the desired properties
    MPMglobal = CorrelateData2Pos_extractMPM2L(trackinfopath, restrictions, tvector, framerate, 1);
    
    % if the intensity belongs to a second color channel that is shifted
    % with respect to the tracking channel - i.e. if a field 
    % colorShiftVector exists - then the event positions in
    % MPMglobal have to be shifted by the shift vector
    % NOTE: What do the dimensions of the shiftvector mean? For example, a 
    % shift of shiftvec=[-10,-5]) means that image2 is shifted in such a 
    % way that the point im1(1,1) in image1 (e.g. clathrin) overlays point 
    % im2(11,6) in image2 (e.g. actin). Thus, the shift has to be subtracted
    % from the positions in the clathrin channel to obtain the correct
    % coordinates in the actin channel. Positions outside the shifted image
    % dimensions have to be set to nan
    
    if isfield(data,'colorShiftVector')
        if ~isempty(data(i).colorShiftVector)
            [msx,msy,msz] = size(MPMglobal);
            [isx,isy] = size(image);
            shiftx = data(i).colorShiftVector(1);
            shifty = data(i).colorShiftVector(2);
            MPMshiftx = MPMglobal(:,1:2:msy)-shifty;
            MPMshifty = MPMglobal(:,2:2:msy)-shiftx;
                       
            badpos = find( (MPMshiftx>isy) | (MPMshiftx<1) | (MPMshifty>isx) | (MPMshifty<1));
            MPMshiftx(badpos) = nan;
            MPMshifty(badpos) = nan;
                  
            MPMglobal(:,1:2:msy,1) = MPMshiftx;
            MPMglobal(:,2:2:msy,1) = MPMshifty;
        end
    end
            
                
    % compact MPMglobal (get rid of nan rows)
    projpos_global = nanmean(MPMglobal(:,:,1),2);
    fpos_global = find(isfinite(projpos_global));
    MPMglobal = MPMglobal(fpos_global,:,:);
        
    % make random positions for reference measurement; every event in
    % MPMglobal will be supplemented by 10 points randomly chosen in the
    % cell area (but in the same time frame as the event)
    MPMbas = MPMglobal(:,:,1);
    goodpos = find( MPMglobal(:,:,2) == refplane );
    MPMbas_use = nan*MPMbas;
    MPMbas_use(goodpos) = MPMbas(goodpos);
    MPMbas = [];
    
    for k=1:nr
        MPMrandom_curr = makeRandomMPM(MPMbas_use, mask, 1);
        if k==1
            [srx,sry] = size(MPMrandom_curr);
            MPMrandom_all = zeros(10*srx,sry);
        end
        MPMrandom_all((k-1)*srx+1:k*srx,:) = MPMrandom_curr;
    end
       
    % time-shift random positions
    MPMrandomshift = CorrelateData2Pos_timeshiftMPM2L(MPMrandom_all, tvector);
    MPMrandomshift(MPMrandomshift==0) = nan;
            
    [lx1,ly1,lz1] = size(MPMglobal);
    [lx2,ly2,lz2] = size(MPMrandomshift);
    
    MPMall = [MPMglobal ; MPMrandomshift];    
        
    imageStackList = Ch2stackList(i).list;    
    total_frame_num = length(imageStackList);
    
    % read all images for this movies
    for k=1:total_frame_num
        
        cimage = imread(imageStackList{k});
        if k==1
            intensityImageStack = zeros(size(cimage,1),size(cimage,2),total_frame_num);
        end
        intensityImageStack(:,:,k) = cimage;
        
    end
    
    % perform the correlation
    CORRresults_all = CorrelateData2Pos(MPMall, intensityImageStack, dvector, 'intensity');
    
    % the results contain a matrix of intensity values for all points     
    ResultsGlobalData(i).intmatrix = CORRresults_all(1:lx1,:,:);
    ResultsRandomData(i).intmatrix = CORRresults_all(lx1+1:lx1+lx2,:,:);
    ResultsNumTraj(i) = lx1;
        
end % of for i-loop

% return to original directory
cd(od);


%% the raw data is now processed for averaging

% NOTE: the function will average the results taking into account possibly 
% variable framerates; if necessary, the slower movie results are 
% interpolated to the fastest available framerate
for i=1:length(data)
    cfr = data(i).framerate;
    if i==1
        ffr = cfr;
    else
        ffr = min(ffr,cfr);
    end
    framerates(i) = cfr;
end
% ffr ist now fastest framerate in data set


figure; hold on;

for i=1:length(data)
    
    % event results
    ResultsGlobalCurr       = ResultsGlobalData(i).intmatrix;
    % results of the random reference points
    ResultsRandomCurr       = ResultsRandomData(i).intmatrix;
    % results of the random reference points are averaged and written into
    % a matrix matching ResultsGlobalCurr for subtraction
    ResultsRandomCurrAV     = nanmean(ResultsRandomCurr,1);
    ResultsRandomCurrAVMAT  = repmat(ResultsRandomCurrAV,ResultsNumTraj(i),1);
    ResultsGlobalCurrCORR   = ResultsGlobalCurr - ResultsRandomCurrAVMAT;
    
    % interpolate for fastest framerate is necessary
    if framerates(i)>ffr
        tff = ffr * tvector;
        tcf = framerates(i) * tvector;
        
        ResultsGlobalCurr_EP = interp1(tcf,ResultsGlobalCurr',tff);
        ResultsRandomCurrAV_EP = interp1(tcf,ResultsRandomCurrAV,tff);
        ResultsRandomCurrAVMAT  = repmat(ResultsRandomCurrAV_EP,ResultsNumTraj(i),1);
        ResultsGlobalCurrCORR   = ResultsGlobalCurr_EP' - ResultsRandomCurrAVMAT;
    end      
 
    % add the background-subtracted results to a global results matrix    
    if i==1
        ResultsGlobalAll = ResultsGlobalCurrCORR;
    else
        ResultsGlobalAll = [ResultsGlobalAll;ResultsGlobalCurrCORR];
    end
    
    [rx,ry,rz] = size(ResultsRandomCurrAV);
    ResultsGlobalCurrAV = nanmean(ResultsGlobalCurr,1);
    
    if rz>1
        for r=1:rz
            plot(scvec*tvector,ResultsRandomCurrAV(:,:,r),'r-');
            plot(scvec*tvector,ResultsGlobalCurrAV(:,:,r),'b-');
        end
    else
        plot(scvec*tvector,ResultsRandomCurrAV,'r-');
        plot(scvec*tvector,nanmean(ResultsGlobalCurr,1),'b-');
    end
    
end



% average the results of individual movies
ResultsGlobalAVE = nanmean(ResultsGlobalAll,1);
ResultsGlobalError = nanstd(ResultsGlobalAll,[],1)/sqrt(sum(ResultsNumTraj));

% write results
result.parfield = ResultsGlobalAVE;
result.parerror = ResultsGlobalError;

result.restrictions = restvector;
result.distvector = dvector;
result.timevector = tvector;
result.numtraj = sum(ResultsNumTraj);


% display results
figure
[rx,ry,rz] = size(ResultsGlobalAVE);
if rz==1
    errorbar(scvec*tvector,ResultsGlobalAVE,ResultsGlobalError);
else
    for r=1:rz
        errorbar(scvec*tvector,ResultsGlobalAVE(:,:,r),ResultsGlobalError(:,:,r));
    end
end

amin = min(scvec*tvector)-0.5;
amax = max(scvec*tvector)+0.5;
axis([amin amax -2.5e-7 2.5e-7]);
xlabel('time point relative to CCP event (sec)');
ylabel('relative kinetic score');



end % of function



    
    