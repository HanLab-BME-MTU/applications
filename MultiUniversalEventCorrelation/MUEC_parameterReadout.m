function [intRes] = MUEC_parameterReadout(data, restvector, channel, reference, sframerate, align, voidAlignPoint, parname)
% Multichannel Universal Event Correlator (MUEC) to read out image
% intensities or parameter values at specified object locations
%
% SYNOPSIS  [intRes] = MUEC_parameterReadout(data, restvector, channel,...
%                           reference, sframerate, align, voidAlignPoint, parname);
%
% INPUT:    data:       data structure, containinng the fields
%                       .source (which points to the data directory
%                       containing lifetime data)
%                       .framerate
%           restvector: restriction vector to define objects of interest;
%                       this vector traditionally has the format
%                       [stat da minfr minsec maxsec]
%                       Definitions see further below
%           channel:    OPTIONAL - designates the path where to look for
%                       parameter data; DEFAULT = directory in .source
%           reference:  OPTIONAL - determines whether or not to upload any
%                       reference data (e.g. background reference
%                       intensities); DEFAULT = 0 ('no')
%           sframerate: OPTIONAL - standard framerate to which all other
%                       framerates are interpolated; DEFAULT = 2sec
%           align:      OPTIONAL - indicate how data are aligned:
%                       'left' (appearance)
%                       'right' (disappearance)
%                       'edge' (right AND left with sigmoid-weighted
%                       superposition);
%                       DEFAULT = 'egde'
%           voidAlignPoint:      OPTIONAL  - indicates whether or not (1/0) to void
%                       left and right alignment point, as in the case of
%                       intensity, where this point has an artifact;
%                       DEFAULT  = 0 (no voiding)
%           parname:    OPTIONAL - if you know the parameter name, you can
%                       enter the name stump, e.g. 'ClathrinInt', and the
%                       function will look for a parameter of that
%                       name - otherwise, the user is asked to specifiy the
%                       parameter file in a GUI
%
% OUTPUT:   intRes:     intensity Results file with the fields
%                       .intMat     = intensity Matrix (number of rows
%                                     equals number of movies)
%                       .errorMat   = error Matrix
%                       .sn         = number of individual trajectories in
%                                     each movie for this condition
%                       .tvec       = time vector;
%                       .framerate  = reference framerate;
%                       .intAVE     = average intensity
%                       .intSEM     = standard error of the mean (for all
%                                     combined data)
%
%
%   DEFINITIONS for restriction Vector:
%
%   stat =  trajectory status [1,2,3]
%           (1: trajectory appears and disappears, 2: trajectory is cut off
%           at the beginning OR end, 3: trajectory is present for entire
%           movie, i.e. cut off at beginning AND end)
%   da =    disappearance status [-1,0,1]
%           (-1: last frame of trajectory, +1, first frame of trajectory,
%           0 all other frames)
%   minfr =  minimum lft in frames, e.g. 4
%   minsec = minimum lifetime in seconds
%   maxsec = maximum lifetime in seconds
%
%   NOTE: the value for da is not really used in any way in this function,
%   since the parameter/intensity is read out for the entire object
%   trajectory anyway, but for historic reasons we keep the restriction
%   vector the way it is for now
%
% Dinah Loerke, March 18, 2009
% Francois Aguet, Feb 2010

if nargin<4 || isempty(reference)
    reference = 0;
end
if nargin<5 || isempty(sframerate)
    sframerate = 2;
end

if nargin>=6
    if strcmp(align,'left')
        alignvar = 1;
    elseif strcmp(align,'right')
        alignvar = 2;
    else
        alignvar = 3;
    end
else
    alignvar = 3;
end
if nargin<7 || isempty(voidAlignPoint)
    voidAlignPoint = 0;
end
if nargin>7 && ischar(parname)
    paramName = [parname '.mat'];
    paramNameRef = [parname '_Ref.mat'];
else
    paramName = 'parameterMat.mat';
    paramNameRef = 'parameterMat_Ref.mat';
end

% ========================================================================
%  load parameter files for all movies
% =========================================================================
nMovies = length(data);
fileName = cell(1:nMovies);
fileNameRef = cell(1:nMovies);

for i = 1:nMovies
    % set search path for parameter/intensity file
    if nargin>2 && ~isempty(channel)
        sourcePath = data(i).(channel);
    else
        sourcePath = data(i).source;
    end
        
    if (exist(sourcePath, 'dir')==7)        
        fileName{i} = getParameterFileName(sourcePath, paramName);
        if reference==1
            fileNameRef{i} = getParameterFileName(sourcePath, paramNameRef);
        end
        % if path is file name, use this file
    elseif (exist(sourcePath, 'file')==2)
        fileName{i} = sourcePath;
    else
        error('no parameter file found');
    end
end

% ========================================================================
%  extract parameter information from files for all movies
% =========================================================================

% trajectory min/max length, standard time vector
padding = 15;
tlen_min = restvector(4)/sframerate;
tlen_max = restvector(5)/sframerate;
tlen_st  = floor((tlen_min+tlen_max)/2);
tvec_standard(1,:) = -padding:tlen_st+padding-1;
tvec_standard(2,:) = -(tlen_st+padding-1):padding;


% loop over all movies
figure;
for i = 1:nMovies
    
    fprintf('Movie no. %d\n', i);

    % extract the positions with the desired parameters (e.g. lifetime range) from lft info data
    posvec = extractPos_vecFromLftRest([data(i).source 'LifetimeInfo'], restvector, data(i).framerate);
    
    % load parameter files
    paramMat = load(fileName{i});
    paramMat = paramMat.iMat_obj;
    
    paraInt = paramMat(posvec,:,1);
    paraT1  = paramMat(posvec,:,2);
    paraT2  = paramMat(posvec,:,3);
    px = size(paraInt, 1);
    
    %%===================================================================
    % if reference parameter values are used, load and calculate them here
    if reference==1
        paramMatRef = load(fileNameRef{i});
        paramMatRef = paramMatRef.iMat_ref;
        
        % determine how many reference points exist for each data point
        % (e.g. 8 or 10 background points per data point)
        rfac = size(paramMatRef,1)/size(paramMat,1);
        % if rfac>1 (e.g. rfac=8), then the assumption is that the 8 ref
        % points for data point number one are in the rows 1-8, those for
        % data point number 2 in rows 9-16, etc.
        % for easier use in the following function, the reference data are
        % immediately averaged here, i.e. they are compressed from rfac>1
        % to a matrix the same size as paramMat_use(:,:,1)
        paraInt2 = NaN*paraInt;
        if rfac>1
            for k=1:length(posvec)
                % current data point
                cpos = posvec(k);
                % all reference points for this data point
                rpos = (rfac*(cpos-1)+1):(rfac*cpos);
                % parameter values for these reference points
                rpar = paramMatRef(rpos,:,1);
                % average
                paraInt2(k,:) = nanmean(rpar);
            end
        else
            paraInt2 = paramMatRef(posvec,:,1);
        end
    end
    
    
    % ===================================================================
    % data alignment:
    % "chaotic evil means never having to say you're sorry"
    %
    % the data are a conglomerate of intensity trajectories of different
    % lengths; they are averaged by first aligning them to the appearance time
    % point, then to the disappearance time point, and then forming a weighted
    % average - this procedure conserves the sharp increase in e.g. intensity
    % during the nucleation and internalization phase
    
    
    
    % ==========================
    % APPEARANCE-aligned
    tvec_app = tvec_standard(1,:);
    emat_app = nan*zeros(px,length(tvec_app));
    if reference==1, emat_appRef = emat_app; end
    for p=1:px
        tloc = paraT1(p,:);
        for t=1:length(tvec_app)
            ctau = tvec_app(t);
            tfind = find(tloc==ctau);
            if ~isempty(tfind)
                emat_app(p,t) = paraInt(p,tfind);
                % if reference parameter values exist, read these and
                % later subtract from data parameter value
                if reference==1
                    emat_appRef(p,t) = paraInt2(p,tfind);
                end
            end
        end
    end
    % if reference positions exist, subtract their signal here
    if reference==1, emat_app = emat_app-emat_appRef; end
    
    % if necessary, interpolate to standard framerate here before
    % proceeding
    if data(i).framerate ~= sframerate
        tvec_or         = data(i).framerate * tvec_app;
        tvec_interpol   = sframerate * tvec_standard(1,:);
        ivec_or         = emat_app;
        ivec_interpol   = interp1(tvec_or,ivec_or',tvec_interpol);
        emat_app = ivec_interpol';
    end
    
    % determine average and error
    evec_app = nanmean(emat_app);
    svec_app = nanstd(emat_app);
    
    % replace the alignment point value if desired
    if voidAlignPoint==1
        p0_app = find(tvec_app==0);
        evec_app(p0_app) = (evec_app(p0_app-1)+evec_app(p0_app+1))/2;
        svec_app(p0_app) = (svec_app(p0_app-1)+svec_app(p0_app+1))/2;
    end
    
    sn(i) = px;
    
    
    
    % ==========================
    % DISAPPEARANCE-aligned
    tvec_disapp = tvec_standard(2,:);
    emat_disapp = nan*zeros(px,length(tvec_disapp));
    if reference==1, emat_disappRef = emat_disapp; end
    for p=1:px
        tloc = paraT2(p,:);
        for t=1:length(tvec_disapp)
            ctau = tvec_disapp(t);
            tfind = find(tloc==ctau);
            if ~isempty(tfind)
                emat_disapp(p,t) = paraInt(p,tfind);
                % if reference parameter values exist, read these and
                % later subtract from data parameter value
                if reference==1
                    emat_disappRef(p,t) = paraInt2(p,tfind);
                end
            end
        end
    end
    % if reference positions exist, subtract their signal here
    if reference==1, emat_disapp = emat_disapp-emat_disappRef; end
    
    % if necessary, interpolate to standard framerate here before
    % proceeding
    if data(i).framerate ~= sframerate
        tvec_or         = data(i).framerate * tvec_disapp;
        tvec_interpol   = sframerate * tvec_standard(2,:);
        ivec_or         = emat_disapp;
        ivec_interpol   = interp1(tvec_or,ivec_or',tvec_interpol);
        emat_disapp = ivec_interpol';
    end
    
    % determine average and error
    evec_disapp = nanmean(emat_disapp);
    svec_disapp = nanstd(emat_disapp);
    
    % replace the alignment point value if desired
    if voidAlignPoint==1
        p0_disapp = find(tvec_disapp==0);
        evec_disapp(p0_disapp) = (evec_disapp(p0_disapp-1)+evec_disapp(p0_disapp+1))/2;
        svec_disapp(p0_disapp) = (svec_disapp(p0_disapp-1)+svec_disapp(p0_disapp+1))/2;
    end
    
    adlen = length(evec_disapp);
    
    
    % ==========================
    % weighting vector
    tvec_sig = tvec_standard(1,:)+tvec_standard(2,:);
    lambda = (restvector(4)/data(i).framerate)/3;
    sig_fun = 1./(1 + exp(-(tvec_sig/lambda)));
    slen = length(sig_fun);
    
    
    % ==========================
    % final result is aligned as desired
    if alignvar==3
        % weighted trace
        int_weighted    = sig_fun.*evec_disapp(adlen-slen+1:adlen) + (1-sig_fun).*evec_app(1:slen);
        error_weighted  = sig_fun.*svec_disapp(adlen-slen+1:adlen) + (1-sig_fun).*svec_app(1:slen);
    elseif alignvar==1
        int_weighted    = evec_app(1:slen);
        error_weighted  = svec_app(1:slen);
    elseif alignvar==2
        int_weighted    = evec_disapp(adlen-slen+1:adlen);
        error_weighted  = svec_disapp(adlen-slen+1:adlen);
    end
    
    
    % ==========================
    % ASSIGN RESULTS
    intMat(i,:)     = int_weighted;
    errorMat(i,:)   = error_weighted;
    tvec            = tvec_standard;
    
    % ==========================
    % plot result
    hold on;
    errorbar(data(i).framerate*tvec_standard(1,:),int_weighted,error_weighted/sqrt(px),'r');
    
    
end % of for i-loop


% write results to data structure
intRes.intMat   = intMat;
intRes.errorMat = errorMat;
intRes.sn       = sn;
intRes.tvec     = tvec;
intRes.framerate = sframerate;

% determine final average and standard error of the mean:

% total error variance sum: std of each movie squared, multiplied by number
% of trajs in the respective movie
varmat = (errorMat.^2).*repmat(sn',1,size(errorMat,2));

% final error of mean: root of total variance sum divided by total n,
% divided by sqrt(n)
varvec = sqrt(nansum(varmat,1)/sum(sn))/sqrt(sum(sn));

intRes.intAVE = nanmean(intMat);
intRes.intSEM = varvec;


% plot final results results
if alignvar==2
    tvec = intRes.tvec(2,:);
else
    tvec = intRes.tvec(1,:);
end
ivec    = intRes.intAVE;
evec    = intRes.intSEM;
fr      = intRes.framerate;

figure
errorbar(fr*tvec, ivec, evec, 'k-');
hold on;
plot(fr*tvec, ivec, 'r-', 'LineWidth', 2);
xlabel('time [s]');
ylabel('parameter intensity above BG');


function fileName = getParameterFileName(sourcePath, pfName)
if (exist([sourcePath  pfName], 'file')==2)
    fileName = [sourcePath pfName];
else
    [fileName, filePath] = uigetfile('*.mat', ['Select ' pfName], sourcePath);
    fileName = [filePath fileName];
end;