function [intRes] = MUEC_parameterReadoutCargo(data, restvector, channel, reference, framerate, align, voidAlignPoint, parname_c0, parname_c1)
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
%           framerate: OPTIONAL - standard framerate to which all other
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
%                       .nt         = number of individual trajectories in
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
% Last modified: Francois Aguet, 12/14/2009


%==================================
% 0. Test inputs and set defaults
%==================================
if (nargin <= 3)
    reference = 0;
end;
if (nargin <= 4)
    framerate = 2;
end;
if (nargin >= 6)
    if strcmp(align,'left')
       alignvar = 1; 
    elseif strcmp(align,'right')
       alignvar = 2;
    else
       alignvar = 3;
    end; 
else % check
    alignvar = 3;
end;
if (nargin < 7)
    voidAlignPoint = 0;
end;
if (nargin>7 && ischar(parname_c0))
    paramName_c0 = [parname_c0 '.mat'];
    paramNameRef_c0 = [parname_c0 '_Ref.mat'];
else
    paramName_c0 = 'parameterMat.mat';
    paramNameRef_c0 = 'parameterMat_Ref.mat';
end;
if (nargin>8 && ischar(parname_c1))
    paramName_c1 = [parname_c1 '.mat'];
    paramNameRef_c1 = [parname_c1 '_Ref.mat'];
else
    paramName_c1 = 'SA_red_ParameterMat_.mat';
    paramNameRef_c1 = 'SA_red_ParameterMat__Ref.mat';
end;


% =======================================================================
% 0. Load parameter files for all movies
% =======================================================================
nMovies = length(data);
fileName_c0 = cell(1:nMovies);
fileNameRef_c0 = cell(1:nMovies);
fileName_c1 = cell(1:nMovies);
fileNameRef_c1 = cell(1:nMovies);

for i=1:nMovies
    % set search path for parameter/intensity file
    if nargin>2 && ~isempty(channel)
        sourcePath = [getfield(data, {1,i}, channel) filesep];
    else
        sourcePath = [data(i).source filesep];
    end;
    
    if exist(sourcePath, 'dir')==7
        
        % Load parameter files for both channels
        fileName_c0{i} = getParameterFileName(sourcePath, paramName_c0);
        fileName_c1{i} = getParameterFileName(sourcePath, paramName_c1);
        
        % Load reference parameter files for both channels
        if reference==1
            fileNameRef_c0{i} = getParameterFileName(sourcePath, paramNameRef_c0);
            fileNameRef_c1{i} = getParameterFileName(sourcePath, paramNameRef_c1);
        end;

    % if path is file name, use this file
    elseif exist(sourcePath, 'file')==2
        fileName_c0{i} = sourcePath;             
    else
        error('No parameter file found.');
    end;
end;

% ========================================================================
% extract parameter information from files for all movies
% ========================================================================

% loop over all movies
h1 = figure;
h2 = figure;
h3 = figure;
for i=1:nMovies
    
    fprintf('Sequence #%d\n', i);
            
    % extract the positions with the desired parameters
    % retrieves the indexes of the positions corresponding to the parameters (e.g. lifetime range from 'lftInfo.mat')
    posvec = extractPos_vecFromLftRest([data(i).source filesep 'LifetimeInfo'], restvector, data(i).framerate);
    nt = length(posvec); % number of CCPs in cohort
    fprintf('Number of tracks in cohort: %d\n', nt);
    % load parameters stored as iMat_obj in parameterMat.mat
    paramMat_c0 = load(fileName_c0{i});
    paramMat_c0 = paramMat_c0.iMat_obj;
    paramMat_c0 = paramMat_c0(:,:,:); % extract parameter values
%     paramMat_c1 = load(fileName_c1{i});
%     paramMat_c1 = paramMat_c1.iMat_obj;
%     paramMat_c1 = paramMat_c1(posvec,:,1);
    
    paraInt = paramMat_c0(posvec,:,1);
    paraT1  = paramMat_c0(posvec,:,2);
    paraT2  = paramMat_c0(posvec,:,3);
    %intRes(i).nt = nt;
    %intRes(i).framerate = data(i).framerate; % useless, remove -> need new format to store prms        
        
    %[px,pf] = size(paraInt);
    %val_max = max(paraT1(:));
    %val_min = min(paraT2(:));
    
    %%===================================================================
    % if reference parameter values are used, load and calculate them here
    if reference==1
        paramMatRef_c0 = load(fileNameRef_c0{i});
        paramMatRef_c0 = paramMatRef_c0.iMat_ref;
        
        %paramMatRef_c1 = load(fileNameRef_c1{i});
        %paramMatRef_c1 = paramMatRef_c1.iMat_ref;
        
        nbgPixels = size(paramMatRef_c0,1)/size(paramMat_c0, 1);
        % determine how many reference points exist for each data point
        % (e.g. 8 or 10 background points per data point)

        % if rfac>1 (e.g. rfac=8), then the assumption is that the 8 ref 
        % points for data point number one are in the rows 1-8, those for
        % data point number 2 in rows 9-16, etc.
        % for easier use in the following function, the reference data are
        % immediately averaged here, i.e. they are compressed from rfac>1
        % to a matrix the same size as paramMat_use(:,:,1)
        if nbgPixels>1
            bgMean_c0 = NaN(size(paraInt));
            %bgMean_c1 = NaN(size(paraInt));
            for k = 1:nt
                %posvec(k)
                %nbgPixels
                % all reference points for this data point
                rpos = (nbgPixels*(posvec(k)-1)+1):(nbgPixels*posvec(k));
                % parameter values for these reference points
                bgMean_c0(k,:) = nanmean(paramMatRef_c0(rpos,:,1));
                %bgMean_c1(k,:) = nanmean(paramMatRef_c1(rpos,:,1));
            end
        else
            bgMean_c0 = paramMatRef_c0(posvec,:,1);
            %bgMean_c1 = paramMatRef_c1(posvec,:,1);
        end
    end
 
    
    %load([sourcePath 'TrackInfoMatrices/trackInfo.mat']);
    %tracks = trackInfo(posvec,:);
    %load([sourcePath 'TrackInfoMatrices' filesep 'trackedFeatures.mat']);
    %tracks = trackedFeatureInfo(posvec,:);    
    %ints = tracks(:,4:8:end);  
    
    
    
    
    % ===================================================================
    % data alignment:
    % "chaotic evil means never having to say you're sorry"
    %
    % the data are a conglomerate of intensity trajectories of different
    % lengths; they are averaged by first aligning them to the appearance time
    % point, then to the disappearance time point, and then forming a weighted
    % average - this procedure conserves the sharp increase in e.g. intensity
    % during the nucleation and internalization phase
    
    % trajectory min/max length, standard time vector
    padding = 15;
    tlen_min = restvector(4)/framerate;
    tlen_max = restvector(5)/framerate;
    tlen_st  = floor((tlen_min+tlen_max)/2);
    tvec_standard(1,:) = -padding:tlen_st+padding-1;
    tvec_standard(2,:) = -(tlen_st+padding-1):padding;

    load([data(i).source 'TrackInfoMatrices' filesep 'cargoStatus.mat']);
    cargoStatus = cargoStatus(posvec)';
    valid = valid(posvec)';

    cidx = cargoStatus==1 & valid==1;
    nidx = cargoStatus==0 & valid==1;
    sum(cidx)
    sum(nidx)
    length(posvec)
    
    paraT1_cargo = paraT1(cidx,:);
    paraT2_cargo = paraT2(cidx,:);
    paraInt_cargo = paraInt(cidx,:);
    bgMean_c0_cargo = bgMean_c0(cidx,:);
    paraT1_nocargo = paraT1(nidx,:);
    paraT2_nocargo = paraT2(nidx,:);
    paraInt_nocargo = paraInt(nidx,:);
    bgMean_c0_nocargo = bgMean_c0(nidx,:);

    [int_weighted, error_weighted] = alignCohort(data(i), restvector, reference, framerate, voidAlignPoint, alignvar, paraT1, paraT2, paraInt, bgMean_c0, tvec_standard);
    [int_weighted_cargo, error_weighted_cargo] = alignCohort(data(i), restvector, reference, framerate, voidAlignPoint, alignvar, paraT1_cargo, paraT2_cargo, paraInt_cargo, bgMean_c0_cargo, tvec_standard);
    [int_weighted_nocargo, error_weighted_nocargo] = alignCohort(data(i), restvector, reference, framerate, voidAlignPoint, alignvar, paraT1_nocargo, paraT2_nocargo, paraInt_nocargo, bgMean_c0_nocargo, tvec_standard);
    
    sn(i) = nt;
    sn_c(i) = length(cidx);
    sn_n(i) = length(nidx);
    
    % ==========================
    % ASSIGN RESULTS
    intMat(i,:)     = int_weighted;
    errorMat(i,:)   = error_weighted;
    tvec            = tvec_standard;
    
    intMat_cargo(i,:)     = int_weighted;
    errorMat_cargo(i,:)   = error_weighted;
    tvec_cargo           = tvec_standard;
    
    intMat_nocargo(i,:)     = int_weighted;
    errorMat_nocargo(i,:)   = error_weighted;
    tvec_nocargo            = tvec_standard;
    
    %statusV{i} = cargoStatus;    
    % ==========================
    
    % separate based on cargo
%     whos
%     int_cargo = int_weighted(cargoStatus==1 & valid==1);
%     int_nocargo = int_weighted(cargoStatus==0 & valid==1);
%     
    
    
    
    % plot result
    figure(h1)
    errorbar(framerate*tvec_standard(1,:),int_weighted,error_weighted/sqrt(nt),'r');
    hold on;

    
    figure(h2) % cargo
    errorbar(framerate*tvec_standard(1,:),int_weighted_cargo,error_weighted_cargo/sqrt(length(cidx)),'r');
    hold on;

    figure(h3) % no cargo
    errorbar(framerate*tvec_standard(1,:),int_weighted_nocargo,error_weighted_nocargo/sqrt(length(nidx)),'r');
    hold on;
    
end % of for i-loop

figure(h1)
title('All pits');
figure(h2)
title('Pits containing cargo');
figure(h3)
title('Empty pits');



% write results to data structure
intRes.intMat   = intMat;
intRes.errorMat = errorMat;
intRes.sn       = sn;
intRes.tvec     = tvec;
intRes.framerate = framerate;
%intRes.cargoStatus = statusV;

% determine final average and standard error of the mean:

% total error variance sum: std of each movie squared, multiplied by number
% of trajs in the respective movie
varmat = (errorMat.^2).*repmat(sn',1,size(errorMat,2));

varmat_cargo = (errorMat_cargo.^2).*repmat(sn_c',1,size(errorMat_cargo,2));
varmat_nocargo = (errorMat_nocargo.^2).*repmat(sn_n',1,size(errorMat_nocargo,2));


% final error of mean (over all movies):
% root of total variance sum divided by total n, divided by sqrt(n)
%varvec = sqrt(nansum(varmat,1)/sum(nt))/sqrt(sum(nt));
varvec = sqrt(nansum(varmat,1)/sum(sn))/sqrt(sum(sn));
varvec_cargo = sqrt(nansum(varmat_cargo,1)/sum(sn_c))/sqrt(sum(sn_c));
varvec_nocargo = sqrt(nansum(varmat_nocargo,1)/sum(sn_n))/sqrt(sum(sn_n));



intRes.intAVE   = nanmean(intMat);
intRes.intSEM   = varvec;


% plot final results
if alignvar==2
    tvec = intRes.tvec(2,:);
    tvec_cargo = tvec_cargo(2,:);
    tvec_nocargo = tvec_nocargo(2,:);
else
    tvec = intRes.tvec(1,:);
    tvec_cargo = tvec_cargo(1,:);
    tvec_nocargo = tvec_nocargo(1,:);
end
ivec    = intRes.intAVE;
evec    = intRes.intSEM;
fr      = intRes.framerate;

% TO DO: remove redunant variables
evec_cargo = varvec_cargo;
evec_nocargo = varvec_nocargo;
ivec_cargo = nanmean(intMat_cargo);
ivec_nocargo = nanmean(intMat_nocargo);


figure
errorbar(fr*tvec,ivec,evec,'k-');
hold on;  
plot( fr*tvec,ivec,'r-','LineWidth',2 );
xlabel('time (sec)');
ylabel('parameter intensity above BG');
title('All pits');

figure
errorbar(fr*tvec_cargo,ivec_cargo,evec_cargo,'k-');
hold on;  
plot( fr*tvec_cargo,ivec_cargo,'r-','LineWidth',2 );
xlabel('time (sec)');
ylabel('parameter intensity above BG');
title('Pits containing cargo');

figure
errorbar(fr*tvec_nocargo,ivec_nocargo,evec_nocargo,'k-');
hold on;  
plot( fr*tvec_nocargo,ivec_nocargo,'r-','LineWidth',2 );
xlabel('time (sec)');
ylabel('parameter intensity above BG');
title('Empty pits');




function fileName = getParameterFileName(sourcePath, pfName)
if (exist([sourcePath  pfName], 'file')==2)
    fileName = [sourcePath pfName];
else
    [fileName, filePath] = uigetfile('*.mat', ['Select ' pfName], sourcePath);
    fileName = [filePath fileName];
end;



function [int_weighted, error_weighted] = alignCohort(data, restvector, reference, framerate, voidAlignPoint, alignvar, paraT1, paraT2, paraInt, bgMean_c0, tvec_standard)

% ==========================
% APPEARANCE-aligned
tvec_app = tvec_standard(1,:);
nt = size(paraT1,1); % #tracks in cohort
emat_app = NaN(nt,length(tvec_app));
if reference==1
    emat_appRef = emat_app;
end;
for k = 1:nt % loop through tracks
    tloc = paraT1(k,:);
    for t = 1:length(tvec_app)
        tfind = find(tloc == tvec_app(t));
        if ~isempty(tfind)
            emat_app(k,t) = paraInt(k,tfind);
            % if reference parameter values exist, read these and
            % later subtract from data parameter value
            if reference==1
                emat_appRef(k,t) = bgMean_c0(k,tfind);
            end
        end
    end
end
% if reference positions exist, subtract their signal here
if reference==1
    emat_app = emat_app - emat_appRef;
end;

% if necessary, interpolate to standard framerate here before proceeding
if data.framerate ~= framerate
    tvec_or         = data.framerate * tvec_app;
    tvec_interpol   = framerate * tvec_standard(1,:);
    ivec_or         = emat_app;
    ivec_interpol   = interp1(tvec_or, ivec_or', tvec_interpol);
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

%sn(i) = nt;



% ==========================
% DISAPPEARANCE-aligned
tvec_disapp = tvec_standard(2,:);
emat_disapp = NaN(nt,length(tvec_disapp));
if reference==1, emat_disappRef = emat_disapp; end
for p=1:nt
    tloc = paraT2(p,:);
    for t=1:length(tvec_disapp)
        tfind = find(tloc == tvec_disapp(t));
        if ~isempty(tfind)
            emat_disapp(p,t) = paraInt(p,tfind);
            % if reference parameter values exist, read these and
            % later subtract from data parameter value
            if reference==1
                emat_disappRef(p,t) = bgMean_c0(p,tfind);
            end
        end
    end
end
% if reference positions exist, subtract their signal here
if reference==1, emat_disapp = emat_disapp-emat_disappRef; end

% if necessary, interpolate to standard framerate here before
% proceeding
if data.framerate ~= framerate
    tvec_or         = data.framerate * tvec_disapp;
    tvec_interpol   = framerate * tvec_standard(2,:);
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
lambda = (restvector(4)/framerate)/3;
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

