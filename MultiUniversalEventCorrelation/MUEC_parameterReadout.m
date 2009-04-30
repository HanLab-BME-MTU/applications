function [intRes] = MUEC_parameterReadout(data, restvector, channel, reference, sframerate, align, voidp, parname); 
% Multichannel Universal Event Correlator (MUEC) to read out image
% intensities or parameter values at specified object locations
% 
% SYNOPSIS  [intRes] = MUEC_parameterReadout(data, restvector, channel,...
%                           reference, sframerate, align, voidp, parname);
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
%           voidp:      OPTIONAL  - indicates whether or not (1/0) to void 
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
% last modified: Dinah Loerke, March 18, 2009


warning off;
od = cd;

%% test inputs and set defaults
if isfield(data,'framerate')
    scvec = data(1).framerate;
else
    scvec = 1;
end

lref = 0;
if nargin>3
    if reference==1
        lref = 1;
    end
end

sframe = 2;
if nargin>4
    if ~isempty(sframerate)
        sframe = sframerate; 
    end
end

alignvar = 3;
if nargin>5
    if strcmp(align,'left')
       alignvar = 1; 
    elseif strcmp(align,'right')
       alignvar = 2;
    end
end

voidAlignPoint = 0;
if nargin>6
    if voidp==1
        voidAlignPoint = 1;
    end
end

pname = 0;
if nargin>7
    if ischar(parname)
        fileNameStump = parname;
        paramName = [parname,'.mat'];
        paramRefName = [parname,'_Ref.mat'];
        pname = 1;
    end
end


% trajectory min/max length, standard time vector
padding = 15;
tlen_min = restvector(4)/sframe;
tlen_max = restvector(5)/sframe;
tlen_st  = floor((tlen_min+tlen_max)/2);
tvec_standard(1,:) =  [-padding:tlen_st+padding-1];
tvec_standard(2,:) =  [-(tlen_st+padding-1):padding];




%% ========================================================================
%  load parameter files for all movies
% =========================================================================

for i=1:length(data)
    
    % set path where to look for parameter/intensity file
    path = data(i).source;
    if nargin>2
        if ~isempty(channel)
            path = getfield(data,{1,i},channel);
        end
    end
            
    % if path is directory, change to this directory and uiget the files
    if exist(path)==7
        cd(path); 
        
        % if parameter file name has been specified, and if a file of that
        % name exists at the specified location, use this file name
        haveName = 0;
        if (pname==1) 
            currCompleteFileName = strcat(path, filesep, paramName); 
            if exist(currCompleteFileName)==2
                completeFileName(i).name = currCompleteFileName;
                haveName = 1;
            end
        end
        
        % else define file name via user input
        if ~haveName
            [fileName, filePath] = uigetfile('.mat','Select parameter file'); 
            completeFileName(i).name = strcat(filePath, fileName);  
        end
        
        % same procedure for reference parameter file name, if reference
        % parameter is used
        if lref==1
            haveNameRef = 0;
            if (pname==1) 
                currCompleteFileNameRef = strcat(path, filesep, paramRefName); 
                if exist(currCompleteFileNameRef)==2
                    completeFileName2(i).name = currCompleteFileNameRef;
                    haveNameRef = 1;
                end
            end
            
            if ~haveNameRef
                [fileName2, filePath2] = uigetfile('.mat','Select reference file'); 
                completeFileName2(i).name = strcat(filePath2, fileName2); 
            end
        end

    % if path is file name, use this file
    elseif exist(path)==2
        completeFileName(i).name = path;             
    else
        error('no parameter file found');
    end
            
end





%% ========================================================================
%  extract parameter information from files for all movies
% =========================================================================

% loop over all movies
for i=1:length(data)
    
    fprintf('movie #%02d',i);
    fprintf('\n');
    
    % load trackinfo file
    cpath = [data(i).source,filesep,'TrackInfoMatrices'];
    lftinfopath = [data(i).source,filesep,'LifetimeInfo'];
        
    framerate = data(i).framerate;
    
    % extract the positions with the desired parameters (e.g. lifetime
    % range) from lft info data
    posvec = extractPos_vecFromLftRest(lftinfopath, restvector, framerate);
    
    % load parameter files
    filename = completeFileName(i).name; 
    loadfile = load(filename);
    loadfieldname = fieldnames(loadfile);
    paramMat = getfield(loadfile,loadfieldname{1});  
    
    % extract parameter values
    paramMat_use = paramMat(posvec,:,:);   
    
    paraInt = paramMat_use(:,:,1);
    paraT1  = paramMat_use(:,:,2);
    paraT2  = paramMat_use(:,:,3);
    [px,pf] = size(paraInt);
    val_max = max(paraT1(:));
    val_min = min(paraT2(:));
    
    %%===================================================================
    % if reference parameter values are used, load and calculate them here
    if lref==1
        filename2 = completeFileName2(i).name; 
        loadfile2 = load(filename2);
        loadfieldname2 = fieldnames(loadfile2);
        paramMat2 = getfield(loadfile2,loadfieldname2{1});  
        % determine how many reference points exist for each data point
        % (e.g. 8 or 10 background points per data point)
        [m1x,m1y] = size(paramMat);
        [m2x,m2y] = size(paramMat2);
        rfac = m2x/m1x;
        % if rfac>1 (e.g. rfac=8), then the assumption is that the 8 ref 
        % points for data point number one are in the rows 1-8, those for
        % data point number 2 in rows 9-16, etc.
        % for easier use in the following function, the reference data are
        % immediately averaged here, i.e. they are compressed from rfac>1
        % to a matrix the same size as paramMat_use(:,:,1)
        paraInt2 = nan*paraInt;
        if rfac>1
            for k=1:length(posvec)
                % current data point
                cpos = posvec(k);
                % all reference points for this data point
                rpos = [(rfac*(cpos-1)+1):(rfac*cpos)];
                % parameter values for these reference points
                rpar = paramMat2(rpos,:,1);
                % average
                paraInt2(k,:) = nanmean(rpar);
            end
        else
            paraInt2 = paramMat2(posvec,:,1);
        end
    end
 
    
%% ===================================================================
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
    if lref==1, emat_appRef = emat_app; end
    for p=1:px
        tloc = paraT1(p,:);
        for t=1:length(tvec_app)
            ctau = tvec_app(t);
            tfind = find(tloc==ctau);
            if ~isempty(tfind)
                emat_app(p,t) = paraInt(p,tfind);
                % if reference parameter values exist, read these and
                % later subtract from data parameter value
                if lref==1
                    emat_appRef(p,t) = paraInt2(p,tfind); 
                end
            end
        end
    end
    % if reference positions exist, subtract their signal here
    if lref==1, emat_app = emat_app-emat_appRef; end
    
    % if necessary, interpolate to standard framerate here before
    % proceeding
    if framerate ~= sframe
        tvec_or         = framerate * tvec_app;
        tvec_interpol   = sframe * tvec_standard(1,:);
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
    if lref==1, emat_disappRef = emat_disapp; end
    for p=1:px
        tloc = paraT2(p,:);
        for t=1:length(tvec_disapp)
            ctau = tvec_disapp(t);
            tfind = find(tloc==ctau);
            if ~isempty(tfind)
                emat_disapp(p,t) = paraInt(p,tfind);
                % if reference parameter values exist, read these and
                % later subtract from data parameter value
                if lref==1
                    emat_disappRef(p,t) = paraInt2(p,tfind); 
                end
            end
        end
    end
    % if reference positions exist, subtract their signal here
    if lref==1, emat_disapp = emat_disapp-emat_disappRef; end
    
    % if necessary, interpolate to standard framerate here before
    % proceeding
    if framerate ~= sframe
        tvec_or         = framerate * tvec_disapp;
        tvec_interpol   = sframe * tvec_standard(2,:);
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
    
    
    % ==========================
    % ASSIGN RESULTS
    intMat(i,:)     = int_weighted;
    errorMat(i,:)   = error_weighted;
    tvec            = tvec_standard;
    
    % ==========================
    % plot result
    hold on;
    errorbar(framerate*tvec_standard(1,:),int_weighted,error_weighted/sqrt(px),'r');
    
    
end % of for i-loop




% return to original directory
cd(od);

% write results to data structure
intRes.intMat   = intMat;
intRes.errorMat = errorMat;
intRes.sn       = sn;
intRes.tvec     = tvec;
intRes.framerate = sframe;

% determine final average and standard error of the mean:

% total error variance sum: std of each movie squared, multiplied by number
% of trajs in the respective movie
varmat = (errorMat.^2).*repmat(sn',1,size(errorMat,2));

% final error of mean: root of total variance sum divided by total n,
% divided by sqrt(n)
varvec = sqrt(nansum(varmat,1)/sum(sn))/sqrt(sum(sn));
    
intRes.intAVE   = nanmean(intMat);
intRes.intSEM   = varvec;


% plot final results results
tvec    = intRes.tvec(1,:);
if alignvar==2, tvec = intRes.tvec(2,:); end
ivec    = intRes.intAVE;
evec    = intRes.intSEM;
fr      = intRes.framerate;

figure
errorbar(fr*tvec,ivec,evec,'k-');
hold on;  
plot( fr*tvec,ivec,'r-','LineWidth',2 );
xlabel('time (sec)');
ylabel('parameter intensity above BG');




end % of function
    
      
