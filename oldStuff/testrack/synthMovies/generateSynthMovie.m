function [synthMovie, r3dMovieHeader, dataProperties, inputSlist, background] = generateSynthMovie(slist,dataProperties,SNR,options,synthMovie)
%generateSyntMovie generates a synthetic yeast movie from coordinates
%
%SYNOPSIS [synthMovie, movieHeader, dataProperties, slist, idlist] = generateSynthMovie(slist,size,SNR,options,synthMovie)
%
%INPUT    slist (opt): real or synthetic slist containing at least the fields slist(t).sp(i).cord
%                       [1x3,pix] and slist(t).sp(i).amp (opt) [betw. 0 and 1 in steps of 2e-12]
%                      Alternatively, a cell of length nTimepoints can be
%                      given, which for each t contains an array of size
%                      nSpots by 4: x,y,z,amp. (x,y,z in microns)
%         dataProperties (opt): pixelSize: default 0.048,0.2
%                               movieSize: default 64,64,18,10
%                               frameTime: default 1 stack/sec, 0.03 sec/slice
%                               newName     : will become name, default ['synthMovie-',nowString]
%         SNR (opt):  SNR (signal/sigmaNoise) {5}
%         options (opt): .movie   1: generate movie, 2: add noise, {3}: generateMovie&addNoise
%                        .psf     1: gauss, {2}: theoretical psf
%                        .microscope [wavelength, NA, magnification] default: [0.525, 1.4, 100]
%                        .background [value] default: 0. if -1, the min of the movie won't be 0, but
%                                    (probably) negative
%                        .saveDir   directory into which the movie should be saved. [directory/'' or 0 (no save)/{'ask'}]
%                        .spots      .nsp number of spots {2}
%                                    .intRatio intensity ratios (relative to first spot, descending) {[1,0.8]}
%                        .sampling  size of the sampling cube in microns {0.001}
%
%OUTPUT   synthMovie: synthetic data (.r3c)
%         movieHeader: header of .r3c-movie
%         dataProperties: dataProperties of .r3c-movie
%         inputSlist: slist containing location and amplitude of generated spots
%         background: background value of the movie
%         TO BE ADDED: slist, idlist
%
%c: jonas, 7/03, based on generateGaussMovieData by dT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maxMem: maximum size of a frame. Currently, 500MB
MAXMEM = 524288000;
% whether to resample or not
resample =1;

def_pixelSize = [0.05, 0.05, 0.05]; % [0.048, 0.048, 0.2]
def_movieSize = [10, 10, 30, 1]; % [64, 64, 18, 10]
def_sampling = 0.001;

%-------------------test input and assign defaults-------------------

if nargin < 1 || isempty(slist)
    calcSlist = 1;
elseif isstruct(slist)
    calcSlist = 2;
elseif iscell(slist)
    calcSlist = 0;
else
    error('wrong format for slist')
end

%dataProperties (do this first because we need movieSize for slist)
if nargin<2|isempty(dataProperties)
    dataProperties = [];
end
%now fill all fields that haven't been filled yet
if isfield(dataProperties,'pixelSize')
    dataProperties.PIXELSIZE_XY = dataProperties.pixelSize(1);
    dataProperties.PIXELSIZE_Z  = dataProperties.pixelSize(2);
    rmfield(dataProperties,'pixelSize')
end
if ~isfield(dataProperties,'PIXELSIZE_XY')
    dataProperties.PIXELSIZE_XY = def_pixelSize(1);
end
if ~isfield(dataProperties,'PIXELSIZE_Z')
    dataProperties.PIXELSIZE_Z = def_pixelSize(3);
end
% store pixelsize in a reasonably short variable
pixelSize = [dataProperties.PIXELSIZE_XY, ...
    dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z];

if ~isfield(dataProperties,'movieSize')
    dataProperties.movieSize = def_movieSize;
end


%define movieLength, stackSize
movieLength = dataProperties.movieSize(4);
stackSize = dataProperties.movieSize(1:3);

if ~isfield(dataProperties,'frameTime')
    %1 sec steps
    meanTime = [1:movieLength]'*1;
    %center acquisition timepoints around the mean, allow 0.03 sec per slice
    stackTime = ((stackSize(3)-1)/2-[stackSize(3)-1:-1:0])*0.03;
    frameTime = repmat(meanTime,1,stackSize(3))+repmat(stackTime,movieLength,1);
    %frameTime starts with 0
    dataProperties.frameTime = frameTime-frameTime(1);
end

if ~isfield(dataProperties,'newName')
    dataProperties.name = ['synthMovie-',nowString];
else
    dataProperties.name = dataProperties.newName;
end

%fill in the rest
if ~isfield(dataProperties,'CH_MAXSLOPE')
    dataProperties.LENSID=12003;
    dataProperties.NA=1.4000;
    dataProperties.WVL=0.5250;
    dataProperties.timeLapse=1;
    dataProperties.expTime=0.0010;
    dataProperties.NDfilter=0.5000;
    dataProperties.cellCycle=0;
    dataProperties.strains=0;
    dataProperties.drugs=0;
    dataProperties.temperature={'30°C'};
    dataProperties.crop=[];
    dataProperties.CH_MAXSLOPE=0.0010;
    dataProperties.F_TEST_PROB=0.9000;
    dataProperties.IDopt= [];
    dataProperties.PATCHSIZE=7;
    dataProperties.CH_MAXNUMINTERV=1000;
    dataProperties.OVERLPSIZE=[15 15 15];
    dataProperties.sigmaCorrection=[1 1];
    dataProperties.split=[];
    dataProperties.MAXSPOTS=5;
    dataProperties.T_TEST_PROB=0.0500;
    dataProperties.help=['synthMovie'];
end

%options
%.movie   1: generate movie, 2: add noise, {3}: generateMovie&addNoise
%.psf     1: gauss, {2}: theoretical psf
%.microscope [wavelength, NA] default: [0.525, 1.4]
%.background [value] default: 0
%.saveDir   directory into which the movie should be saved. [directory/'' (no save)/{'ask'}]
%.spots      .nsp number of spots {2} .intRatio intensity ratios (relative to first spot, descending) [1,0.8]
%.sampling  size of the sampling cube in microns {0.001}

if nargin<4|isempty(options)
    options = [];
end

if ~isfield(options,'movie')
    options.movie = 3;
elseif options.movie<1|options.movie>3
    error('options.movie out of bounds')
end
jobList = bsum2bvec(options.movie);

if ~isfield(options,'psf')
    options.psf = 2;
elseif options.psf<1|options.psf>2
    error('options.psf out of bounds')
end


if isfield(options,'microscope')
    if length(options.microscope)==3
        wvl = options.microscope(1);
        NA = options.microscope(2);
        magnification = options.microscope(3);
        if wvl>10
            error('wavelength has to be in microns!')
        else
            dataProperties.WVL = wvl;
        end
        if NA<=0
            error('NA out of bounds')
        else
            dataProperties.NA = NA;
        end
    else
        error('wrong format of options.microscope')
    end
else %take default values
    wvl = dataProperties.WVL;
    NA = dataProperties.NA;
    magnification = 100;
end

% now that we know NA and wvl, we can calculate filterparms
[FT_XY, FT_Z] = calcFilterParms(...
    wvl,NA,1.51,'gauss',[], pixelSize(2:3));
patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];

if ~isfield(options,'background')
    options.background = 0;
elseif (options.background<0 & options.background~=-1)|options.background>0.9
    error('options.background out of bounds')
end

if ~isfield(options,'save')
    options.save = 'ask';
elseif ~isstr(options.save)|isempty(options.save)
    options.save = '';
elseif ~strcmp(options.save,'ask')
    if ~isdir(options.save)
        error('options.save is not a directory')
    end
end

if ~isfield(options,'spots')
    options.spots.nsp = 2;
    options.spots.intRatio = [1 0.8];
elseif ~isfield(options.spots,'nsp')
    options.spots.nsp = 2;
elseif ~isfield(options.spots,'intRatio')
    options.spots.intRatio = [1:-0.2:1.2-options.spots.nsp*0.2];
end
nsp = options.spots.nsp;
intRatio = options.spots.intRatio;
if nsp<1|nsp>9
    error('options.spots.nsp out of bounds')
end
if nsp>dataProperties.MAXSPOTS
    dataProperties.MAXSPOTS = nsp;
end

% sampling
if ~isfield(options,'sampling')
    sampling = def_sampling;
elseif any(options.sampling > pixelSize)
    error('with the current options.sampling, the movie would be undersampled during generation')
else
    sampling = options.sampling;
end




%----- slist ---------
% all tags will be listed as {[x,y,z,a]}. In case the input is in the form
% of this list, we do not generate an slist from it, because it is then
% likely to contain a large number of tags that can not be handled by the
% linker. If no list is given, we either take the slist or generate one.

%place tags
if calcSlist % only do something here if we get/produce an slist
    if calcSlist == 1 % and not 2
        for t = 1:movieLength
            for i = 1:nsp
                %place spots randomly, make sure they are never too close to border (whithin 5 5 3 pix)
                slist(t).sp(i).cord = (rand(1,3).*(stackSize-[10 10 6]))+[5 5 3];
            end
        end
    else
        nsp = length(slist(1).sp);
        options.spots.intRatio = [1:-0.2:1.2-options.spots.nsp*0.2];
        movieLength = length(slist);
        dataProperties.movieSize(4) = movieLength;
    end

    %make sure slist does not start empty
    tstart = 1;
    while isempty(slist(tstart)) && tstart<=length(slist)
        tstart = tstart+1;
    end
    if tstart > length(slist)
        error('invalid slist (no entry)')
    end
    slist = slist(tstart:end);
    movieLength = movieLength - tstart + 1;
    dataProperties.movieSize(4) = movieLength;

    %set amplitudes
    if ~isfield(slist(1).sp(1),'amp')
        DYNAMICRANGE = 0.003;
        for t = 1:movieLength
            for i = 1:nsp
                %set amplitudes; currently no bleaching or random fluctuations
                slist(t).sp(i).amp = options.spots.intRatio(i)*DYNAMICRANGE;
            end
        end
    end




    % read slist into pointCell
    pointCell = cell(movieLength,1);
    for t = 1:movieLength
        tmpPt = [];
        for i = 1:length(slist(t).sp)
            % convert to matlab coordinates and to microns
            % (independent of sampling)
            tmpPt = [tmpPt; ...
                slist(t).sp(i).cord([2,1,3]).*pixelSize,...
                slist(t).sp(i).amp];
        end
        pointCell{t} = tmpPt;
    end

else
    pointCell = slist;
    % read movieLength (assume user knows what he/she's doing)
    movieLength = length(pointCell);

    % check amp here
    allPoints = cat(1,pointCell{:});
    ampList = allPoints(:,4);
    if any(ampList>1)
        error('not normalized amplitudes!')
    else
        DYNAMICRANGE = max(ampList(:));
    end
    clear allPoints ampList

end % if calcSlist

%snr
if nargin<3 | isempty(SNR)
    SNR = 10;
elseif SNR<1
    error('snr below 1')
end

if exist('synthMovie')
    if ~isempty(synthMovie)
        jobList = 2;
        options.save = '';
        movieLength = size(synthMovie,5);
    else
        synthMovie=zeros([stackSize 1 movieLength]);
    end
end



%---------------end test input-----------------------------------------------------------



%----------------------create movie with psf------------------------

%init parms

%noise variance
noisevar=DYNAMICRANGE^2/(SNR^2);

%psf size (size where psf is 5% maxVal)
psfSze=roundOddOrEven(2*sqrt(-2*dataProperties.FILTERPRM(1:3).^2*log(0.05)),'odd','inf');
hPsf=floor(psfSze/2);

switch resample
    case 1

        % find a good sampling rate
        retry = 1;

        % size of output image in microns
        outSize = stackSize.*pixelSize;

        n = 1;
        while retry
            % try to generate two movies of tmpStackSize
            tmpStackSize = ceil(outSize/(sampling*n));
            %memSize = prod(tmpStackSize + floor(psfSze.*pixelSize/(sampling*n)));
            memSize = prod(tmpStackSize)*8;
            if memSize < MAXMEM
                retry = 0;
            else
                retry = 1;
                n = n+1;
            end
        end



        % warn user if sampling changed
        if n>1
            sampling = sampling*n;
            disp(sprintf('new sampling: %f',sampling))
        end

        % calculate new psf
        [sxy,sz] = calcFilterParms(...
            wvl,NA,1.51,'gauss',[], [sampling, sampling]);
        psfSze=roundOddOrEven(2*sqrt(-2*[sxy,sxy,sz].^2*log(0.05)),'odd','inf');
        hPsf=floor(psfSze/2);
        
        % WARNING: the origin is outside the image by half a pixel.
        % Changing the pixelSize therefore changes the placement of the
        % origin!
        correctOrigin = 0.5*(pixelSize-sampling*ones(1,3));
    case 0
        sampling = pixelSize;
        tmpStackSize = stackSize;
        sxy = dataProperties.FILTERPRM(1);
        sz  = dataProperties.FILTERPRM(3);
        correctOrigin = zeros(1,3);
end
%init psf
if options.psf == 2
    disp('sigmas of real psf could be a bit off')
    % from fitparmspos.mat @ c:\usr\thomann\matlab\dtOwn\wideFieldPsf\psfDataDir\psfdri1003
    psfparms=[ 0.9926    1.5032   25.2041   41.9328    0.0778    0.1049   -0.1882];
    psfob=psf(NA,magnification,wvl*1000,10000,1,psfparms(1),psfparms(2));
end


%loop time points
for t=1:movieLength

    % initialize Tmp block. Allow for psf out of bounds
    %% synthStackTmp = zeros([tmpStackSize+2*hPsf]);
    synthStackTmp = zeros(tmpStackSize);

    % generate raw movie

    if any(jobList==1)

        % read list of positions for current frame
        positions = pointCell{t};

        % switch according to psf. If gaussian, calc and position in stack. If
        % bessel, calc only once, then position in stack. possibly
        % interpolate for better results

        if options.psf == 2
            % calculate bessel psf. For later: generate with 1 pixel more
            % around, then interpolate within loop
            psfData=...
                makePSF(psfob,psfSze + ones(1,3),[1000*sampling(1),1000*sampling(end)],...
                [0 0 0]);
        end


        for i = 1:size(positions,1)

            % read position, transform to pix
            pos = (positions(i,1:3)-correctOrigin)./sampling;

            % calculate position within pixel
            posfullpix = floor(pos);
            % [0,0,0] is in the center of the center pixel
            subpixelShift = pos - posfullpix;
            switch options.psf
                case 1 %gaussian psf
                    % generate gauss
                    psfData = ...
                        GaussMask3D([sxy,sxy,sz],psfSze,subpixelShift);
                case 2
                    % interpolate bessel, set psf the correct size
                    xVector = 0.5:1:psfSze(1);
                    yVector = 0.5:1:psfSze(2);
                    zVector = 0.5:1:psfSze(3);
                    psfData = interp3(xVector,yVector,zVector,...
                        psfData,...
                        xVector+subpixelShift(1),...
                        yVector+subpixelShift(2),...
                        zVector+subpixelShift(3));
                    psfData = psfData(2:end-1,2:end-1,2:end-1);
            end

            % set amplitude
            psfData = psfData * positions(i,4);

            % place in stack
            synthStackTmp = ...
                pushstamp3d(synthStackTmp, psfData, posfullpix, 0, 1);

        end % for i=1:size(positions,1)

    end % if any jobList == 1


    % downsample and assign
    switch resample
        case 1
            synthMovie(:,:,:,1,t) = ...
                resampleMean(synthStackTmp,sampling .* ones(1,3), pixelSize);
        case 0
            synthMovie(:,:,:,1,t) = synthStackTmp;
    end

    % clear tmp block
    clear synthStackTmp;


% add dark noise (no shot noise)

if any(jobList == 2)
    %add noise
    nbd=synthMovie(:,:,:,1,t);
    nbd = nbd + sqrt(noisevar)*randn(size(nbd));
    synthMovie(:,:,:,1,t)=nbd;
end

end %for t=1:movieLength

%make sure the minimum value is zero (may introduce/change background)
if options.background == -1
    background = 0;
else
    minInt = min(synthMovie(:));
    if minInt+options.background>=0
        background = options.background;
    else
        background = options.background - min(synthMovie(:));
    end
end
synthMovie = synthMovie + background;


%---------------------end create movie-----------------------------------------------

%---------------------assign output--------------------------------------------------
%synthMovie: already assigned
%movieHeader
r3dMovieHeader.pixelX=dataProperties.PIXELSIZE_XY;
r3dMovieHeader.pixelY=dataProperties.PIXELSIZE_XY;
r3dMovieHeader.pixelZ=dataProperties.PIXELSIZE_Z;
r3dMovieHeader.lensID=dataProperties.LENSID;
r3dMovieHeader.numCols=stackSize(1);
r3dMovieHeader.numRows=stackSize(2);
r3dMovieHeader.numZSlices=stackSize(3);
r3dMovieHeader.numTimepoints=movieLength;
r3dMovieHeader.numWvs=1;
r3dMovieHeader.wvl=dataProperties.WVL;
time = dataProperties.frameTime';
r3dMovieHeader.Time=time(:)';
r3dMovieHeader.expTime=dataProperties.expTime;
r3dMovieHeader.ndFilter=dataProperties.NDfilter;
r3dMovieHeader.cropInfo=dataProperties.crop;
r3dMovieHeader.correctInfo=[];
%dataProperties: already assigned
%inputSlist
inputSlist = slist;

%save if selected
oldDir = pwd;
switch options.save
    case ''
        return
    case 'ask'
        saveDir = uigetdir(pwd,'select save directory (will create subDir)');
        if isempty(saveDir)|saveDir==0
            return
        end
        fsIdx = strfind(saveDir,fileSep);
        if ~strcmp(saveDir(fsIdx+1:end),dataProperties.name)
            cd(saveDir);
            mkdir(dataProperties.name);
            cd(dataProperties.name);
        else
            cd(saveDir);
        end
    otherwise
        cd(options.save);
end

writemat([dataProperties.name,'.r3c'],synthMovie);
save('r3dMovieHeader','r3dMovieHeader');
save('tmpDataProperties','dataProperties');

cd(oldDir);
