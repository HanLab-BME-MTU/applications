function [synthMovie, r3dMovieHeader, dataProperties, inputSlist, background] = generateSynthMovie(slist,dataProperties,SNR,options,synthMovie)
%generateSyntMovie generates a synthetic yeast movie from coordinates
%
%SYNOPSIS [synthMovie, movieHeader, dataProperties, slist, idlist] = generateSynthMovie(slist,size,SNR,options,synthMovie)
%
%INPUT    slist (opt): real or synthetic slist containing at least the fields slist(t).sp(i).cord
%                       [1x3,pix] and slist(t).sp(i).amp (opt) [betw. 0 and 1 in steps of 2e-12] 
%                       default: see code
%         dataProperties (opt): pixelSize: default 0.0515,0.2
%                               movieSize: default 64,64,18,.,10
%                               frameTime: default 1 stack/sec, 0.03 sec/slice
%                               newName     : will become name, default ['synthMovie-',nowString]
%         SNR (opt):  SNR (signal/sigmaNoise) {5}
%         options (opt): .movie   1: generate movie, 2: add noise, {3}: generateMovie&addNoise
%                        .psf     1: gauss, {2}: theoretical psf 
%                        .microscope [wavelength, NA] default: [525, 1.4]
%                        .background [value] default: 0. if -1, the min of the movie won't be 0, but
%                                    (probably) negative
%                        .saveDir   directory into which the movie should be saved. [directory/'' or 0 (no save)/{'ask'}]
%                        .spots      .nsp number of spots {2} 
%                                    .intRatio intensity ratios (relative to first spot, descending) {[1,0.8]}
%
%         TO BE IMPLEMENTED: synthMovie (opt): synthetic movie to which to add noise
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


%-------------------test input and assign defaults-------------------

%dataProperties (do this first because we need movieSize for slist)
if nargin<2|isempty(dataProperties)
    dataProperties = [];
end
%now fill all fields that haven't been filled yet
if ~isfield(dataProperties,'PIXELSIZE_XY')
    dataProperties.PIXELSIZE_XY = 0.0515;
end
if ~isfield(dataProperties,'PIXELSIZE_Z')
    dataProperties.PIXELSIZE_Z = 0.2;
end
if ~isfield(dataProperties,'movieSize')
    dataProperties.movieSize = [64,64,18,10];
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
    dataProperties.FT_SIGMA=[1.5289 1.5289 1.1756];
    dataProperties.FILTERPRM=[1.5289 1.5289 1.1756 7 7 5];
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
    if length(options.microscope)==2
        wvl = options.microscope(1);
        NA = options.microscope(2);
        if wvl>1
            error('wavelength in microns!')
        else
            dataProperties.WVL = wvl;
        end
        if NA<1|NA>1.5
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
end    

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


%slist
%place tags
if nargin<1 | isempty(slist)
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
while isempty(slist(tstart))&tstart<=length(slist)
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

%make sure that we do not have too large amplitudes
ampList = [];
for t = 1:movieLength
    ampList = [ampList;cat(2,slist(t).sp.amp)];
end
if any(ampList>1)
    error('not normalized amplitudes!')
else
    DYNAMICRANGE = max(ampList(:));
end

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

%psf size
psfSze=[17 17 11];
hPsf=floor(psfSze/2);

%init psf
if options.psf == 2
    disp('sigmas of real psf could be a bit off')
    magnification = 100;
    % from fitparmspos.mat @ c:\usr\thomann\matlab\dtOwn\wideFieldPsf\psfDataDir\psfdri1003
    psfparms=[ 0.9926    1.5032   25.2041   41.9328    0.0778    0.1049   -0.1882];
    psfob=psf(NA,magnification,wvl*1000,10000,1,psfparms(1),psfparms(2));
end

% init movie. tmp: too large movie, that we cut back later
synthStackTmp = zeros([stackSize+2*hPsf]);

%loop time points
for t=1:movieLength
    if any(jobList==1)
        for i=1:nsp
            pos=slist(t).sp(i).cord;
            %chage to matlab coords an:
            pos=[pos(2) pos(1) pos(3)]+hPsf;
            %we can only assign values to full pixels
            posfullpix = floor(pos);
            subpixel_shift=pos-posfullpix;
            
            switch options.psf
                case 1 %gaussian psf
                    psfdata=multiGaussFit(psfSze,[subpixel_shift 1 0],dataProperties);
                    %make sure that our maximum == 1*slist(t).sp(i).amp
                    psfdata=slist(t).sp(i).amp*psfdata/max(psfdata(:));
                case 2 %calculated psf
                    psfdata=makePSF(psfob,psfSze,[1000*dataProperties.PIXELSIZE_XY,1000*dataProperties.PIXELSIZE_Z],[subpixel_shift 0 0]);
                    %make sure that our maximum == 1*slist(t).sp(i).amp
                    psfdata=slist(t).sp(i).amp*psfdata/max(psfdata(:));
            end
            
            %place psf in tmp block     
            synthStackTmp(posfullpix(1)-hPsf(1):posfullpix(1)+hPsf(1), posfullpix(2)-hPsf(2):posfullpix(2)+hPsf(2),posfullpix(3)-hPsf(3):posfullpix(3)+hPsf(3))=...
                synthStackTmp(posfullpix(1)-hPsf(1):posfullpix(1)+hPsf(1), posfullpix(2)-hPsf(2):posfullpix(2)+hPsf(2),posfullpix(3)-hPsf(3):posfullpix(3)+hPsf(3))+psfdata;
            
        
        
        end;
        
        %cut tmp block
        synthMovie(:,:,:,1,t)=...
            synthStackTmp(hPsf(1)+1:end-hPsf(1),hPsf(2)+1:end-hPsf(2),hPsf(3)+1:end-hPsf(3));
        %empty tmp block
        synthStackTmp = zeros([stackSize+2*hPsf]);
        
    end
    
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
    