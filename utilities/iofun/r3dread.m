function [image, filename, r3dMovieHeader]= r3dread(filename,start,nTimes,MOVIERANGE,correctionDataFromVar)
%R3DREAD reads the SOFTWORX *.r3d files into matlab (UNIX byte ordering)
%
% SYNOPSIS [image, filename, header]= r3dread(filename,start,nTimes,MOVIERANGE,correctionDataFromVar)
%
% INPUT (all optional)
%       filename :              string, if the filename is omitted,
%                               a file selection dialog is opened
%       start    :              frame with which to start reading {default: 1}
%       nTimes   :              number of frames to read {default: length(movie)}
%       MOVIERANGE :            possible range of values (256 for 8-bit-movies) {2^12}
%       correctionDataFromVar : struct for background correction OR 0. 
%                                   fields: .info.correctFile: filenames of
%                                                   dark movies
%                                           .info.correctFrames: [firstN,
%                                                   lastN] frames used for
%                                                   correction
%                                           .image: image used for
%                                                   correction
%                               if zero, raw movie is loaded. If not, the
%                               program will return a corrected movie, if
%                               there is either correctionDataFromVar or
%                               correctionDataFromFile. 
%
%
% OUTPUT image : image series =[X,Y,Z,WVL,TIME]
%        fname : (optional) string 
%        r3dMovieHeader: structure containing the information from the movie header

% c: 6/3/00	dT
% @ 7/3/01 : the header contains a lot information that is ignored for the moment

%define constant
DV_MAGIC = -16224;
HEADER_SIZE = 1024;

%test input and assign defaults

%set start to 1 if not exist
if nargin < 2 | isempty(start)
    start=1;
end

%assign nTimes below, after movie has been read an # of tp are known

%assume 12-bit movie
if nargin<4 | isempty(MOVIERANGE)
    MOVIERANGE = 2^12;
end

%if no filename, open file selection dialog
if(nargin==0 | isempty(filename))
    [fname,path]=uigetfile('*.r3d;*_R3D.dv','select movie file');
    if(fname(1)==0)
        image=[];
        fname=[];
        return;
    end;
    cd(path);
    filename=[path fname];
else
    path = [];
end;

%fill in path
if isempty(path)
    % is there a path in the filename?
    fileSepIdx = strfind(filename,filesep);
    if isempty(fileSepIdx)
        % no, we are in the right directory
        path = [pwd,filesep];
    else
        % take path from filename
        path = filename(1:fileSepIdx(end));
    end
end

% add filenameextension if neccessary
if isempty(findstr(filename,'.'))
    filename=[filename,'.r3d'];
end;

[file, message] = fopen(filename,'r','b');
if file == -1
    image=[];
    filename =['file ',filename,' not found.'];
    return;
end

% read full header for output if needed
if nargout > 2
    r3dMovieHeader = readr3dheader(filename);
else
    r3dMovieHeader = [];
end


% read in a number which identifies the correct byte ordering
fseek(file,96,-1);
nDVID = fread(file,1,'short');
if (nDVID~=DV_MAGIC)
    %try little endian
    fclose(file);
    [file, message] = fopen(filename,'r','l');
    fseek(file,96,-1);
    nDVID = fread(file,1,'short');
    %if still not good ->  corrupt file
    if (nDVID~=DV_MAGIC)
        error('file is corrupt..')
    end;
end;

% rewind
fseek(file,0,-1);

%-------------
% read header
%--------------

%first BLOCK
block = fread(file,24,'int32');
%extract interesting information
numCol = block(1);
numRow = block(2);
numImages = block(3);
firstImage = block(24);


%skip values to number of time points
fseek(file,180,-1);
numTimes=fread(file,1,'short'); 

% set to all images if not exist
if nargin < 3 | isempty(nTimes)
    nTimes=numTimes;
elseif nTimes == 0
    image = [];
    return
end

% max nTimes to read
nTimes=min(nTimes,numTimes-start+1);

%skip values to number of wavelengths
fseek(file,196,-1);
numWvl=fread(file,1,'short');
numZ=numImages/(numTimes*numWvl);

%-------------------------------
% check for backgroundCorrection
%-------------------------------

% test input
if nargin < 5 | isempty(correctionDataFromVar)
    checkForBGFile = 1;
    correctionDataFromVar = [];
    correctBG = 0; % we don't know yet whether to correct
else
    % if struct: still check for BG, if not, don't correct BG
    if isstruct(correctionDataFromVar)
        checkForBGFile = 1;
        correctBG = 1; % we will check whether there is a file or not
    else % if not zero: bad user, plz read help
        checkForBGFile = 0;
        correctBG = 0;
    end
end

% look for file
if checkForBGFile
    if exist([path,'correctionData.mat'],'file')
        loadStruct = load([path,'correctionData.mat']);
        correctionDataFromFile = loadStruct.correctionData;
        correctBG = 1;
    else
        % there is no file
        correctionDataFromFile = [];
        % leave correctBG as is
    end
end

if correctBG
    %calculate background correction image
    [backgroundCorrectionImage,minStart,maxEnd] = r3dreadCorrectBackground(filename,numRow,numCol,numZ,numWvl,numTimes,correctionDataFromFile,correctionDataFromVar,MOVIERANGE);
    
    %adjust start and nTimes if necessary
    start = start + minStart - 1;
    nTimes=min(nTimes,maxEnd-start+1);
else
    backgroundCorrectionImage = zeros(numRow,numCol,numZ,numWvl);
end

%------------
% load movie
%------------



%allocate mem
image=zeros(numRow,numCol,numZ,numWvl,nTimes);
readIm=zeros(numCol,numRow,numZ);
fseek(file,HEADER_SIZE+firstImage+(start-1)*numCol*numRow*numZ*2,-1);
for t=1:nTimes
    for w=1:numWvl
        readIm(:)=fread(file,numCol*numRow*numZ,'int16');
        %and rotate and normalize
        for z=1:size(readIm,3)
            % do backgroundCorrection in loop because we 
            image(:,:,z,w,t) = readIm(:,:,z)';
        end;
        % adjust image range and correct background
        image(:,:,:,w,t) = (image(:,:,:,w,t)/MOVIERANGE)-backgroundCorrectionImage(:,:,:,w);
    end;
end;

fclose(file);



% header structure:

%struct dv_head {
% BLOCK1
%	long   numCol,numRow,numImages;			   /* nsec +AD0- nz-nw+ACo-nt */
%	long   mode;
%	long   nxst, nyst, nzst;
%	long   mx, my, mz;
%	float xlen, ylen, zlen;
%	float alpha, beta, gamma;
%	long   mapc, mapr, maps;
%	float min1, max1, amean;
%	long   ispg, next;
% END BLOCK1                   offset 96
%	short nDVID,nblank;			 /* nblank preserves byte boundary */
%	char  ibyte[28];
% BLOCK2
%	short nint,nreal;
%	short nres,nzfact;
% END BLOCK2
% BLOCK3
%	float min2,max2,min3,max3,min4,max4;
% END BLOCK3
% BLOCK4
%	short filetype, lens, n1, n2, v1, v2;
% END BLOCK4
% BLOCK5
%	float min5,max5;
% END BLOCK5
% BLOCK6                      offset 180
%	short numtimes;
%	short imagesequence;
% END BLOCK6
% BLOCK7
%	float tiltx, tilty, tiltz;
% END BLOCK7
% BLOCK8                      offset 196
%	short NumWaves, iwav1, iwav2, iwav3, iwav4, iwav5;
% END BLOCK8
% BLOCK9
%	float zorig, xorig, yorig;
% END BLOCK9
%	long   nlab;
%	char  label[800];
%};

