% Reads a TEM movie file or image sequence, computes shift of image vs first frame,
%corrects by wrapping or padding, saves as uncompressed avi
% Josh Sugar and Dave Robinson, Sandia National Labs
%Copyright 2014 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
%there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
%Export of this program may require a license from the United States Government.

%Please cite the Microscopy Today article _blank_ if you use this script.
%***********************************************************************************************
InputForm = questdlg('What is your input file format?','Input format','PNG','TIF','Movie','Movie');

inmovie=0;
if strncmp(InputForm,'PNG',3);
    filetype = 'png';
elseif strncmp(InputForm,'TIF',3);
    filetype = 'tif';
else
    filetype = 'Movie';
    inmovie=1;
end

% Image sequence can be generated from quicktime pro.
% TIFF, when uncompressed, runs faster in both quicktime and matlab.
% TIFF Packbits compression in quicktime doesn't do much for TEM images.
% PNG uses lossless compression, so is more efficient with disk space,
% especially useful if you run quicktime pro on another computer and
% transfer files via USB drive.
% Quicktime tends to choke when reassembling long image sequences into a movie.
% This version generates uncompressed AVI output.

% For movies:
% Containers supported: .avi, .mpg, .mj2 (extremely obscure)
% Windows only: .wmv, .asf, .asx
% Mac only: .mp4, .mov
% Codecs supported: good luck - generally only obsolete ones.
% Microsoft DirectShow codecs include:
% Cinepak, DV, H.264 (MS MP2), ISO MP4 1.0, MS MP4 v3, MJPEG, MP1, MP2

% Bug: sometimes the shifts are off by one width or height, resulting in
% excessive padding in 'pad' mode (innocuous in 'wrap' mode).

% x domain shift f(x-a) transforms to F(w)exp(-iaw)
% Find a, given f(x) and f(x-a)
% F(w)*conj(F(w)exp(-iaw)) gives |F(w)|^2*exp(iaw) (keeping shift but discarding phase info in image)
% ifft of |F(w)| is maximum at zero because all components are in phase there.
% ifft of the whole thing is shifted by a, so identify the location of the maximum.

%% Step 1: Calculate shift x,y pairs for each frame

% Open File
if inmovie;
[moviefilename,inputpath]=uigetfile('*.*','Choose a movie file to drift correct (cancel skips to next step)');

else % image sequence
inputpath=uigetdir('*.*','Choose a folder containing the image sequence (cancel skips to next step)');
olddir=pwd;
cd(inputpath);
wildcard=strcat('*.',filetype);
files=dir(wildcard);
moviefilename = files(1).name;
cd(olddir)
end

if inputpath~=0; % skip to step 2 if canceled

 tic; % start stopwatch

 moviefile=fullfile(inputpath,moviefilename);

 % Set up input file
 if inmovie;
  readobject=VideoReader(moviefile);
  get(readobject);
  % Get number of frames
  if isempty(readobject.NumberofFrames)==1;
    disp('I need to read through this variable bitrate movie to count the frames.')
    read(readobject,inf);
  end
  nFrames = readobject.NumberOfFrames;

 else % image sequence
  files=dir(fullfile(inputpath,strcat('*.',filetype)));
  nFrames = length(files);
 end % inmovie

 shifty = zeros(nFrames,1);
 shiftx = zeros(nFrames,1);

 % Get reference frame (first one in the sequence)

 if inmovie;
  frameref=read(readobject,1);
 else % image sequence
  frameref=imread(fullfile(inputpath,files(1).name),filetype);
 end %inmovie

 imageref=frameref(:,:,1);
 fft_ref=fft2(imageref); 

 [vidHeight vidWidth blank] = size(frameref); % The blank variable here gets rid of extra padding !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 centery=(vidHeight/2)+1;
 centerx=(vidWidth/2)+1;

 wb=waitbar(0,'Please Wait... Calculating Drift...');
 for i=1:nFrames;
    if inmovie;
     imagei=read(readobject,i);
    else % image sequence
     imagei=imread(fullfile(inputpath,files(i).name),filetype);
    end % inmovie
    framei=imagei(:,:,1);
    fft_frame=fft2(framei);
    prod=fft_ref.*conj(fft_frame);
    cc=ifft2(prod);
    [maxy,maxx]=find(fftshift(cc)==max(max(cc)));
    % fftshift moves corners to center; max(max()) gives largest element; find returns indices of that point
 
    shifty(i)=maxy(1)-centery;
    shiftx(i)=maxx(1)-centerx;
    % previous version didn't subtract center point here
    if i>1 % Checks to see if there is an ambiguity problem with FFT because of the periodic boundary in FFT
        if abs(shifty(i)-shifty(i-1))>vidHeight/2;
            shifty(i)=shifty(i)-sign(shifty(i)-shifty(i-1))*vidHeight;
        end
        if abs(shiftx(i)-shiftx(i-1))>vidWidth/2;
            shiftx(i)=shiftx(i)-sign(shiftx(i)-shiftx(i-1))*vidWidth;
        end
    end

    if mod(i,500)==0;
        disp(sprintf('%d frames of %d done',i,nFrames))
        waitbar(i/nFrames,wb);
    end
 end % i loop
 close(wb);
 % save results to file
 shiftdatafilename=fullfile(inputpath,strcat(moviefilename,'_driftdata.mat'));
 if inmovie;
  inputpath = moviefile;
 end
 save(shiftdatafilename,'nFrames','shifty','shiftx','inputpath','filetype','moviefilename');

 disp('I finished calculating the drift data.  It was saved to:')
 disp(shiftdatafilename)
 toc % report time elapsed
 load gong
 player=audioplayer(y,Fs);
 play(player);
end % if didn't cancel drift correct

%% Step 2: apply shift to movie frames

clear all
% Get parameters and input file from user

shiftmode = questdlg('Pad frame with empty space, or wrap within same bounds?','Pad or Wrap','Pad','Wrap','Wrap');

avgframes = str2num(questdlg('Average together how many frames?','Averaging','5','3','1','1'));

bigIter = str2num(questdlg('Max number of input frames per output file','Output file size','15000','5000','2000','2000'));

[shiftfilename,shiftpathname]=uigetfile('*.mat','Please locate the drift data for the movie to correct...');

if shiftfilename~=0; % skip to end if canceled

 tic % start stopwatch

 shiftfile=fullfile(shiftpathname,shiftfilename);
 load(shiftfile);
 inmovie = 0;
 if strncmp(filetype,'Mov',3);
  inmovie = 1;
 end

% Read folder to get basic info

 if inmovie;
  readobject=VideoReader(inputpath);
  get(readobject);
  vidHeight = readobject.Height;
  vidWidth = readobject.Width;
  framerate = readobject.Framerate;

 else % image sequence
  files=dir(fullfile(inputpath,strcat('*.',filetype)));
  nFrames = length(files);
  framerate = 29.97;

  % get height and width from first image
  frameref=imread(fullfile(inputpath,files(1).name),filetype);
  frameref=uint8(frameref(:,:,1));
  [vidHeight vidWidth] = size(frameref); % The blank variable here gets rid of extra padding !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end %inmovie

 if strncmp(shiftmode,'Pad',3);
     newsizey=2*max(abs(shifty))+vidHeight;
	newsizex=2*max(abs(shiftx))+vidWidth;
	% assumes max positive shift = max negative shift; centers reference frame
	
	midindexy=(newsizey-vidHeight)/2+1;
	midindexx=(newsizex-vidWidth)/2+1;
 else
	% wrap
	newsizey=vidHeight;
	newsizex=vidWidth;
	midindexy=1; % these won't be used
	midindexx=1;
 end

 for j=1:floor(nFrames/bigIter)+1;
    disp(sprintf('\nMovie file %d\n',j))
 if floor(nFrames/bigIter)+1==1;
    k=1;
    bigIter=nFrames;
 end
 if j==1;
 k=1;
 elseif j==floor(nFrames/bigIter)+1;
    k=(j-1)*bigIter+1;
    bigIter=mod(nFrames,2000);
 else
    k=(j-1)*bigIter+1;
 end

 % Create array of 8-bit grayscale frame structures
 clear mov_shift
 numofframes = floor(bigIter/avgframes)-1;
 mov_shift(1:numofframes)=struct('cdata',zeros(newsizey,newsizex,1,'uint8'),...
    'colormap',[linspace(0,1,256)',linspace(0,1,256)',linspace(0,1,256)']);
 disp('I allocated the new movie structure...')

 % Process each frame
 for i = 1:numofframes;
    frame_sum=zeros(newsizey,newsizex,'uint16');
    for z = 1:avgframes;
     zz=k+(i-1)*avgframes+z-1;
     if zz>nFrames;
         continue;
     end
     if inmovie;
      imagei=read(readobject,zz);
     else
      imagei=imread(fullfile(inputpath,files(zz).name),filetype);
     end
     framei=uint8(imagei(:,:,1));
	 if strncmp(shiftmode,'Pad',3);
        frame_shift=zeros(newsizey,newsizex,'uint8');
	    frame_shift(midindexy+shifty(zz):midindexy+shifty(zz)+(vidHeight-1),...
			midindexx+shiftx(zz):midindexx+shiftx(zz)+(vidWidth-1))=framei;
		frame_sum=frame_sum+uint16(frame_shift);
	 else
		% Wrap
		frame_sum=frame_sum+uint16(circshift(framei,[shifty(zz) shiftx(zz)]));
     end % if pad or wrap
    end % z loop
    mov_shift(i).cdata=uint8(round(frame_sum/z));
    if mod(i,500)==0;
       disp(sprintf('%d frames of %d done',i,numofframes))
    end
 end % i loop
 k=k+bigIter;

 % Write Movie to new AVI File
 disp('I am writing the movie now...')
 if inmovie;
  [moviepath,moviename,movieext]=fileparts(inputpath);
  filewritename=strcat(fullfile(moviepath,moviename),'corrected_',num2str(j),'.avi');
 else % image sequence
  filewritename=strcat(fullfile(inputpath,moviefilename),'corrected_',num2str(j),'.avi');
 end % inmovie
 movie2avi(mov_shift,filewritename,'fps',framerate,'compression','None')

 end % j loop (bigIter blocks) 

 disp('The corrected movies are complete ...enjoy...')
 toc % report time elapsed
 load gong
 player=audioplayer(y,Fs);
 play(player);
end % if didn't cancel dialog box

%NOTICE:
%For five (5) years from 01/24/2014, the United States Government is granted
%for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
%worldwide license in this data to reproduce, prepare derivative works, and perform
%publicly and display publicly, by or on behalf of the Government. There is provision
%for the possible extension of the term of this license. Subsequent to that period or
%any extension granted, the United States Government is granted for itself and others
%acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
%data to reproduce, prepare derivative works, distribute copies to the public, perform
%publicly and display publicly, and to permit others to do so. The specific term of the
%license can be identified by inquiry made to Sandia Corporation or DOE.
 
%NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY,
%NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR
%IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS
%OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE
%WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
 
%Any licensee of this software has the obligation and responsibility to abide by the
%applicable export control laws, regulations, and general prohibitions relating to the
%export of technical data. Failure to obtain an export control license or other authority
%from the Government may result in criminal liability under U.S. laws.
 
%(End of Notice)